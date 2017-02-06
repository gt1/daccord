/*
    daccord
    Copyright (C) 2016-2017 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <cstdlib>
#include <libmaus2/util/LineBuffer.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/bambam/BamFlagBase.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamAlignmentDecoderBase.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

static bool startsWith(std::string const & line, std::string const & prefix)
{
	return line.size() >= prefix.size() &&
		line.substr(0,prefix.size()) == prefix;
}

static std::vector<libmaus2::fastx::FastAReader::pattern_type> readText(std::string const & textfn)
{
	libmaus2::fastx::FastAReader FA(textfn);
	std::vector<libmaus2::fastx::FastAReader::pattern_type> V;
	libmaus2::fastx::FastAReader::pattern_type pat;
	while ( FA.getNextPatternUnlocked(pat) )
	{
		V.push_back(pat);
	}

	return V;
}

static std::vector<std::string> split(std::string const & line)
{
	std::vector<std::string> tok;
	uint64_t l = 0;

	while ( l < line.size() )
	{
		while ( l < line.size() && ::std::isspace(line[l]) )
			++l;

		uint64_t h = l;
		while ( h < line.size() && !::std::isspace(line[h]) )
			++h;

		if ( l != h )
			tok.push_back(line.substr(l,h-l));

		l = h;
	}

	return tok;
}

struct MatchLine
{
	std::string name;
	int64_t start;
	uint64_t size;
	bool strand;
	uint64_t seqsize;
	std::string text;

	MatchLine()
	{}
	MatchLine(
		std::string const & rname,
		int64_t const & rstart,
		uint64_t const & rsize,
		bool const & rstrand,
		uint64_t const & rseqsize,
		std::string const & rtext
	) : name(rname), start(rstart), size(rsize), strand(rstrand), seqsize(rseqsize), text(rtext) {}
};

struct Match
{
	std::vector < MatchLine > ML;

	void handle(
		libmaus2::bambam::BamBlockWriterBase & writer,
		std::map < std::string, uint64_t> const & M, std::vector<libmaus2::fastx::FastAReader::pattern_type> const & /* Vref */,
		std::map<std::string,std::string> const & replmap
		)
	{
		if ( ML.size() == 2 )
		{
			for ( uint64_t i = 1; i < ML.size(); ++i )
			{
				bool const ok = ( ML[i].text.size() == ML[0].text.size() );
				if ( ! ok )
				{
					std::cerr << "[E] malformed alignment" << std::endl;
					return;
				}
			}

			if ( !ML[0].strand )
			{
				std::cerr << "[E] cannot handle RC on refseq" << std::endl;
				return;
			}

			std::vector < libmaus2::bambam::BamFlagBase::bam_cigar_ops > Vop;
			std::ostringstream readstr;
			std::ostringstream refstr;
			for ( uint64_t i = 0; i < ML[0].text.size(); ++i )
			{
				if ( ML[1].text[i] != '-' )
					readstr.put(ML[1].text[i]);
				if ( ML[0].text[i] != '-' )
					refstr.put(ML[0].text[i]);

				if ( ML[0].text[i] == ML[1].text[i] )
				{
					assert ( ML[0].text[i] != '-' );
					#if defined(CIGARDEBUG)
					Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH);
					#else
					Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL);
					#endif
				}
				else
				{
					assert ( ML[0].text[i] != ML[1].text[i] );

					if ( ML[0].text[i] == '-' )
						Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS);
					else if ( ML[1].text[i] == '-' )
						Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL);
					else
					{
						#if defined(CIGARDEBUG)
						Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH);
						#else
						Vop.push_back(libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF);
						#endif
					}
				}
			}


			uint64_t l = 0;
			std::ostringstream cigstream;
			if ( ML[1].start )
			{
				cigstream << ML[1].start << 'H';
			}
			while ( l < Vop.size() )
			{
				uint64_t h = l+1;
				while ( h < Vop.size() && Vop[h] == Vop[l] )
					++h;

				cigstream << (h-l) << Vop[l];

				l = h;
			}

			if ( ML[1].start + ML[1].size < ML[1].seqsize )
			{
				cigstream << (ML[1].seqsize-(ML[1].start + ML[1].size)) << 'H';
			}

			std::string const readdata = readstr.str();
			std::string const cigstr = cigstream.str();
			std::string const qual(readdata.size(),255);

			if ( replmap.find(ML[0].name) != replmap.end() )
				ML[0].name = replmap.find(ML[0].name)->second;

			if ( M.find(ML[0].name) == M.end() )
			{
				std::cerr << "[E] unknown ref seq " << ML[0].name << std::endl;
				return;
			}

			uint64_t const refid = M.find(ML[0].name)->second;

			::libmaus2::fastx::UCharBuffer buffer;
			libmaus2::bambam::BamSeqEncodeTable seqenc;

			#if 0
			std::cerr << "refstr=" << refstr.str() << std::endl;
			std::cerr << "compare=" << Vref[refid].spattern.substr(ML[0].start,ML[0].size) << std::endl;
			#endif

			libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment(
				buffer,seqenc,
				ML[1].name,
				refid,
				ML[0].start,
				255 /* map */,
				ML[1].strand ? 0 : libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE,
				cigstr,
				-1,
				-1,
				0,
				readdata,
				qual,
				0,
				true
			);

			::libmaus2::bambam::MdStringComputationContext mdcontext;
			std::string const ref = refstr.str();

			// std::cerr << ref << std::endl;

			libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(buffer.buffer,buffer.length,mdcontext,ref.begin());

			libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(buffer,"MD",mdcontext.md.get());
			libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(buffer,"NM",'i',mdcontext.nm);

			#if 0
			::libmaus2::bambam::BamFormatAuxiliary aux;
			libmaus2::bambam::BamAlignmentDecoderBase::formatAlignment(
				std::cout,
				buffer.buffer,
				buffer.length,
				header,
				aux
			);

			std::cout << std::endl;
			#endif

			writer.writeBamBlock(buffer.buffer,buffer.length);
		}
		else if ( ML.size() )
		{
			std::cerr << "[E] cannot handle multi alignment" << std::endl;
		}
	}
};

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		libmaus2::util::ArgParser const arg(argc,argv);

		std::map<std::string,std::string> replmap;
		if ( arg.size() > 1 )
		{
			std::string replname = arg[1];
			libmaus2::aio::InputStreamInstance ISI(replname);
			while ( ISI )
			{
				std::string line;
				std::getline(ISI,line);
				if ( line.size() && line.find('\t') != std::string::npos )
				{
					uint64_t const p = line.find('\t');
					std::string const key = line.substr(0,p);
					std::string const value = line.substr(p+1);
					replmap[key] = value;
				}
			}
		}

		libmaus2::util::LineBuffer LB(std::cin);

		std::vector<libmaus2::fastx::FastAReader::pattern_type> Vref = readText(arg[0]);
		std::map < std::string, uint64_t> M;
		std::ostringstream headerstream;
		headerstream << "@HD\tVN:1.5\tSO:unknown\n";
		for ( uint64_t i = 0; i < Vref.size(); ++i )
		{
			M [ Vref[i].getShortStringId() ] = i;
			headerstream << "@SQ\tSN:" << Vref[i].getShortStringId() << "\tLN:" << Vref[i].spattern.size() << "\n";
		}

		std::string const headertext = headerstream.str();
		//std::cout << headertext;

		libmaus2::bambam::BamHeader bamheader(headertext);
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(libmaus2::bambam::BamBlockWriterBaseFactory::construct(bamheader, arginfo));

		char const * a = NULL;
		char const * e = NULL;
		Match match;

		while ( LB.getline(&a,&e) )
		{
			// std::cerr << std::string(a,e) << std::endl;
			std::string const line(a,e);

			if ( startsWith(line,"s") || startsWith(line,"a") )
			{
				std::vector<std::string> tokens = split(line);

				if ( tokens.size() >= 7 && tokens[0] == "s" )
				{
					std::string const readname = tokens[1];
					int64_t const start = atol(tokens[2].c_str());
					uint64_t const size = atol(tokens[3].c_str());
					std::string const strand = tokens[4];

					if ( strand != "+" && strand != "-" )
					{
						std::cerr << "[E] unknown strand " << strand << std::endl;
						continue;
					}

					uint64_t const seqsize = atol(tokens[5].c_str());
					std::string const text = tokens[6];

					MatchLine ML(readname,start,size,strand == "+",seqsize,text);
					match.ML.push_back(ML);
				}
				if ( tokens.size() >= 1 && tokens[0] == "a" )
				{
					match.handle(*writer,M,Vref,replmap);
					match.ML.resize(0);
				}
			}
		}

		if ( match.ML.size() )
			match.handle(*writer,M,Vref,replmap);

		writer.reset();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
