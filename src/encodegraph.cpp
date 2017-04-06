/*
    daccord
    Copyright (C) 2017 German Tischler

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
#include <libmaus2/dazzler/align/GraphDecoder.hpp>
#include <libmaus2/dazzler/align/GraphEncoder.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <config.h>

static void checkGraph(std::istream & istr, std::vector<std::string> const & arg)
{
	libmaus2::dazzler::align::GraphDecoder GD(istr);

	libmaus2::dazzler::align::GraphDecoderContext::shared_ptr_type scontext = GD.getContext();
	libmaus2::dazzler::align::GraphDecoderContext & context = *scontext;

	for ( uint64_t i = 0; i < arg.size(); ++i )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type AF(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(arg[i]));
		libmaus2::dazzler::align::Overlap OVL;

		while ( AF->peekNextOverlap(OVL) )
		{
			int64_t const aid = OVL.aread;
			std::vector < libmaus2::dazzler::align::Overlap > V;

			while ( AF->peekNextOverlap(OVL) && OVL.aread == aid )
			{
				AF->getNextOverlap(OVL);
				V.push_back(OVL);
			}

			GD.decode(istr,aid,context);

			assert ( context.size() == V.size() );

			for ( uint64_t j = 0; j < V.size(); ++j )
			{
				assert ( V[j].isInverse() == context[j].inv );
				assert ( V[j].aread == context[j].aread );
				assert ( V[j].bread == context[j].bread );
				assert ( V[j].path.abpos == context[j].abpos );
				assert ( V[j].path.aepos == context[j].aepos );
				assert ( V[j].path.bbpos == context[j].bbpos );
				assert ( V[j].path.bepos == context[j].bepos );
				assert ( V[j].path.diffs == context[j].diffs );
			}
		}
	}

	GD.returnContext(scontext);

}

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

static int getDefaultVerbose()
{
	return 0;
}

static int getDefaultCheck()
{
	return 0;
}

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

static std::string helpMessage(libmaus2::util::ArgParser const & arg)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("verbose", formatRHS("verbosity",getDefaultVerbose())));
	optionMap . push_back ( std::pair < std::string, std::string >("check", formatRHS("check graph after creating it",getDefaultCheck())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));

	uint64_t maxlhs = 0;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		assert ( ita->first.size() );

		if ( ita->first.size() == 1 )
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+1));
		else
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+2));
	}

	std::ostringstream messtr;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		std::string const key = ita->first;

		messtr << "\t";
		messtr << std::setw(maxlhs) << std::setfill(' ');
		if ( key.size() == 1 )
			messtr << (std::string("-")+key);
		else
			messtr << (std::string("--")+key);

		messtr << std::setw(0);

		messtr << ": ";

		messtr << ita->second;
		messtr << "\n";
	}

	return messtr.str();
}


int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 1 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.graph in.las ...\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::timing::RealTimeClock rtc;
			rtc.start();

			std::string const out = arg[0];
			std::vector<std::string> Vin;
			for ( uint64_t i = 1; i < arg.size(); ++i )
				Vin.push_back(arg[i]);

			int const verbose = arg.uniqueArgPresent("verbose");
			libmaus2::dazzler::align::GraphEncoder::encodegraph(out,Vin,getTmpFileBase(arg),verbose);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			int const check = arg.uniqueArgPresent("check");

			if ( check )
			{
				std::cerr << "[V] checking graph...";
				libmaus2::aio::InputStreamInstance istr(out);
				checkGraph(istr,Vin);
				std::cerr << "done." << std::endl;
			}

			return EXIT_SUCCESS;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
