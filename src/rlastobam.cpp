/*
    daccord
    Copyright (C) 2018 German Tischler-Höhle

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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/lcs/NPLinMem.hpp>
#include <libmaus2/bambam/CigarStringParser.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamBlockWriterBase.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/parallel/NumCpus.hpp>

void process(
	std::vector < libmaus2::dazzler::align::Overlap > & VOVL,
	libmaus2::dazzler::db::DatabaseFile & DBa,
	libmaus2::dazzler::db::DatabaseFile & DBb,
	int64_t const tspace,
	libmaus2::bambam::BamBlockWriterBase & writer,
	uint64_t const numthreads
)
{
	std::vector < ::libmaus2::fastx::UCharBuffer::shared_ptr_type > VB(VOVL.size());

	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
	#endif
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
	{
		libmaus2::dazzler::align::Overlap const & OVL = VOVL[i];

		libmaus2::lcs::NPLinMem np;
		libmaus2::lcs::AlignmentTraceContainer ATC;
		libmaus2::autoarray::AutoArray< std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
                libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Aop;
                ::libmaus2::fastx::UCharBuffer::shared_ptr_type sbuffer(new ::libmaus2::fastx::UCharBuffer);
                ::libmaus2::fastx::UCharBuffer & buffer = *sbuffer;
                VB[i] = sbuffer;
                libmaus2::bambam::BamSeqEncodeTable seqenc;

		std::string const a = DBa.decodeRead(OVL.aread,false);
		std::string const b = DBb.decodeRead(OVL.bread,OVL.isInverse());
		std::string const s = DBb.getReadName(OVL.bread);
		std::string const h = std::string(b.size(),'H');

		OVL.computeTrace(
			reinterpret_cast<uint8_t const *>(a.c_str()),
			reinterpret_cast<uint8_t const *>(b.c_str()),
			tspace,
			ATC,
			np
		);

		uint64_t const ciglen = libmaus2::bambam::CigarStringParser::traceToCigar(ATC,Aopblocks,Aop,0,OVL.path.bbpos,b.size() - OVL.path.bepos,0);

		libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment(
			buffer,seqenc,
			s.c_str(),
			s.size(),
			OVL.aread,
			OVL.path.abpos,
			255 /* mapq */,
			OVL.isInverse() ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0,
			Aop.begin(),ciglen,
			-1 /* nextrefid */,
			-1 /* nextpos */,
			0 /* template length */,
			b.begin(),
			b.size(),
			h.begin()
		);

		libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(buffer,"RG","RGID");
		libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(buffer,"tr",'i',OVL.isTrue());
	}

	for ( uint64_t i = 0; i < VOVL.size(); ++i )
	{
		::libmaus2::fastx::UCharBuffer & buffer = *(VB[i]);
		writer.writeBamBlock(buffer.buffer,buffer.length);
	}

	VOVL.resize(0);
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		libmaus2::util::ArgInfo arginfo(argc,argv);
		
		if ( arg.size() < 4 )
		{
			std::cerr << "usage: " << argv[0] << " [<-tnumthreads>] <a.db> <b.db> <ref.fasta> <in.las> >out.bam";
			return EXIT_FAILURE;
		}
		
		std::string const dba = arg[0];
		std::string const dbb = arg[1];
		std::string const ref = arg[2];
		std::string const las = arg[3];

		uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::parallel::NumCpus::getNumLogicalProcessors();

		arginfo.replaceKey("level","0");

		std::map<uint64_t,std::string> M;
		{
			libmaus2::fastx::FastAReader R(ref);
			libmaus2::fastx::FastAReader::pattern_type pattern;
			for ( uint64_t i = 0; R.getNextPatternUnlocked(pattern); ++i )
				M [ i ] = pattern.getShortStringId();
		}

		libmaus2::dazzler::db::DatabaseFile DBa(dba);
		DBa.computeTrimVector();
		libmaus2::dazzler::db::DatabaseFile DBb(dbb);
		DBb.computeTrimVector();

		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type tptr(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(las));
		libmaus2::dazzler::align::Overlap OVL;
		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(las);

		std::ostringstream headerstream;
		headerstream << "@HD\tVN:1.5\tSO:unknown\n";
		for ( uint64_t i = 0; i < DBa.size(); ++i )
		{
			uint64_t const l = DBa[i].size();
			std::string const s = M.find(i)->second;
			headerstream << "@SQ\tSN:" << s << "\tLN:" << l << "\n";
		}
		headerstream << "@RG\tID:RGID\tSM:SAMPLE\n";

		libmaus2::bambam::BamHeader bamheader(headerstream.str());
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(libmaus2::bambam::BamBlockWriterBaseFactory::construct(bamheader, arginfo));

		std::vector < libmaus2::dazzler::align::Overlap > VOVL;

		while ( tptr->getNextOverlap(OVL) )
		{
			VOVL.push_back(OVL);

			if ( VOVL.size() >= 1024 )
			{
				process(VOVL,DBa,DBb,tspace,*writer,numthreads);
			}

		}

		process(VOVL,DBa,DBb,tspace,*writer,numthreads);

		writer.reset();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
