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
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/bambam/BamAccessor.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/TrueOverlap.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>
#include <libmaus2/dazzler/align/DalignerIndexDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <config.h>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

/*
 parameters:

 -t : default number of logical cores, threads
 */

static std::string helpMessage(libmaus2::util::ArgParser const & arg)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("verbose", formatRHS("verbosity",false)));
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

struct BamAlignmentTypeInfo
{
	typedef libmaus2::bambam::BamAlignment element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		return pointer_type();
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct BamAlignmentAllocator
{
	typedef libmaus2::bambam::BamAlignment element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct BamAlignmentContainer
{
	typedef BamAlignmentContainer this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string const bamfn;
	libmaus2::bambam::BamNumericalIndexDecoder indexdec;
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::bambam::BamAlignment,BamAlignmentAllocator,BamAlignmentTypeInfo> & bamAlignmentFreeList;
	std::map < uint64_t, libmaus2::bambam::BamAlignment::shared_ptr_type > M;

	BamAlignmentContainer(
		std::string const & rbamfn,
		libmaus2::parallel::LockedGrowingFreeList<libmaus2::bambam::BamAlignment,BamAlignmentAllocator,BamAlignmentTypeInfo> & rbamAlignmentFreeList
	) : bamfn(rbamfn), indexdec(libmaus2::bambam::BamNumericalIndexBase::getIndexName(bamfn)), bamAlignmentFreeList(rbamAlignmentFreeList)
	{

	}

	void clear()
	{
		for ( std::map < uint64_t, libmaus2::bambam::BamAlignment::shared_ptr_type >::iterator it = M.begin(); it != M.end(); ++it )
			bamAlignmentFreeList.put(it->second);
		M.clear();
	}

	libmaus2::bambam::BamAlignment const & operator[](uint64_t const i)
	{
		if ( M.find(i) == M.end() )
		{
			libmaus2::bambam::BamAccessor A(bamfn,indexdec,i);
			libmaus2::bambam::BamAlignment::shared_ptr_type P = bamAlignmentFreeList.get();
			libmaus2::bambam::BamAlignment & B = A[i];
			P->swap(B);
			M[i] = P;
		}

		assert ( M.find(i) != M.end() );

		return *(M.find(i)->second);
	}
};

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

int marktrue(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const db = arg[1];
	std::string const bam = arg[2];
	std::vector<std::string> VI(arg.restargs.begin()+3,arg.restargs.end());

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

	libmaus2::dazzler::db::DatabaseFile DB(db);
	DB.computeTrimVector();

	libmaus2::parallel::LockedGrowingFreeList<libmaus2::bambam::BamAlignment,BamAlignmentAllocator,BamAlignmentTypeInfo> bamAlignmentFreeList;
	libmaus2::autoarray::AutoArray < BamAlignmentContainer::unique_ptr_type > ABAC(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::TrueOverlap::unique_ptr_type > ATO(numthreads);
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(VI);
	libmaus2::dazzler::align::TrueOverlapStats TOS;

	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		BamAlignmentContainer::unique_ptr_type t(new BamAlignmentContainer(bam,bamAlignmentFreeList));
		ABAC[i] = UNIQUE_PTR_MOVE(t);

		libmaus2::dazzler::align::TrueOverlap::unique_ptr_type o(new libmaus2::dazzler::align::TrueOverlap(TOS,DB,tspace));
		ATO[i] = UNIQUE_PTR_MOVE(o);
	}

	std::string const tmpprefix = getTmpFileBase(arg);
	libmaus2::dazzler::align::AlignmentWriterArray AWA(tmpprefix + "AWA",numthreads,tspace);

	for ( uint64_t i = 0; i < VI.size(); ++i )
	{
		std::string const & las = VI[i];
		std::string const indalignerlasindexname = libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(las);

		if (
			! libmaus2::util::GetFileSize::fileExists(indalignerlasindexname)
			||
			libmaus2::util::GetFileSize::isOlder(indalignerlasindexname,las)
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(las,&std::cerr);
		}

		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(numthreads);

		for ( uint64_t i = 0; i < Adalindex.size(); ++i )
		{
			libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
				new libmaus2::dazzler::align::DalignerIndexDecoder(las,indalignerlasindexname)
			);
			Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
		}


		int64_t const minaread = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(las);
		int64_t const maxaread = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(las);

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for ( int64_t j = minaread; j <= maxaread; ++j )
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			libmaus2::dazzler::align::AlignmentFileDecoder::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileAt(las,j,j+1,*Adalindex[tid]));

			BamAlignmentContainer & BAC = *ABAC[tid];
			BAC.clear();
			libmaus2::dazzler::align::TrueOverlap & TO = *ATO[tid];

			libmaus2::bambam::BamAlignment const & bam_a = BAC[j];
			libmaus2::dazzler::align::AlignmentWriter & AW = AWA[tid];

			libmaus2::dazzler::align::Overlap OVL;
			while ( pdec->getNextOverlap(OVL) )
			{
				libmaus2::bambam::BamAlignment const & bam_b = BAC[OVL.bread];
				bool const istrue = TO.trueOverlap(OVL,bam_a,bam_b);

				if ( istrue )
					OVL.flags |= OVL.getTrueFlag();
				else
					OVL.flags &= (~static_cast<uint64_t>(OVL.getTrueFlag()));

				AW.put(OVL);
			}

			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] " << j << std::endl;
			}
		}
	}

	AWA.merge(outfn,tmpprefix+"_AWA_merge");

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las in.db in.bam in.las ...\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return marktrue(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
