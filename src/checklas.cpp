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
#include <libmaus2/bambam/BamAccessor.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/aio/TempFileArray.hpp>
#include <config.h>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/parallel/LockedBool.hpp>
#include <libmaus2/parallel/TerminatableSynchronousQueue.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static uint64_t getDefaultMinLen()
{
	return 1000;
}

struct UPair
{
	uint64_t key;
	uint64_t value;

	UPair() : key(0), value(0) {}
	UPair(uint64_t const rkey, uint64_t const rvalue) : key(rkey), value(rvalue) {}
	UPair(std::istream & in)
	:
		key(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		value(libmaus2::util::NumberSerialisation::deserialiseNumber(in))
	{

	}

	std::istream & deserialise(std::istream & in)
	{
		*this = UPair(in);
		return in;
	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,key);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,value);
		return out;
	}

	bool operator<(UPair const & M) const
	{
		if ( key != M.key )
			return key < M.key;
		else
			return value < M.value;
	}
};


bool isBlock(std::string const & fn, int64_t & blockid)
{
	std::deque<std::string> const tokens = libmaus2::util::stringFunctions::tokenize(fn,std::string("."));

	// std::cerr << "checking " << fn << " tokens " << tokens.size() << std::endl;

	if ( tokens.size() < 3 )
		return false;

	std::istringstream istr(tokens[tokens.size()-2]);
	uint64_t i;
	istr >> i;
	blockid = i;

	return istr && istr.peek() == std::istream::traits_type::eof();
}

void printResultLine(
	int64_t const aread,
	uint64_t const miscnt,
	uint64_t const gotcnt,
	uint64_t const extracnt,
	libmaus2::bambam::BamHeader const & header,
	libmaus2::bambam::BamAlignment const & algn,
	double const mintrue,
	double const avgtrue,
	double const maxtrue,
	uint64_t const dcnt,
	std::string const & comment
	)
{
	double const erate = algn.getErrorRate();
	libmaus2::math::IntegerInterval<int64_t> const I = algn.getReferenceInterval();
	int32_t const refid = algn.getRefID();

	libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::coutlock);
	// if ( mintrue <= 1 && avgtrue < 0.5 )
	std::cout << "[R] " << aread << " missing " << miscnt << " got " << gotcnt << " extra " << extracnt << " erate " << erate << " "
		<< header.getRefIDName(refid) << ":" << I.from << "," << I.to << " " << ((miscnt+gotcnt) ? (static_cast<double>(miscnt)/(miscnt+gotcnt)) : 0)
		<< " " << mintrue << " " << avgtrue << " " << maxtrue << " " << dcnt << (comment.size() ? (std::string(" ") + comment) : comment) << "\n";
}

/*
 * checklas <dbfn> <perfpile.las> <reads.bam> <in.las>
 *
 * -t32: numthreads
 * --minsiglen1000: min signif align len
 * --seqdepth30: sequencing depth (for detecting how many lowest error number alignments are true)
 * --verbose: verbosity
 * --tlow: read id start (inclusive)
 * --thigh: read id end (inclusive)
 * --mark: generate marked file (using flag for marking true overlaps)
 * --mis: generate missed true overlaps file
 * --keep: generate kept true overlaps file
 * -T: tempfilebase
 *
 */

static uint64_t getDefaultMinSigLen()
{
	return 1000;
}

static uint64_t getDefaultSeqDepth()
{
	return 30;
}

int checklas(libmaus2::util::ArgParser const & arg)
{
	std::string const bamfn = arg[2];
	std::string const genfn = arg[1];
	std::string const dbfn = arg[0];
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	int volatile failed = 0;
	libmaus2::parallel::PosixSpinLock failedlock;
	int64_t const minlen = arg.uniqueArgPresent("minlen") ? arg.getUnsignedNumericArg<uint64_t>("minlen") : 1000;

	int64_t const gentspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(genfn);

	libmaus2::bambam::BamHeader::unique_ptr_type Pheader;
	{
		libmaus2::bambam::BamDecoder dec(bamfn);
		libmaus2::bambam::BamHeader::unique_ptr_type Theader(dec.getHeader().uclone());
		Pheader = UNIQUE_PTR_MOVE(Theader);
	}
	libmaus2::bambam::BamHeader const & header = *Pheader;

	uint64_t const genindexmod = 1024;
	libmaus2::bambam::BamNumericalIndexGenerator::indexFileCheck(bamfn,genindexmod,1 /* numthreads */);

	if (
		! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(genfn))
		||
		libmaus2::util::GetFileSize::isOlder(libmaus2::dazzler::align::OverlapIndexer::getIndexName(genfn),genfn)
	)
	{
		libmaus2::dazzler::align::OverlapIndexer::constructIndex(genfn,&std::cerr);
	}


	uint64_t const minsiglen = arg.uniqueArgPresent("minsiglen") ? arg.getUnsignedNumericArg<uint64_t>("minsiglen") : getDefaultMinSigLen();
	uint64_t const seqdepth = arg.uniqueArgPresent("seqdepth") ? arg.getUnsignedNumericArg<uint64_t>("seqdepth") : getDefaultSeqDepth();
	bool const verbose = arg.uniqueArgPresent("verbose");

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();

	int64_t const tlow = arg.uniqueArgPresent("tlow") ? arg.getUnsignedNumericArg<uint64_t>("tlow") : std::numeric_limits<int64_t>::min();
	int64_t const thigh = arg.uniqueArgPresent("thigh") ? arg.getUnsignedNumericArg<uint64_t>("thigh") : std::numeric_limits<int64_t>::max();

	bool const flagmark = arg.uniqueArgPresent("mark");
	bool const flagmis = arg.uniqueArgPresent("mis");
	bool const flagkeep = arg.uniqueArgPresent("keep");

	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	for ( uint64_t i = 3; (!failed) && i < arg.size(); ++i )
	{
		std::string const infn = arg[i];

		if (
			! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(infn))
			||
			libmaus2::util::GetFileSize::isOlder(libmaus2::dazzler::align::OverlapIndexer::getIndexName(infn),infn)
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(infn,&std::cerr);
		}

		int64_t blockid = 0;
		isBlock(infn,blockid);
		std::pair<uint64_t,uint64_t> const BI = DB.getTrimmedBlockInterval(blockid);

		// std::cerr << "BI=" << BI.first << "," << BI.second << std::endl;

		libmaus2::aio::TempFileArray TFA(tmpfilebase+"_tpnum",numthreads);

		int64_t minaread = std::max(tlow, libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(infn));
		int64_t maxaread = std::min(thigh,libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(infn));

		if ( maxaread < minaread )
			continue;

		int64_t const toparead = maxaread + 1;
		int64_t const numreads = toparead - minaread;

		std::cerr << "[V] processing read interval [" << minaread << "," << toparead << ")" << std::endl;

		uint64_t const tnumpacks = 8 * numthreads;
		uint64_t const packsize = (numreads + tnumpacks - 1)/tnumpacks;
		uint64_t const numpacks = (numreads + packsize - 1)/packsize;
		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(infn);

		if ( tspace != gentspace )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "tspace mismatch between gen " << gentspace << " and provided las file " << tspace << std::endl;
			lme.finish();
			throw lme;
		}

		uint64_t maxreq = 0;

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t packid = 0; packid < numpacks; ++packid )
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			std::ostream & out = TFA[tid];

			uint64_t const packlow = minaread + packid * packsize;
			uint64_t const packhigh = std::min(static_cast<int64_t>(packlow+packsize),toparead);
			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(infn,packlow,packhigh));
			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pgen(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(genfn,packlow,packhigh));

			libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;

			int64_t prevaread = std::numeric_limits<int64_t>::min();
			bool prevvalid = false;
			libmaus2::dazzler::align::Overlap OVL;
			libmaus2::dazzler::align::Overlap OVLgen;

			while ( Plas->getNextOverlap(OVL) )
			{
				if ( OVL.aread != prevaread )
				{
					if ( prevvalid )
					{
						while ( Pgen->peekNextOverlap(OVLgen) && OVLgen.aread < OVL.aread )
							Pgen->getNextOverlap(OVLgen);

						uint64_t tsum = 0;
						while ( Pgen->peekNextOverlap(OVLgen) && OVLgen.aread == OVL.aread )
						{
							Pgen->getNextOverlap(OVLgen);
							uint64_t const o = OVL.getTracePoints(tspace, 0 /* trace id */, TPV);
							tsum += o;
						}

						UPair(prevaread,tsum).serialise(out);
						// std::cerr << "U" << std::endl;
						maxreq = std::max(maxreq,tsum);
					}

					prevvalid = true;
					prevaread = OVL.aread;
				}

				// uint64_t const o = OVL.getTracePoints(tspace, 0 /* trace id */, TPV);
			}

			if ( prevvalid )
			{
				while ( Pgen->peekNextOverlap(OVLgen) && OVLgen.aread < OVL.aread )
					Pgen->getNextOverlap(OVLgen);

				uint64_t tsum = 0;
				while ( Pgen->peekNextOverlap(OVLgen) && OVLgen.aread == OVL.aread )
				{
					Pgen->getNextOverlap(OVLgen);
					uint64_t const o = OVL.getTracePoints(tspace, 0 /* trace id */, TPV);
					tsum += o;
				}

				UPair(prevaread,tsum).serialise(out);
				// std::cerr << "U" << std::endl;
				maxreq = std::max(maxreq,tsum);
			}
		}

		if ( failed )
		{
			return EXIT_FAILURE;
		}

		uint64_t const memtotal = 1024*1024*1024;

		uint64_t const memfull = std::max(memtotal,maxreq);
		std::string const tpmergefn = tmpfilebase+"_tpnummerge";
		libmaus2::util::TempFileRemovalContainer::addTempFile(tpmergefn);
		std::less<UPair> upcomp;
		TFA.merge<UPair>(tpmergefn,upcomp);

		{
			libmaus2::aio::SerialisedPeeker<UPair>::unique_ptr_type UPP(new libmaus2::aio::SerialisedPeeker<UPair>(tpmergefn));
			UPair U;

			libmaus2::bambam::BamAccessor::unique_ptr_type tptr(new libmaus2::bambam::BamAccessor(bamfn,BI.first));
			libmaus2::bambam::BamAccessor & BD = *tptr;
			libmaus2::dazzler::align::Overlap OVLgen;

			uint64_t nextid = BI.first;

			while ( UPP->getNext(U) )
			{
				while ( nextid < U.key )
				{
					libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pgen(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(genfn,nextid,nextid+1));
					uint64_t miscnt = 0;
					while ( Pgen->getNextOverlap(OVLgen) )
					{
						int64_t const aspan = OVLgen.path.aepos - OVLgen.path.abpos;
						if ( aspan >= minlen )
							++miscnt;
					}
					printResultLine(nextid, miscnt, 0 /* gotcnt */, 0 /* extracnt */, header, BD[nextid], 0 /* mintrue */, 0 /* avgtrue */, 0 /* maxtrue */, 0, std::string("complete miss"));

					// nextid is missing
					++nextid;
				}

				assert ( nextid == U.key );
				++nextid;
			}
			while ( nextid < BI.second )
			{
				libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pgen(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(genfn,nextid,nextid+1));
				uint64_t miscnt = 0;
				while ( Pgen->getNextOverlap(OVLgen) )
				{
					int64_t const aspan = OVLgen.path.aepos - OVLgen.path.abpos;
					if ( aspan >= minlen )
						++miscnt;
				}
				printResultLine(nextid, miscnt, 0 /* gotcnt */, 0 /* extracnt */, header, BD[nextid], 0 /* mintrue */, 0 /* avgtrue */, 0 /* maxtrue */, 0, "complete miss");

				// nextid is missing

				++nextid;
			}
		}

		libmaus2::aio::SerialisedPeeker<UPair>::unique_ptr_type UPP(new libmaus2::aio::SerialisedPeeker<UPair>(tpmergefn));
		libmaus2::parallel::PosixSpinLock UPPlock;
		UPair U;
		uint64_t memfree = memfull;
		libmaus2::parallel::LockedBool running(false);
		libmaus2::parallel::TerminatableSynchronousQueue<UPair> Q;


		while ( UPP->peekNext(U) )
		{
			bool const ok = U.value <= memfull;
			if ( !ok )
			{
				std::cerr << "[E] total memory " << memfull << " too small for " << U.value << std::endl;
				assert ( ok );
			}

			if ( U.value <= memfree )
			{
				UPP->getNext(U);
				memfree -= U.value;
				Q.enque(U);
				running.set(true);
			}
			else
			{
				break;
			}
		}

		libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type Amark;
		libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type Amis;
		libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type Akeep;

		if ( flagmark )
		{
			libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type tptr(new libmaus2::dazzler::align::AlignmentWriterArray(tmpfilebase + "_marktmp", numthreads, tspace));
			Amark = UNIQUE_PTR_MOVE(tptr);
		}
		if ( flagmis )
		{
			libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type tptr(new libmaus2::dazzler::align::AlignmentWriterArray(tmpfilebase + "_mistmp", numthreads, tspace));
			Amis = UNIQUE_PTR_MOVE(tptr);
		}
		if ( flagkeep )
		{
			libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type tptr(new libmaus2::dazzler::align::AlignmentWriterArray(tmpfilebase + "_keeptmp", numthreads, tspace));
			Akeep = UNIQUE_PTR_MOVE(tptr);
		}

		libmaus2::autoarray::AutoArray < libmaus2::bambam::BamAccessor::unique_ptr_type > ABA(numthreads);
		for ( uint64_t i = 0; i < ABA.size(); ++i )
		{
			libmaus2::bambam::BamAccessor::unique_ptr_type tptr(new libmaus2::bambam::BamAccessor(bamfn,minaread));
			ABA[i] = UNIQUE_PTR_MOVE(tptr);
		}

		uint64_t const minsigtp = (minsiglen + tspace - 1) / tspace;

		#pragma omp parallel num_threads(numthreads)
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			libmaus2::dazzler::align::Overlap OVLgen;
			libmaus2::dazzler::align::Overlap OVLin;
			libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
			libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPVQ;
			libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TraceBlock> TBV;
			libmaus2::bambam::BamAccessor & BD = *(ABA[tid]);

			typedef libmaus2::util::FiniteSizeHeap<
				std::pair<uint64_t,uint64_t>,
				std::greater< std::pair<uint64_t,uint64_t> > > fsh_type;
			typedef fsh_type::shared_ptr_type fsh_ptr_type;

			libmaus2::autoarray::AutoArray < fsh_ptr_type > Aheap;
			uint64_t Aheapo = 0;

			while ( running.get() )
			{
				UPair U;

				try
				{
					U = Q.deque();

					#if 0
					{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << U.key << "," << U.value << std::endl;
					}
					#endif

					libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(infn,U.key,U.key+1));
					libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pgen(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(genfn,U.key,U.key+1));
					uint64_t gentraceid = 0;
					uint64_t qtraceid = 0;
					std::vector<uint64_t> qtrue;

					for ( uint64_t i = 0; i < Aheapo; ++i )
						Aheap[i]->clear();

					int64_t const aread = U.key;
					std::vector<bool> Bgen;

					uint64_t extracnt = 0;

					while (
						Plas->peekNextOverlap(OVLin)
						||
						Pgen->peekNextOverlap(OVLgen)
					)
					{
						if (
							(!(Pgen->peekNextOverlap(OVLgen)))
							||
							(
								Plas->peekNextOverlap(OVLin)
								&&
								OVLin.bread < OVLgen.bread
							)
						)
						{
							// extra
							bool const ok = Plas->getNextOverlap(OVLin);
							assert ( OVLin.aread == aread );
							assert ( ok );

							if ( verbose )
							{
								libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
								std::cerr << "[V] extra " << aread << " -> " << OVLin.bread << " span " << (OVLin.path.aepos - OVLin.path.abpos) << std::endl;
							}

							uint64_t const tbvo = OVLin.getTraceBlocks(tspace, TBV, 0, true /* only full blocks */);

							uint64_t const id = qtraceid++;

							if ( tbvo >= minsigtp )
							{
								uint64_t errsum = 0;

								for ( uint64_t i = 0; i < minsigtp-1; ++i )
									errsum += TBV[i].err;

								for ( uint64_t i = 0; i < (tbvo - minsigtp + 1); ++i )
								{
									errsum += TBV[i + minsigtp - 1].err;

									assert ( TBV[i].A.first % tspace == 0 );
									uint64_t const p = TBV[i].A.first / tspace;

									while ( !(p < Aheapo) )
									{
										fsh_ptr_type tptr(new fsh_type(seqdepth));
										Aheap.push(Aheapo,tptr);
									}
									assert ( p < Aheapo );

									if ( Aheap[p]->full() )
									{
										assert ( ! Aheap[p]->empty() );
										if ( errsum < Aheap[p]->top().first )
										{
											Aheap[p]->pop();
											Aheap[p]->push(std::pair<uint64_t,uint64_t>(errsum,id));
										}
									}
									else
									{
										Aheap[p]->push(std::pair<uint64_t,uint64_t>(errsum,id));
									}

									#if 0
									uint64_t errcheck = 0;
									for ( uint64_t j = 0; j < minsigtp; ++j )
										errcheck += TBV[i+j].err;
									assert ( errcheck == errsum );
									#endif

									errsum -= TBV[i].err;
								}
							}

							extracnt++;

							// put in mark file
							if ( Amark )
								(*Amark)[tid].put(OVLin);
						}
						else if (
							(!(Plas->peekNextOverlap(OVLin)))
							||
							(
								Pgen->peekNextOverlap(OVLgen)
								&&
								OVLgen.bread < OVLin.bread
							)
						)
						{
							// missing
							bool const ok = Pgen->getNextOverlap(OVLgen);
							assert ( OVLgen.aread == aread );
							assert ( ok );
							gentraceid += 1;

							int64_t const aspan = OVLgen.path.aepos - OVLgen.path.abpos;

							if ( aspan >= minlen )
							{
								Bgen.push_back(false);

								if ( verbose )
								{
									libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
									std::cerr << "[V] missing " << aread << " -> " << OVLgen.bread << " span " << aspan << " no match" << std::endl;
								}
							}
							else
							{
								Bgen.push_back(true);
							}

						}
						else
						{
							assert ( OVLin.aread == aread );
							assert ( OVLgen.aread == aread );

							assert ( Plas->peekNextOverlap(OVLin) && Pgen->peekNextOverlap(OVLgen) );

							int64_t const bread = OVLin.bread;

							uint64_t o = 0;
							for ( ; Pgen->peekNextOverlap(OVLgen) && OVLgen.bread == bread; ++gentraceid )
							{
								o = OVLgen.getTracePoints(tspace, gentraceid /* trace id */, TPV, o);
								bool const ok = Pgen->getNextOverlap(OVLgen);
								assert ( ok );

								int64_t const aspan = OVLgen.path.aepos - OVLgen.path.abpos;

								if ( aspan >= minlen )
									Bgen.push_back(false);
								else
									Bgen.push_back(true);
							}

							std::sort(TPV.begin(),TPV.begin()+o);

							#if 0
							for ( uint64_t i = 0; i < o; ++i )
								std::cerr << "TPV[]=" << TPV[i].apos << "," << TPV[i].bpos << std::endl;
							#endif

							//std::cerr << "o=" << o << std::endl;
							uint64_t qfound = 0;


							for ( ; Plas->peekNextOverlap(OVLin) && OVLin.bread == bread; ++qtraceid )
							{
								bool const ok = Plas->getNextOverlap(OVLin);
								assert ( ok );

								uint64_t const oo = OVLin.getTracePoints(tspace, qtraceid /* trace id */, TPVQ, 0);
								uint64_t const tbvo = OVLin.getTraceBlocks(tspace, TBV, 0, true /* only full blocks */);

								if ( tbvo >= minsigtp )
								{
									uint64_t errsum = 0;

									for ( uint64_t i = 0; i < minsigtp-1; ++i )
										errsum += TBV[i].err;

									for ( uint64_t i = 0; i < (tbvo - minsigtp + 1); ++i )
									{
										errsum += TBV[i + minsigtp - 1].err;

										assert ( TBV[i].A.first % tspace == 0 );
										uint64_t const p = TBV[i].A.first / tspace;

										while ( !(p < Aheapo) )
										{
											fsh_ptr_type tptr(new fsh_type(seqdepth));
											Aheap.push(Aheapo,tptr);
										}
										assert ( p < Aheapo );

										if ( Aheap[p]->full() )
										{
											assert ( ! Aheap[p]->empty() );
											if ( errsum < Aheap[p]->top().first )
											{
												Aheap[p]->pop();
												Aheap[p]->push(std::pair<uint64_t,uint64_t>(errsum,qtraceid));
											}
										}
										else
										{
											Aheap[p]->push(std::pair<uint64_t,uint64_t>(errsum,qtraceid));
										}

										#if 0
										uint64_t errcheck = 0;
										for ( uint64_t j = 0; j < minsigtp; ++j )
											errcheck += TBV[i+j].err;
										assert ( errcheck == errsum );
										#endif

										errsum -= TBV[i].err;
									}
								}

								#if 0
								std::cerr << "oo=" << oo << std::endl;

								for ( uint64_t i = 0; i < oo; ++i )
									std::cerr << "TPVQ[]=" << TPVQ[i].apos << "," << TPVQ[i].bpos << std::endl;
								#endif

								struct TracePointPosComparator
								{
									bool operator()(libmaus2::dazzler::align::TracePoint const & A, libmaus2::dazzler::align::TracePoint const & B) const
									{
										if ( A.apos != B.apos )
											return A.apos < B.apos;
										else
											return A.bpos < B.bpos;
									}
								};

								bool found = false;
								uint64_t genid = 0;
								for ( uint64_t i = 0; !found && i < oo; ++i )
								{
									typedef libmaus2::dazzler::align::TracePoint tp_type;
									std::pair<tp_type const *, tp_type const *> E = ::std::equal_range(TPV.begin(),TPV.begin()+o,TPVQ[i],TracePointPosComparator());

									bool const lfound = E.second != E.first;

									if ( lfound )
									{
										found = true;
										genid = E.first->id;
										Bgen[genid] = true;
									}
								}
								if ( found )
								{
									if ( Akeep )
										(*Akeep)[tid].put(OVLin);

									OVLin.flags |= libmaus2::dazzler::align::Overlap::getTrueFlag();

									if ( Amark )
										(*Amark)[tid].put(OVLin);

									qfound++;
									qtrue.push_back(qtraceid);
								}
								else
								{
									if ( Amark )
										(*Amark)[tid].put(OVLin);

									++extracnt;

									#if 0
									if ( verbose )
									{
										libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
										std::cerr << "ref " << aread << " " << OVLin.bread << " traceblocks\n";
										for ( uint64_t i = 0; i < o; ++i )
											std::cerr << "(" << TPV[i].apos << "," << TPV[i].bpos << ")";
										std::cerr << std::endl;
										#if 0
										std::cerr << "query traceblocks\n";
										for ( uint64_t i = 0; i < oo; ++i )
											std::cerr << "(" << TPVQ[i].apos << "," << TPVQ[i].bpos << ")";
										std::cerr << std::endl;
										#endif
									}
									#endif
									// not found
								}

								// ZZZ fill
							}

							#if 0
							{
								libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
								std::cerr << OVLgen.aread << " " << OVLgen.bread << " qfound=" << qfound << " gen " << gentraceid << std::endl;
							}
							#endif
						}
					}

					double maxtrue = 0;
					double mintrue = std::numeric_limits<double>::max();
					double avgtrue = 0.0;
					double avgcnt = 0;

					for ( uint64_t i = 0; i < Aheapo; ++i )
					{
						assert ( Aheap[i] );

						if ( ! Aheap[i]->empty() )
						{
							uint64_t truecnt = 0;
							uint64_t falsecnt = 0;

							while ( ! Aheap[i]->empty() )
							{
								std::pair<uint64_t,uint64_t> const P = Aheap[i]->pop();

								#if 0
								std::cerr << "Aheap[" << i << "] " << P.first << "," << P.second << std::endl;
								#endif

								std::vector<uint64_t>::const_iterator it =
									::std::lower_bound(qtrue.begin(),qtrue.end(),P.second);
								if ( it != qtrue.end() && *it == P.second )
								{
									truecnt++;
								}
								else
								{
									falsecnt++;
								}
							}

							double const truerate = static_cast<double>(truecnt) / (truecnt+falsecnt);

							maxtrue = std::max(maxtrue,truerate);
							mintrue = std::min(mintrue,truerate);
							avgtrue += truerate;
							avgcnt += 1;

							// std::cerr << "aread=" << aread << " window=" << i << " truecnt=" << truecnt << " falsecnt=" << falsecnt << std::endl;
						}
					}

					avgtrue = avgcnt ? (avgtrue/avgcnt) : avgtrue;

					#if 0
					if ( avgcnt )
					{
						std::cerr << "aread=" << aread << " mintrue=" << mintrue << " maxtrue=" << maxtrue << " avgtrue=" << avgtrue/avgcnt << std::endl;
					}
					#endif


					{
						libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Pgen(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(genfn,U.key,U.key+1));

						uint64_t miscnt = 0;
						uint64_t gotcnt = 0;
						for ( uint64_t i = 0; i < Bgen.size(); ++i )
						{
							bool const ok = Pgen->getNextOverlap(OVLgen);
							assert ( ok );

							if ( ! Bgen[i] )
							{
								miscnt++;
								if ( Amis )
									(*Amis) [ tid ] . put(OVLgen);
								// std::cerr << OVLgen.aread << " " << OVLgen.bread << " missing " << i << std::endl;

								if ( verbose )
								{
									libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
									std::cerr << "[V] missing " << aread << " -> " << OVLgen.bread << " span " << (OVLgen.path.aepos - OVLgen.path.abpos) << " wrong match " << OVLgen << std::endl;

									uint64_t const o = OVLgen.getTracePoints(tspace, gentraceid /* trace id */, TPV, 0);
									std::cerr << "ref " << aread << " " << OVLgen.bread << " traceblocks\n";
									for ( uint64_t i = 0; i < o; ++i )
										std::cerr << "(" << TPV[i].apos << "," << TPV[i].bpos << ")";
									std::cerr << std::endl;

								}
							}
							else
							{
								gotcnt++;
							}
						}
						printResultLine(aread, miscnt, gotcnt, extracnt, header, BD[aread], mintrue, avgtrue, maxtrue, avgcnt,"regular");
					}

					#if 0
					for ( uint64_t i = 0; i < qtrue.size(); ++i )
					{
						std::cerr << "qtrue[" << i << "]=" << qtrue[i] << std::endl;
					}
					#endif


					{
						libmaus2::parallel::ScopePosixSpinLock slock(UPPlock);
						memfree += U.value;

						while ( UPP->peekNext(U) )
						{
							bool const ok = U.value <= memfull;
							if ( !ok )
							{
								libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
								std::cerr << "[E] total memory " << memfull << " too small for " << U.value << std::endl;
								assert ( ok );
							}


							if ( U.value <= memfree )
							{
								UPP->getNext(U);
								memfree -= U.value;
								Q.enque(U);
							}
							else
							{
								break;
							}
						}
						if ( ! UPP->peekNext(U) )
						{
							Q.terminate();
						}
					}
				}
				catch(std::exception const & ex)
				{
					running.set(false);
				}

			}
		}

		UPP.reset();
		libmaus2::aio::FileRemoval::removeFile(tpmergefn);

		if ( Amis )
			Amis->merge(infn + ".mis.las", tmpfilebase + ".mis.las.tmp");
		if ( Akeep )
			Akeep->merge(infn + ".keep.las", tmpfilebase + ".keep.las.tmp");
		if ( Amark )
			Amark->merge(infn + ".mark.las", tmpfilebase + ".mark.las.tmp");

		#if 0
		while ( UPP.getNext(U) )
		{
			std::cerr << U.key << "\t" << U.value << std::endl;
		}
		#endif
	}

	return failed ? EXIT_FAILURE : EXIT_SUCCESS;
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
	optionMap . push_back ( std::pair < std::string, std::string >("minsiglen", formatRHS("minimum length of ground truth alignments considered",getDefaultMinSigLen())));
	optionMap . push_back ( std::pair < std::string, std::string >("seqdepth", formatRHS("sequencing depth",getDefaultSeqDepth())));
	optionMap . push_back ( std::pair < std::string, std::string >("verbose", formatRHS("verbosity",false)));
	optionMap . push_back ( std::pair < std::string, std::string >("tlow", formatRHS("minimum read id considered",std::numeric_limits<int64_t>::min())));
	optionMap . push_back ( std::pair < std::string, std::string >("thigh", formatRHS("minimum read id considered",std::numeric_limits<int64_t>::max())));
	optionMap . push_back ( std::pair < std::string, std::string >("mark", formatRHS("write marked true alignment file",false)));
	optionMap . push_back ( std::pair < std::string, std::string >("mis", formatRHS("write missed true alignment file",false)));
	optionMap . push_back ( std::pair < std::string, std::string >("keep", formatRHS("write kept true alignment file",false)));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("minlen", formatRHS("minimum length for alignments",getDefaultMinLen())));

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
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] in.db gen.las in.bam in.las ...\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return checklas(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
