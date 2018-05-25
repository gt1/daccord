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

#include <config.h>

#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/OverlapProperCheck.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <RepeatIdComparator.hpp>
#include <libmaus2/util/IntervalTree.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static double getDefaultTermVal()
{
	return 0.35;
}

static uint64_t getDefaultDepthThres()
{
	return 10;
}

struct RepeatPosComparator
{
        bool operator()(Repeat const & A, Repeat const & B) const
        {
        	if ( A.id != B.id )
        		return A.id < B.id;
		else
	                return A.abpos < B.abpos;
        }
};

struct RepeatIndexer : public libmaus2::sorting::SerialisingSortingBufferedOutputFile<Repeat,RepeatPosComparator>::IndexCallback
{
	std::vector< std::pair<uint64_t,uint64_t> > V;
	uint64_t cnt;
	int64_t previd;

	RepeatIndexer() : cnt(0), previd(-1) {}
	RepeatIndexer(std::istream & in) : previd(-1)
	{
		deserialise(in);
	}
	RepeatIndexer(std::string const & fn) : previd(-1)
	{
		deserialise(fn);
	}

	void serialise(std::string const & fn) const
	{
		libmaus2::aio::OutputStreamInstance OSI(fn);
		serialise(OSI);
	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,cnt);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,V.size());
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			libmaus2::util::NumberSerialisation::serialiseNumber(out,V[i].first);
			libmaus2::util::NumberSerialisation::serialiseNumber(out,V[i].second);
		}

		return out;
	}

	void deserialise(std::string const & fn)
	{
		libmaus2::aio::InputStreamInstance ISI(fn);
		deserialise(ISI);
	}

	std::istream & deserialise(std::istream & in)
	{
		cnt = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		V.resize(libmaus2::util::NumberSerialisation::deserialiseNumber(in));
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			V[i].first = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
			V[i].second = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		}

		return in;
	}

	void operator()(Repeat const & A, uint64_t const p)
	{
		++cnt;

		if ( static_cast<int64_t>(A.id) != previd )
		{
			V.push_back(std::pair<uint64_t,uint64_t>(A.id,p));
			// std::cerr << A.id << "," << A.abpos << "," << A.aepos << "\t" << p << std::endl;

			previd = A.id;
		}
	}

	std::vector < std::pair<uint64_t,uint64_t> > split(uint64_t const tparts)
	{
		uint64_t const partsize = std::max ( static_cast<uint64_t>(1), (V.size() + tparts - 1)/tparts );
		uint64_t const parts = (V.size() + partsize - 1)/partsize;
		std::vector < std::pair<uint64_t,uint64_t> > P;
		std::vector < std::pair<uint64_t,uint64_t> > Q;

		for ( uint64_t i = 0; i < parts; ++i )
		{
			uint64_t const low = i*partsize;
			uint64_t const high = std::min(low+partsize,static_cast<uint64_t>(V.size()));
			P.push_back(std::pair<uint64_t,uint64_t>(low,high));
		}

		for ( uint64_t i = 0; i < P.size(); ++i )
		{
			uint64_t const id = V[P[i].first].first;
			uint64_t const p = V[P[i].first].second;
			Q.push_back(std::pair<uint64_t,uint64_t>(id,p));
		}
		Q.push_back(std::pair<uint64_t,uint64_t>(std::numeric_limits<uint64_t>::max(),0));

		return Q;
	}
};


struct Handler
{
	typedef Handler this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	std::set<uint64_t> active;
	libmaus2::util::FiniteSizeHeap< std::pair<int64_t,uint64_t> > H;

	Handler() : H(1024) {}

	void handleVector(
		std::ostream & out,
		int64_t const prevread,
		std::vector< std::pair<int64_t,int64_t> > & V,
		uint64_t thres
	)
	{
		if ( V.size() )
		{
			// sort repeat evidence
			std::sort(V.begin(),V.end());

			std::vector < std::pair<int64_t,int64_t> > IV;

			int64_t istart = 0;
			for ( uint64_t i = 0; i < V.size(); ++i )
			{
				// remove inactive repeat regions before start of V[i]
				while ( !H.empty() && H.top().first <= V[i].first )
				{
					std::pair<int64_t,uint64_t> P = H.pop();

					if ( P.first > istart )
					{
						if ( active.size() >= thres )
						{
							// std::cerr << prevread << " " << istart << "," << P.first << "," << active.size() << std::endl;

							if ( IV.size() && IV.back().second == istart )
								IV.back().second = P.first;
							else
								IV.push_back(std::pair<int64_t,int64_t>(istart,P.first));
						}
					}
					istart = P.first;

					active.erase(P.second);
				}

				if ( active.size() && V[i].first > istart )
				{
					if ( active.size() >= thres )
					{
						// std::cerr << prevread << " " << istart << "," << V[i].first << "," << active.size() << std::endl;

						if ( IV.size() && IV.back().second == istart )
							IV.back().second = V[i].first;
						else
							IV.push_back(std::pair<int64_t,int64_t>(istart,V[i].first));
					}
				}

				H.pushBump(std::pair<int64_t,uint64_t>(V[i].second,i));
				active.insert(i);

				istart = V[i].first;
			}

			while ( !H.empty() )
			{
				std::pair<int64_t,uint64_t> P = H.pop();

				if ( P.first > istart )
				{
					if ( active.size() >= thres )
					{
						// std::cerr << prevread << " " << istart << "," << P.first << "," << active.size() << std::endl;

						if ( IV.size() && IV.back().second == istart )
							IV.back().second = P.first;
						else
							IV.push_back(std::pair<int64_t,int64_t>(istart,P.first));
					}
				}
				istart = P.first;

				active.erase(P.second);
			}

			// write repeat regions as triples (readid,from,to)
			for ( uint64_t i = 0; i < IV.size(); ++i )
			{
				libmaus2::util::NumberSerialisation::serialiseNumber(out,prevread);
				libmaus2::util::NumberSerialisation::serialiseNumber(out,IV[i].first);
				libmaus2::util::NumberSerialisation::serialiseNumber(out,IV[i].second);
				// std::cerr << prevread << " " << IV[i].first << "," << IV[i].second << std::endl;
			}

			V.resize(0);
		}
	}
};


static void filterRepeats(std::string const mergetmp)
{
	std::cerr << "[V] filtering" << std::endl;
	{
		libmaus2::aio::SerialisedPeeker<Repeat> SP(mergetmp);
		Repeat R;

		std::string const tmpfn = mergetmp + ".filtertmp";
		libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(tmpfn));

		while ( SP.peekNext(R) )
		{
			uint64_t const id = R.id;

			std::vector < libmaus2::math::IntegerInterval<int64_t> > VI;

			uint64_t incnt = 0;
			while ( SP.peekNext(R) && R.id == id )
			{
				SP.getNext(R);
				VI.push_back(libmaus2::math::IntegerInterval<int64_t>(R.abpos,R.aepos-1));
				VI = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VI);
				++incnt;
			}

			for ( uint64_t i = 0; i < VI.size(); ++i )
				Repeat(id,VI[i].from,VI[i].to+1).serialise(*OSI);
		}

		OSI->flush();
		OSI.reset();

		libmaus2::aio::OutputStreamFactoryContainer::rename(tmpfn,mergetmp);
	}
	std::cerr << "[V] filtered" << std::endl;
}

struct RepeatTranslation
{
	static std::string computeTranslation(
		std::string const tmpfilebase,
		uint64_t const numthreads,
		std::vector<uint64_t> const & Vchange,
		std::vector<Repeat> const & VR,
		libmaus2::dazzler::db::DatabaseFile const & DB,
		libmaus2::util::IntervalTree const & IT,
		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> const & Adalindex,
		std::vector<std::string> const & Vinfn,
		bool const precise,
		std::vector<uint64_t> const & RL,
		std::string const & mergetmp
	)
	{
		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vinfn);

		// temp output files
		std::vector<std::string> Vtmp(numthreads);
		libmaus2::autoarray::AutoArray <
			libmaus2::aio::OutputStreamInstance::unique_ptr_type
		> Aout(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			std::ostringstream ostr;
			ostr << tmpfilebase << "_thread_" << i << ".tmp";
			std::string const fn = ostr.str();
			Vtmp[i] = fn;
			libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
			libmaus2::aio::OutputStreamInstance::unique_ptr_type tptr(new libmaus2::aio::OutputStreamInstance(fn));
			Aout[i] = UNIQUE_PTR_MOVE(tptr);
		}

		uint64_t volatile ufinished = 0;
		libmaus2::parallel::PosixSpinLock ufinishedlock;

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t q = 0; q < Vchange.size(); ++q )
		{
			uint64_t const z = Vchange[q];

			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			std::ostream & out = *(Aout[tid]);

			typedef std::vector<Repeat>::const_iterator it;
			std::pair<it,it> IP = std::equal_range(VR.begin(),VR.end(),Repeat(z),RepeatIdComparator());
			uint64_t const low  = IP.first - VR.begin();
			uint64_t const high = IP.second - VR.begin();

			std::basic_string<uint8_t> const ua = DB.getu(z);

			if ( high-low )
			{
				struct IRange
				{
					uint64_t id;
					uint64_t from;
					uint64_t to;

					IRange() {}
					IRange(
						uint64_t rid,
						uint64_t rfrom,
						uint64_t rto
					) : id(rid), from(rfrom), to(rto) {}

					uint64_t getFrom() const
					{
						return from;
					}
					uint64_t getTo() const
					{
						return to;
					}
				};

				libmaus2::geometry::RangeSet<IRange> RS(VR.back().aepos);
				for ( uint64_t i = low; i < high; ++i )
					RS.insert(IRange(i,VR[i].abpos,VR[i].aepos));

				uint64_t const lasid = IT.find(z);
				// libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(Vinfn[lasid],z,z+1));

				uint64_t const zstart = (*Adalindex[lasid])[z];
				uint64_t const zend = (*Adalindex[lasid])[z+1];

				libmaus2::aio::InputStreamInstance ISI(Vinfn[lasid]);
				ISI.clear();
				ISI.seekg(zstart);

				libmaus2::dazzler::align::SimpleOverlapParser SOP(ISI,tspace,16*1024*1024,libmaus2::dazzler::align::OverlapParser::overlapparser_do_split,zend-zstart);

				bool const small = libmaus2::dazzler::align::AlignmentFile::tspaceToSmall(tspace);
				libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > Apath;

				while ( SOP.parseNextBlock() )
				{
					libmaus2::dazzler::align::OverlapData & data = SOP.getData();

					for ( uint64_t i = 0; i < data.size(); ++i )
					{
						std::pair<uint8_t const *, uint8_t const *> const P = data.getData(i);
						assert ( libmaus2::dazzler::align::OverlapData::getARead(P.first) == static_cast<int64_t>(z) );
						int32_t const abpos = libmaus2::dazzler::align::OverlapData::getABPos(P.first);
						int32_t const aepos = libmaus2::dazzler::align::OverlapData::getAEPos(P.first);

						std::vector<IRange const *> IRV;
						RS.search(IRange(0,abpos,aepos),IRV);

						for ( uint64_t i = 0; i < IRV.size(); ++i )
						{
							Repeat const & R = VR[IRV[i]->id];

							libmaus2::math::IntegerInterval<int64_t> const IA(R.abpos,R.aepos-1);
							libmaus2::math::IntegerInterval<int64_t> const IB(abpos,aepos-1);
							libmaus2::math::IntegerInterval<int64_t> const IC = IA.intersection(IB);

							if ( precise )
							{
								int32_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(P.first);
								bool const inv = libmaus2::dazzler::align::OverlapData::getInverseFlag(P.first);

								libmaus2::lcs::NP np;
								libmaus2::lcs::AlignmentTraceContainer ATC;
								uint64_t const plen = libmaus2::dazzler::align::OverlapData::decodeTraceVector(P.first,Apath,small);
								int32_t const bbpos = libmaus2::dazzler::align::OverlapData::getBBPos(P.first);
								int32_t const bepos = libmaus2::dazzler::align::OverlapData::getBEPos(P.first);

								libmaus2::dazzler::align::Overlap::computeTrace(Apath.begin(),plen,
									abpos,aepos,
									bbpos,bepos,
									ua.c_str(),
									DB.getu(bread,inv).c_str(),
									tspace,
									ATC,
									np);

								// offset of intersection on trace
								uint64_t const aoff = IC.from - IB.from;

								// advance on trace
								std::pair<uint64_t,uint64_t> const advB = libmaus2::lcs::AlignmentTraceContainer::advanceA(ATC.ta,ATC.te,aoff);
								// string length used
								std::pair<uint64_t,uint64_t> const slB  = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ATC.ta,ATC.ta + advB.second);
								// advance by length of IC on trace
								std::pair<uint64_t,uint64_t> const advE = libmaus2::lcs::AlignmentTraceContainer::advanceA(ATC.ta+advB.second,ATC.te,IC.to-IC.from+1);
								// string length used
								std::pair<uint64_t,uint64_t> const slE  = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ATC.ta + advB.second,ATC.ta + advB.second + advE.second);

								// start and end of b
								uint64_t bstart = bbpos + slB.second;
								uint64_t bend   = bstart + slE.second;

								// apply RC if necessary
								if ( inv )
								{
									std::swap(bstart,bend);

									bstart = RL[bread] - bstart;
									bend   = RL[bread] - bend;
								}

								// check if repeat is already present
								std::pair<it,it> const BIP = std::equal_range(VR.begin(),VR.end(),Repeat(bread),RepeatIdComparator());
								bool drop = false;
								for ( it c = BIP.first; c != BIP.second; ++c )
									if ( c->abpos <= bstart && c->aepos >= bend )
										drop = true;

								// write if repeat is not already present
								if ( ! drop )
								{
									// construct and serialise repeat object
									Repeat RN(bread,bstart,bend);
									RN.serialise(out);
								}
								#if 0
								else
								{
									Repeat RN(bread,bstart,bend);
									std::cerr << "dropping Repeat=" << RN.id << "," << RN.abpos << "," << RN.aepos << std::endl;
								}
								#endif

								#if 0
								std::pair<it,it> IP = std::equal_range(VR.begin(),VR.end(),Repeat(OVL.bread),RepeatIdComparator());
								for ( it i = IP.first; i != IP.second; ++i )
									std::cerr << "\t" << i->abpos << "," << i->aepos << std::endl;
								#endif
							}
							else // imprecise
							{
								// align to trace points
								uint64_t const uleft = std::min(static_cast<uint64_t>(((IC.from + tspace - 1)/tspace)*tspace),static_cast<uint64_t>(RL[z]));
								uint64_t const uright = std::min(static_cast<uint64_t>( ((IC.to+1)/tspace)*tspace ),static_cast<uint64_t>(RL[z]));

								if ( uright > uleft )
								{
									assert ( static_cast<int64_t>(uleft) >= abpos );
									assert ( static_cast<int64_t>(uright) <= aepos );
									assert ( uleft % tspace == 0 );
									assert ( (uright % tspace == 0) || (uright == RL[z]) );

									uint64_t const plen = libmaus2::dazzler::align::OverlapData::decodeTraceVector(P.first,Apath,small);

									int32_t const bbpos = libmaus2::dazzler::align::OverlapData::getBBPos(P.first);
									int32_t const bepos = libmaus2::dazzler::align::OverlapData::getBEPos(P.first);
									int32_t const bread = libmaus2::dazzler::align::OverlapData::getBRead(P.first);

									libmaus2::dazzler::align::Overlap::OffsetInfo offleft = libmaus2::dazzler::align::Overlap::getBforAOffset(
										tspace,
										abpos,
										aepos,
										bbpos,
										bepos,
										uleft,
										Apath.begin(),
										plen
									);
									libmaus2::dazzler::align::Overlap::OffsetInfo offright = libmaus2::dazzler::align::Overlap::getBforAOffset(
										tspace,
										abpos,
										aepos,
										bbpos,
										bepos,
										uright,
										Apath.begin(),
										plen
									);

									assert ( offleft.apos == static_cast<int64_t>(uleft) );
									assert ( offright.apos == static_cast<int64_t>(uright) );

									// start and end of b
									uint64_t bstart = offleft.bpos;
									uint64_t bend   = offright.bpos;

									// apply RC if necessary
									if ( libmaus2::dazzler::align::OverlapData::getInverseFlag(P.first) )
									{
										std::swap(bstart,bend);

										bstart = RL[bread] - bstart;
										bend   = RL[bread] - bend;
									}

									// check if repeat is already present
									std::pair<it,it> const BIP = std::equal_range(VR.begin(),VR.end(),Repeat(bread),RepeatIdComparator());
									bool drop = false;
									for ( it c = BIP.first; c != BIP.second; ++c )
										if ( c->abpos <= bstart && c->aepos >= bend )
											drop = true;

									// write if repeat is not already present
									if ( ! drop )
									{
										// construct and serialise repeat object
										Repeat RN(bread,bstart,bend);
										RN.serialise(out);
									}

									// std::cerr << "read " << OVL.bread << " overlaps read " << z << " in " << R.abpos << "," << R.aepos << " with " << OVL.path.abpos << "," << OVL.path.aepos << std::endl;
								}
							}
						}

					}
				}
			}

			{
				bool print = false;
				uint64_t vprint = 0;

				ufinishedlock.lock();

				if ( ++ufinished % 100 == 0 )
				{
					vprint = ufinished;
					print = true;
				}

				ufinishedlock.unlock();

				if ( print )
				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << "[V] p" << vprint << "/" << Vchange.size() << std::endl;
				}
			}
		}

		// flush output files
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			Aout[i]->flush();
			Aout[i].reset();
		}

		// merge to a single file
		std::cerr << "[V] merging files" << std::endl;
		libmaus2::util::TempFileRemovalContainer::addTempFile(mergetmp);
		RepeatIndexer indexer;
		libmaus2::sorting::SerialisingSortingBufferedOutputFile<Repeat,RepeatPosComparator>::reduce(Vtmp,indexer,mergetmp,512*1024*1024 /* blocksize */);
		for ( uint64_t i = 0; i < Vtmp.size(); ++i )
			libmaus2::aio::FileRemoval::removeFile(Vtmp[i]);
		libmaus2::util::TempFileRemovalContainer::addTempFile(mergetmp + ".index");
		indexer.serialise(mergetmp + ".index");
		std::cerr << "[V] merged files, cnt=" << indexer.cnt << std::endl;

		return mergetmp;
	}

	static std::string computeTranslation(libmaus2::util::ArgParser const & arg)
	{
		std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
		uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
		bool const precise = arg.uniqueArgPresent("p");

		std::string const outrep = arg[0];
		std::string const inrep = arg[1];
		std::string const changefn = inrep + ".change";
		std::string const dbfn = arg[2];

		std::vector<std::string> Vinfn;
		for ( uint64_t i = 3; i < arg.size(); ++i )
			Vinfn.push_back(arg[i]);

		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(Vinfn.size());

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < Vinfn.size(); ++i )
		{
			std::ostringstream memindexname;
			memindexname << "mem:index_" << i << ".las.idx";
			std::string const memindexfn = memindexname.str();
			libmaus2::util::GetFileSize::copy(libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(Vinfn[i]),memindexfn);

			libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
				new libmaus2::dazzler::align::DalignerIndexDecoder(Vinfn[i],memindexfn)
			);
			Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
		}

		uint64_t io = 0;
		::libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > IH;
		for ( uint64_t i = 0; i < Vinfn.size(); ++i )
		{
			if (
			! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]))
				||
				libmaus2::util::GetFileSize::isOlder(
					libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]),
					Vinfn[i]
				)
			)
			{
				libmaus2::dazzler::align::OverlapIndexer::constructIndex(Vinfn[i]);
				std::cerr << "[V] constructing index for " << Vinfn[i] << "\n";

			}

			int64_t const mini = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(Vinfn[i]);
			int64_t const maxi = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(Vinfn[i]);

			if ( mini >= 0 )
			{
				assert ( maxi >= mini );
				Vinfn[io] = Vinfn[i];
				IH.push(io,std::pair<uint64_t,uint64_t>(mini,maxi+1));
			}
		}
		Vinfn.resize(io);
		IH.resize(io);
		for ( uint64_t i = 1; i < io; ++i )
		{
			bool const ok = IH[i-1].second <= IH[i].first;
			if ( ! ok )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] LAS files or not given in increasing order of A read intervals" << std::endl;
				lme.finish();
				throw lme;
			}
		}
		for ( uint64_t i = 1; i < io; ++i )
			IH[i-1].second = IH[i].first;

		libmaus2::util::IntervalTree IT(IH,0,IH.size(),true);

		std::cerr << "[V] copying " << dbfn << " to memory...";
		std::vector<std::string> tracklist;
		tracklist.push_back("inqual");
		// libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbfn,"mem:db1prefix",&tracklist));
		libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
			libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbfn,&tracklist)
		);
		std::cerr << "done." << std::endl;
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
		PDB->computeTrimVector();

		libmaus2::dazzler::db::DatabaseFile & DB = *PDB;
		if ( DB.part != 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] partial databases are not supported" << std::endl;
			lme.finish();
			throw lme;
		}

		// get read length vector
		std::vector<uint64_t> RL;
		DB.getAllReadLengths(RL);

		// load change vector
		std::vector<uint64_t> Vchange;
		{
			libmaus2::aio::InputStreamInstance ISI(changefn);
			Vchange = libmaus2::util::NumberSerialisation::deserialiseNumberVector<uint64_t>(ISI);
		}

		uint64_t rlow = 0;
		uint64_t rhigh = Vchange.size();

		if ( arg.uniqueArgPresent("J") )
		{
			std::string const Js = arg["J"];
			std::istringstream istr(Js);
			int64_t Icnt;
			istr >> Icnt;

			if ( ! istr )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse " << Js << std::endl;
				lme.finish();
				throw lme;
			}

			int const c = istr.get();

			if ( ! istr || c == std::istream::traits_type::eof() || c != ',' )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse " << Js << std::endl;
				lme.finish();
				throw lme;
			}

			int64_t Idiv;
			istr >> Idiv;

			if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse " << Js << std::endl;
				lme.finish();
				throw lme;
			}

			if ( ! Idiv )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] divisor in J argument cannot be zero" << std::endl;
				lme.finish();
				throw lme;
			}

			uint64_t const span = rhigh-rlow;
			uint64_t const packsize = std::max ( (span + Idiv - 1)/Idiv, static_cast<uint64_t>(1) );

			rlow  = std::min(Icnt * packsize,rhigh);
			rhigh = std::min(rlow + packsize,rhigh);
		}

		Vchange = std::vector<uint64_t>(Vchange.begin()+rlow,Vchange.begin()+rhigh);

		libmaus2::aio::InputStreamInstance ISI(inrep);
		std::vector<Repeat> VR;
		while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
		{
			Repeat R(ISI);
			VR.push_back(R);
		}

		return computeTranslation(tmpfilebase,numthreads,Vchange,VR,DB,IT,Adalindex,Vinfn,precise,RL,outrep);
	}

	static std::string mergeTranslationParts(std::string const tmpfilebase, std::vector<std::string> const & transtmp)
	{
		std::cerr << "[V] concatenating translation parts" << std::endl;
		typedef libmaus2::aio::SerialisingSortingBufferedOutput<Repeat,RepeatPosComparator>::BlockDescriptor block_type;
		std::vector < block_type > transblocks;
		std::string const transcat = tmpfilebase + "_transcat";
		libmaus2::util::TempFileRemovalContainer::addTempFile(transcat);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type Ptranscat(new libmaus2::aio::OutputStreamInstance(transcat));
		uint64_t ptranscat = 0;
		for ( uint64_t zz = 0; zz < transtmp.size(); ++zz )
		{
			libmaus2::aio::InputStreamInstance::unique_ptr_type pISI(new libmaus2::aio::InputStreamInstance(transtmp[zz]));
			std::istream & ISI = *pISI;
			uint64_t const fs = libmaus2::util::GetFileSize::getFileSize(ISI);

			libmaus2::util::GetFileSize::copy(ISI,*Ptranscat,fs);

			RepeatIndexer indexer(transtmp[zz] + ".index");

			transblocks.push_back(block_type(indexer.cnt,ptranscat));

			ptranscat += fs;

			pISI.reset();
			libmaus2::aio::FileRemoval::removeFile(transtmp[zz]);
			libmaus2::aio::FileRemoval::removeFile(transtmp[zz] + ".index");
		}
		Ptranscat->flush();
		Ptranscat.reset();
		std::cerr << "[V] concatenated translation parts" << std::endl;

		std::cerr << "[V] merging translation parts" << std::endl;
		std::ostringstream mergetmpstr;
		mergetmpstr << tmpfilebase + "_thread_merge.tmp";
		std::string const mergetmp = mergetmpstr.str();

		{
			RepeatPosComparator RPC;
			typename libmaus2::sorting::SerialisingMergingReadBack<Repeat,RepeatPosComparator>::unique_ptr_type mergeptr(
				new libmaus2::sorting::SerialisingMergingReadBack<Repeat,RepeatPosComparator>(transcat,transblocks,RPC,1024)
			);
			Repeat R;

			libmaus2::aio::OutputStreamInstance::unique_ptr_type Pout(new libmaus2::aio::OutputStreamInstance(mergetmp));
			RepeatIndexer indexer;

			while ( mergeptr->getNext(R) )
			{
				uint64_t const p = Pout->tellp();
				indexer(R,p);
				R.serialise(*Pout);
			}

			Pout->flush();
			Pout.reset();

			indexer.serialise(mergetmp + ".index");
		}
		libmaus2::aio::FileRemoval::removeFile(transcat);
		std::cerr << "[V] merged translation parts" << std::endl;

		return mergetmp;
	}

	static void filterTranslationParts(std::string const & mergetmp, uint64_t const numthreads, std::string const & tmpfilebase)
	{
		RepeatIndexer indexer(mergetmp + ".index");
		std::cerr << "[V] loaded indexer, cnt " << indexer.cnt << std::endl;

		std::vector < std::pair<uint64_t,uint64_t> > const S = indexer.split(numthreads);

		uint64_t const slimit = S.size() ? (S.size()-1) : 0;
		std::vector < std::string > Vfiltertmpfn(slimit);

		std::cerr << "[V] filtering in parallel" << std::endl;
		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < slimit; ++i )
		{
			std::ostringstream fnostr;
			fnostr << tmpfilebase << "_filtertmp" << i;
			Vfiltertmpfn[i] = fnostr.str();
			libmaus2::util::TempFileRemovalContainer::addTempFile(Vfiltertmpfn[i]);
			libmaus2::aio::OutputStreamInstance OSI(Vfiltertmpfn[i]);
			libmaus2::aio::InputStreamInstance ISI(mergetmp);
			ISI.clear();
			ISI.seekg(S[i].second);
			libmaus2::aio::SerialisedPeeker<Repeat> SP(ISI);
			Repeat R;

			while ( SP.peekNext(R) && R.id < S[i+1].first )
			{
				uint64_t const id = R.id;

				std::vector < libmaus2::math::IntegerInterval<int64_t> > VI;
				while ( SP.peekNext(R) && R.id == id )
				{
					SP.getNext(R);
					VI.push_back(libmaus2::math::IntegerInterval<int64_t>(R.abpos,R.aepos-1));
					VI = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(VI);
				}

				for ( uint64_t i = 0; i < VI.size(); ++i )
					Repeat(id,VI[i].from,VI[i].to+1).serialise(OSI);
			}

			// std::cerr << "[" << S[i].first << "," << S[i+1].first << ")\t" << S[i].second << std::endl;
		}
		std::cerr << "[V] filtered in parallel" << std::endl;

		libmaus2::aio::FileRemoval::removeFile(mergetmp + ".index");

		std::cerr << "[V] concatenating" << std::endl;
		{
			libmaus2::aio::OutputStreamInstance OSI(mergetmp);

			for ( uint64_t i = 0; i < slimit; ++i )
			{
				libmaus2::aio::InputStreamInstance::unique_ptr_type pISI(new libmaus2::aio::InputStreamInstance(Vfiltertmpfn[i]));
				std::istream & ISI = *pISI;
				libmaus2::aio::SerialisedPeeker<Repeat> SP(ISI);
				Repeat R;
				while ( SP.getNext(R) )
				{
					R.serialise(OSI);
				}
				pISI.reset();
				libmaus2::aio::FileRemoval::removeFile(Vfiltertmpfn[i]);
			}
		}
		std::cerr << "[V] concatenated" << std::endl;
	}
};

int lasrefinerepeats(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const changefn = outfn + ".change";

	while ( libmaus2::util::GetFileSize::getFileSize(changefn) > sizeof(uint64_t) /* length of serialised empty vector */ )
	{
		std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
		uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
		bool const precise = arg.uniqueArgPresent("p");

		std::string const dbfn = arg[1];

		std::vector<std::string> Vinfn;
		for ( uint64_t i = 2; i < arg.size(); ++i )
			Vinfn.push_back(arg[i]);

		std::ostringstream transfnostr;
		transfnostr << tmpfilebase << "_trans_part";
		std::string const transfn = transfnostr.str();
		libmaus2::util::TempFileRemovalContainer::addTempFile(transfn);
		std::vector<std::string> transtmp(1,transfn);

		#if 1
		std::vector < std::string > subargs;

		subargs.push_back(std::string("computeTranslation"));
		subargs.push_back(std::string("-T") + tmpfilebase);
		{
			std::ostringstream argstr;
			argstr << "-t" << numthreads;
			subargs.push_back(argstr.str());
		}
		if ( precise )
		{
			subargs.push_back("-p");
		}

		subargs.push_back(transfn); // output file
		subargs.push_back(outfn); // input repeats
		subargs.push_back(dbfn);
		for ( uint64_t i = 0; i < Vinfn.size(); ++i )
			subargs.push_back(Vinfn[i]);

		libmaus2::util::ArgParser const subarg(subargs);
		RepeatTranslation::computeTranslation(subarg);
		#endif

		std::string const mergetmp = RepeatTranslation::mergeTranslationParts(tmpfilebase, transtmp);
		RepeatTranslation::filterTranslationParts(mergetmp,numthreads,tmpfilebase);

		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(Vinfn.size());

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < Vinfn.size(); ++i )
		{
			std::ostringstream memindexname;
			memindexname << "mem:index_" << i << ".las.idx";
			std::string const memindexfn = memindexname.str();
			libmaus2::util::GetFileSize::copy(libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(Vinfn[i]),memindexfn);

			libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
				new libmaus2::dazzler::align::DalignerIndexDecoder(Vinfn[i],memindexfn)
			);
			Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
		}

		uint64_t io = 0;
		::libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > IH;
		for ( uint64_t i = 0; i < Vinfn.size(); ++i )
		{
			if (
			! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]))
				||
				libmaus2::util::GetFileSize::isOlder(
					libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]),
					Vinfn[i]
				)
			)
			{
				libmaus2::dazzler::align::OverlapIndexer::constructIndex(Vinfn[i]);
				std::cerr << "[V] constructing index for " << Vinfn[i] << "\n";

			}

			int64_t const mini = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(Vinfn[i]);
			int64_t const maxi = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(Vinfn[i]);

			if ( mini >= 0 )
			{
				assert ( maxi >= mini );
				Vinfn[io] = Vinfn[i];
				IH.push(io,std::pair<uint64_t,uint64_t>(mini,maxi+1));
			}
		}
		Vinfn.resize(io);
		IH.resize(io);
		for ( uint64_t i = 1; i < io; ++i )
		{
			bool const ok = IH[i-1].second <= IH[i].first;
			if ( ! ok )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] LAS files or not given in increasing order of A read intervals" << std::endl;
				lme.finish();
				throw lme;
			}
		}
		for ( uint64_t i = 1; i < io; ++i )
			IH[i-1].second = IH[i].first;

		libmaus2::util::IntervalTree IT(IH,0,IH.size(),true);

		std::cerr << "[V] copying " << dbfn << " to memory...";
		std::vector<std::string> tracklist;
		tracklist.push_back("inqual");
		// libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbfn,"mem:db1prefix",&tracklist));
		libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
			libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbfn,&tracklist)
		);
		std::cerr << "done." << std::endl;
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
		PDB->computeTrimVector();

		libmaus2::dazzler::db::DatabaseFile & DB = *PDB;
		if ( DB.part != 0 )
		{
			std::cerr << "Partial databases are not supported." << std::endl;
			return EXIT_FAILURE;
		}

		// get read length vector
		std::vector<uint64_t> RL;
		DB.getAllReadLengths(RL);

		// int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vinfn);

		// load change vector
		std::vector<uint64_t> Vchange;
		{
			libmaus2::aio::InputStreamInstance ISI(changefn);
			Vchange = libmaus2::util::NumberSerialisation::deserialiseNumberVector<uint64_t>(ISI);
		}

		// load repeats
		libmaus2::aio::InputStreamInstance ISI(outfn);
		std::vector<Repeat> VR;
		while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
		{
			Repeat R(ISI);
			VR.push_back(R);
		}

		// sort repeats by id,startpos
		std::sort(VR.begin(),VR.end(),RepeatPosComparator());

		#if 0
		RepeatTranslation::computeTranslation(tmpfilebase,numthreads,Vchange,VR,DB,IT,Adalindex,Vinfn,precise,RL,transfn);
		#endif


		// erase change vector
		bool changed = false;
		Vchange.resize(0);

		{
			std::vector < Repeat > BR;

			filterRepeats(mergetmp);

			// newly created repeats
			std::vector < Repeat > VRN;

			libmaus2::aio::SerialisedPeeker<Repeat> SP(mergetmp);
			Repeat R;
			while ( SP.peekNext(R) )
			{
				uint64_t const id = R.id;
				std::vector < Repeat > LVR;
				while ( SP.peekNext(R) && R.id == id )
				{
					SP.getNext(R);
					LVR.push_back(R);
				}

				typedef std::vector<Repeat>::const_iterator it;
				std::pair<it,it> IP = std::equal_range(VR.begin(),VR.end(),Repeat(id),RepeatIdComparator());

				std::vector < libmaus2::math::IntegerInterval<int64_t> > IVA;
				std::vector < libmaus2::math::IntegerInterval<int64_t> > IVB;
				for ( it i = IP.first; i != IP.second; ++i )
				{
					IVA.push_back(libmaus2::math::IntegerInterval<int64_t>(i->abpos,i->aepos-1));
					IVB.push_back(IVA.back());
				}
				for ( uint64_t i = 0; i < LVR.size(); ++i )
				{
					IVB.push_back(libmaus2::math::IntegerInterval<int64_t>(LVR[i].abpos,LVR[i].aepos-1));
				}

				std::vector < libmaus2::math::IntegerInterval<int64_t> > IVM =
					libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(IVB);

				// remove any repat intervals shorter than minlen
				int64_t const minlen = 500;
				uint64_t io = 0;
				for ( uint64_t i = 0; i < IVM.size(); ++i )
					if ( IVM[i].diameter() >= minlen )
						IVM[io++] = IVM[i];
				IVM.resize(io);

				if ( IVM != IVA )
				{
					bool extend = false;
					for ( uint64_t i = 0; i < IVM.size(); ++i )
					{
						std::vector< libmaus2::math::IntegerInterval<int64_t> > VD =
							libmaus2::math::IntegerInterval<int64_t>::difference(IVM[i],IVA);

						uint64_t diam = 0;
						for ( uint64_t j = 0; j < VD.size(); ++j )
							diam += VD[j].diameter();

						if ( diam >= minlen )
							extend = true;
					}

					if ( extend )
					{
						#if 0
						std::cerr << "list for read " << id << " changed from ";
						for ( uint64_t i = 0; i < IVA.size(); ++i )
							std::cerr << IVA[i];
						std::cerr << " to ";
						for ( uint64_t i = 0; i < IVM.size(); ++i )
							std::cerr << IVM[i];
						std::cerr << std::endl;
						#endif

						#if 0
						for ( uint64_t i = 0; i < IVM.size(); ++i )
						{
							std::vector< libmaus2::math::IntegerInterval<int64_t> > VD =
								libmaus2::math::IntegerInterval<int64_t>::difference(IVM[i],IVA);

							uint64_t diam = 0;
							for ( uint64_t j = 0; j < VD.size(); ++j )
								diam += VD[j].diameter();

							if ( diam >= minlen )
							{
								std::cerr << "\tnew diam " << diam << " for " << IVM[i] << std::endl;
								for ( uint64_t j = 0; j < VD.size(); ++j )
									std::cerr << "\t\t" << VD[j] << std::endl;
							}
						}
						#endif

						for ( uint64_t i = 0; i < IVM.size(); ++i )
							VRN.push_back(Repeat(id,IVM[i].from,IVM[i].to+1));

						changed = true;
						Vchange.push_back(id);
					}
				}
			}

			if ( changed )
			{
				// merge to vector VRM
				std::vector < Repeat > VRM;

				uint64_t ilow = 0, jlow = 0;

				while ( ilow < VR.size() && jlow < VRN.size() )
				{
					if ( VR[ilow].id < VRN[jlow].id )
					{
						uint64_t const id = VR[ilow].id;

						while ( ilow < VR.size() && VR[ilow].id == id )
							VRM.push_back(VR[ilow++]);
					}
					else if ( VRN[jlow].id < VR[ilow].id )
					{
						uint64_t const id = VRN[jlow].id;

						while ( jlow < VRN.size() && VRN[jlow].id == id )
							VRM.push_back(VRN[jlow++]);
					}
					else
					{
						uint64_t const id = VRN[jlow].id;

						while ( jlow < VRN.size() && VRN[jlow].id == id )
							VRM.push_back(VRN[jlow++]);

						while ( ilow < VR.size() && VR[ilow].id == id )
							ilow++;
					}
				}

				while ( ilow < VR.size() )
					VRM.push_back(VR[ilow++]);
				while ( jlow < VRN.size() )
					VRM.push_back(VRN[jlow++]);

				libmaus2::aio::OutputStreamInstance OSI(outfn);
				for ( uint64_t i = 0; i < VRM.size(); ++i )
					VRM[i].serialise(OSI);
				OSI.flush();
			}
		}

		libmaus2::aio::FileRemoval::removeFile(mergetmp);
		libmaus2::aio::FileRemoval::removeFile(mergetmp + ".index");

		std::cerr << "[V] Vchange.size()=" << Vchange.size() << std::endl;

		{
			libmaus2::aio::OutputStreamInstance::unique_ptr_type COSI(new libmaus2::aio::OutputStreamInstance(changefn));
			libmaus2::util::NumberSerialisation::serialiseNumberVector(*COSI,Vchange);
			COSI->flush();
			COSI.reset();
		}
	}

	return EXIT_SUCCESS;
}

int lasdetectsimplerepeats(libmaus2::util::ArgParser const & arg)
{
	double const termval = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultTermVal();
	uint64_t const dthres = arg.uniqueArgPresent("d") ? arg.getParsedArg<uint64_t>("d") : getDefaultDepthThres();
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

	std::string const tmpfn = tmpfilebase + ".reptmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmpfn);

	std::vector<std::string> Vtmpfn;
	libmaus2::autoarray::AutoArray<libmaus2::aio::OutputStreamInstance::unique_ptr_type> AtmpOSI(numthreads);
	libmaus2::autoarray::AutoArray<Handler::unique_ptr_type> Ahandler(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		std::ostringstream ostr;
		ostr << tmpfilebase << ".reptmp." << i;
		Vtmpfn.push_back(ostr.str());

		libmaus2::util::TempFileRemovalContainer::addTempFile(Vtmpfn.back());

		libmaus2::aio::OutputStreamInstance::unique_ptr_type tptr(
			new libmaus2::aio::OutputStreamInstance(Vtmpfn.back())
		);

		AtmpOSI[i] = UNIQUE_PTR_MOVE(tptr);

		Handler::unique_ptr_type hptr(new Handler());
		Ahandler[i] = UNIQUE_PTR_MOVE(hptr);
	}

	std::string const outfn = arg[0];
	std::string const dbfn = arg[1];

	std::vector<std::string> Vinfn;
	for ( uint64_t i = 2; i < arg.size(); ++i )
		Vinfn.push_back(arg[i]);

	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(Vinfn.size());

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		if (

			! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]))
			||
			libmaus2::util::GetFileSize::isOlder(libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]),Vinfn[i])
			||
			! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(Vinfn[i]))
			||
			libmaus2::util::GetFileSize::isOlder(libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(Vinfn[i]),Vinfn[i])
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(Vinfn[i]);
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "[V] constructed index for " << Vinfn[i] << "\n";

		}

		std::ostringstream memindexname;
		memindexname << "mem:index_" << i << ".las.idx";
		std::string const memindexfn = memindexname.str();
		libmaus2::util::GetFileSize::copy(libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(Vinfn[i]),memindexfn);

		libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
			new libmaus2::dazzler::align::DalignerIndexDecoder(Vinfn[i],memindexfn)
		);
		Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
	}

	uint64_t io = 0;
	::libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > IH;
	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		if (
		! libmaus2::util::GetFileSize::fileExists(libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]))
			||
			libmaus2::util::GetFileSize::isOlder(
				libmaus2::dazzler::align::OverlapIndexer::getIndexName(Vinfn[i]),
				Vinfn[i]
			)
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(Vinfn[i]);
			std::cerr << "[V] constructing index for " << Vinfn[i] << "\n";

		}

		int64_t const mini = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(Vinfn[i]);
		int64_t const maxi = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(Vinfn[i]);

		if ( mini >= 0 )
		{
			assert ( maxi >= mini );
			Vinfn[io] = Vinfn[i];
			IH.push(io,std::pair<uint64_t,uint64_t>(mini,maxi+1));
		}
	}
	Vinfn.resize(io);
	IH.resize(io);
	for ( uint64_t i = 1; i < io; ++i )
	{
		bool const ok = IH[i-1].second <= IH[i].first;
		if ( ! ok )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] LAS files or not given in increasing order of A read intervals" << std::endl;
			lme.finish();
			throw lme;
		}
	}
	for ( uint64_t i = 1; i < io; ++i )
		IH[i-1].second = IH[i].first;

	libmaus2::util::IntervalTree IT(IH,0,IH.size(),true);

	std::cerr << "[V] copying " << dbfn << " to memory...";
	std::vector<std::string> tracklist;
	tracklist.push_back("inqual");
	// libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbfn,"mem:db1prefix",&tracklist));
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbfn,&tracklist)
	);
	std::cerr << "done." << std::endl;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
	PDB->computeTrimVector();

	libmaus2::dazzler::db::DatabaseFile & DB = *PDB;
	if ( DB.part != 0 )
	{
		std::cerr << "Partial databases are not supported." << std::endl;
		return EXIT_FAILURE;
	}

	// get read length vector
	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	// load inqual track
	libmaus2::dazzler::db::Track::unique_ptr_type Ptrack(DB.readTrack("inqual",0));
	// overlap checker class
	libmaus2::dazzler::align::OverlapProperCheck const OPC(RL,*Ptrack,termval);

	// output stream
	libmaus2::aio::OutputStreamInstance::unique_ptr_type pOSI(new libmaus2::aio::OutputStreamInstance(tmpfn));
	// std::ostream & out = *pOSI;

	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		std::string const infn = Vinfn[i];
		libmaus2::dazzler::align::DalignerIndexDecoder & dalindex = *Adalindex[i];

		int64_t const mini = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(infn);
		int64_t const maxi = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(infn);

		// skip empty file
		if ( mini < 0 )
			continue;

		uint64_t const range = maxi-mini+1;
		uint64_t const low = mini;
		uint64_t const high = maxi+1;

		uint64_t const readsperthread = (range + numthreads - 1)/numthreads;
		uint64_t const numblocks = (range + readsperthread - 1)/readsperthread;

		int volatile tfailed = 0;
		libmaus2::parallel::PosixSpinLock tfailedlock;

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t t = 0; t < numblocks; ++t )
		{
			try
			{
				#if defined(_OPENMP)
				uint64_t const tid = omp_get_thread_num();
				#else
				uint64_t const tid = 0;
				#endif

				// repeat handler
				Handler & H = *(Ahandler[tid]);

				uint64_t const tlow = low + t * readsperthread;
				uint64_t const thigh = std::min(tlow + readsperthread,high);

				// open LAS file
				libmaus2::dazzler::align::AlignmentFileDecoder::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileAt(infn,tlow,thigh,dalindex));
				// get tspace parameter
				int64_t const tspace = Plas->tspace;

				std::ostream & out = *(AtmpOSI[tid]);

				// overlap
				libmaus2::dazzler::align::Overlap OVL;

				// repeat evidence vector
				std::vector< std::pair<int64_t,int64_t> > V;

				// previous read id
				int64_t prevread = -1;

				// read overlaps from LAS file
				for ( uint64_t c = 0 ; Plas->getNextOverlap(OVL) ; ++c )
				{
					// if this is a new A read id
					if ( OVL.aread != prevread )
					{
						// handle data for previous A read
						H.handleVector(out,prevread,V,dthres);
						if ( (prevread != -1) && (prevread % 1000 == 0) )
						{
							libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
							std::cerr << "[V] " << prevread << "/[" << tlow << "," << thigh << ")" << std::endl;
						}
						prevread = OVL.aread;
					}

					// get if this is a proper overlap
					libmaus2::dazzler::align::OverlapProperCheck::OverlapProperCheckInfo const proper = OPC(OVL,tspace);

					// if improper on left and improper on right
					if ( ! proper.termleft && ! proper.termright )
					{
						// std::cerr << "repeat " << OVL.aread << " " << OVL.path.abpos << "," << OVL.path.aepos << std::endl;

						// put in evidence vector
						V.push_back(std::pair<uint64_t,uint64_t>(OVL.path.abpos,OVL.path.aepos));
					}
				}

				// handle last read
				H.handleVector(out,prevread,V,dthres);

				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << "[V] finished [" << tlow << "," << thigh << ")" << std::endl;
				}
			}
			catch(std::exception const & ex)
			{
				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << "[E] " << ex.what() << std::endl;
				}
				{
					libmaus2::parallel::ScopePosixSpinLock slock(tfailedlock);
					tfailed = 1;
				}
			}
		}

		if ( tfailed )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] parallel processing loop failed" << std::endl;
			lme.finish();
			throw lme;
		}
	}

	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		AtmpOSI[i]->flush();
		AtmpOSI[i].reset();

		{
			libmaus2::aio::InputStreamInstance ISI(Vtmpfn[i]);
			while ( ISI.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::util::NumberSerialisation::serialiseNumber(*pOSI,
					libmaus2::util::NumberSerialisation::deserialiseNumber(ISI)
				);
			}
		}

		{
			libmaus2::aio::FileRemoval::removeFile(Vtmpfn[i]);
		}
	}

	pOSI->flush();
	pOSI.reset();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<Repeat,RepeatPosComparator>::sort(tmpfn,128*1024*1024 /* blocksize */);

	libmaus2::aio::OutputStreamFactoryContainer::rename(tmpfn,outfn);

	std::string const changefn = outfn + ".change";
	{
		libmaus2::aio::OutputStreamInstance::unique_ptr_type COSI(new libmaus2::aio::OutputStreamInstance(changefn));
		std::vector<uint64_t> Vchange;
		for ( uint64_t i = 0; i < RL.size(); ++i )
			Vchange.push_back(i);
		libmaus2::util::NumberSerialisation::serialiseNumberVector(*COSI,Vchange);
		COSI->flush();
		COSI.reset();
	}

	return EXIT_SUCCESS;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.rep in.db in.las\n";
			std::cerr << std::endl;
			std::cerr << "optional parameters:" << std::endl << std::endl;
			std::cerr << " -e: error threshold for proper alignment termination (default: " << getDefaultTermVal() << ")" << std::endl;
			std::cerr << " -d: depth threshold for repeat detection (default: " << getDefaultDepthThres() << ")" << std::endl;
			return EXIT_SUCCESS;
		}
		else
		{
			int r = lasdetectsimplerepeats(arg);
			if ( r == EXIT_SUCCESS )
				r = lasrefinerepeats(arg);
			return r;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
