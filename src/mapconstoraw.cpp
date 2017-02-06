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
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/autoarray/AutoArray.hpp>
#include <libmaus2/math/lowbits.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <queue>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/dazzler/align/TracePoint.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/lcs/NPL.hpp>
#include <libmaus2/parallel/SimpleThreadWorkPackage.hpp>
#include <libmaus2/parallel/SimpleThreadPool.hpp>
#include <libmaus2/parallel/SimpleThreadWorkPackageDispatcher.hpp>
#include <libmaus2/parallel/SimpleThreadPoolWorkPackageFreeList.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/StrCmpNum.hpp>
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/lcs/SuffixArrayLCS.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/parallel/LockedFreeList.hpp>
#include <libmaus2/parallel/LockedBool.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>

struct Match
{
	int64_t off;
	int64_t l;

	Match()
	{

	}
	Match(int64_t const roff, int64_t const rl)
	: off(roff), l(rl) {}

	bool operator<(Match const & M) const
	{
		return off < M.off;
	}

	int64_t getFrom() const
	{
		return off;
	}

	int64_t getTo() const
	{
		return off + l;
	}

	libmaus2::math::IntegerInterval<int64_t> getInterval() const
	{
		return libmaus2::math::IntegerInterval<int64_t>(getFrom(),getTo()-1);
	}
};



static std::string getBase(std::string const & id)
{
	assert ( id.find('_') != std::string::npos );
	return id.substr(0,id.find_last_of('_'));
}

static std::string tolow(std::string s)
{
	for ( uint64_t i = 0; i < s.size(); ++i )
		s[i] = ::tolower(s[i]);
	return s;
}

struct FragmentPatch
{
	std::string base;
	std::pair<uint64_t,uint64_t> P;
	std::string S;

	FragmentPatch()
	{
	}

	FragmentPatch(std::string const & rbase, std::pair<uint64_t,uint64_t> const & rP, std::string const & rS) : base(rbase), P(rP), S(rS) {}

	bool operator<(FragmentPatch const & O) const
	{
		if ( P.first != O.P.first )
			return P.first < O.P.first;
		else
			return P.second > O.P.second;
	}

	uint64_t size() const
	{
		return P.second-P.first;
	}

	bool inter(FragmentPatch const & O) const
	{
		libmaus2::math::IntegerInterval<int64_t> IA(P.first,P.second-1);
		libmaus2::math::IntegerInterval<int64_t> IB(O.P.first,O.P.second-1);

		return ! (IA.intersection(IB).isEmpty());
	}
};

struct PatchSizeComparator
{
	bool operator()(FragmentPatch const & A, FragmentPatch const & B) const
	{
		return A.size() > B.size();
	}
};

std::vector < FragmentPatch > filterPatches(std::vector < FragmentPatch > patches)
{
	std::sort(patches.begin(),patches.end());

	bool ok = true;
	for ( uint64_t i = 1; i < patches.size(); ++i )
		ok = ok && ( patches[i].P.first >= patches[i-1].P.second );

	if ( ok )
		return patches;

	// order by decreasing diameter
	std::sort(patches.begin(),patches.end(),PatchSizeComparator());

	// greedily use biggest patches until nothing fits anymore
	std::vector < FragmentPatch > fpatches;
	for ( uint64_t i = 0; i < patches.size(); ++i )
	{
		bool ok = true;
		for ( uint64_t j = 0; j < fpatches.size(); ++j )
			if ( fpatches[j].inter(patches[i]) )
				ok = false;

		if ( ok )
			fpatches.push_back(patches[i]);
	}

	std::sort(fpatches.begin(),fpatches.end());

	return fpatches;
}

void printStringLines(std::ostream & out, std::string const & s, uint64_t const cols)
{
	char const * ita = s.c_str();
	char const * itc = ita;
	char const * ite = ita + s.size();

	while ( itc != ite )
	{
		uint64_t const rest = ite-itc;
		uint64_t const tocopy = std::min(cols,rest);

		out.write(itc,tocopy);
		itc += tocopy;
		out.put('\n');
	}
}

void handle(std::vector<FragmentPatch> & patches, std::string & ref, std::ostream & out)
{
	if ( patches.size() )
	{
		patches = filterPatches(patches);

		for ( uint64_t i = 1; i < patches.size(); ++i )
			assert ( patches[i].P.first >= patches[i-1].P.second );

		std::ostringstream ostr;
		uint64_t l = 0;

		out << '>' << patches[0].base << " ";

		for ( uint64_t i = 0; i < patches.size(); ++i )
		{
			assert ( l <= patches[i].P.first );

			out << "[" << patches[i].P.first << "," << patches[i].P.second << ")";

			if ( l < patches[i].P.first )
				ostr << tolow(ref.substr(l,patches[i].P.first-l));

			ostr << patches[i].S;

			l = patches[i].P.second;
		}

		if ( l < ref.size() )
			ostr << tolow(ref.substr(l));

		out << "\n";
		// out << ostr.str() << "\n";
		printStringLines(out,ostr.str(),80);
	}

	patches.resize(0);
}


typedef std::vector<FragmentPatch> patch_vec;
typedef libmaus2::util::shared_ptr<patch_vec>::type patch_vec_ptr;

struct PatchVectorTypeInfo
{
	typedef patch_vec element_type;
	typedef patch_vec_ptr pointer_type;

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}

	static pointer_type getNullPointer()
	{
		pointer_type p;
		return p;
	}
};

struct PatchVectorAllocator
{
	typedef patch_vec element_type;

	patch_vec_ptr operator()() const
	{
		patch_vec_ptr V(new patch_vec);
		return V;
	}
};


struct ReadWorkPackage : public libmaus2::parallel::SimpleThreadWorkPackage
{
	typedef ReadWorkPackage this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	ReadWorkPackage() : SimpleThreadWorkPackage() {}
	ReadWorkPackage(uint64_t const rpriority, uint64_t const rdispatcherid) : SimpleThreadWorkPackage(rpriority,rdispatcherid) {}

	char const * getPackageName() const
	{
		return "ReadWorkPackage";
	}
};

struct ProcessWorkPackage : public libmaus2::parallel::SimpleThreadWorkPackage
{
	typedef ProcessWorkPackage this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string ref;
	patch_vec_ptr V;

	ProcessWorkPackage() : SimpleThreadWorkPackage() {}
	ProcessWorkPackage(uint64_t const rpriority, uint64_t const rdispatcherid) : SimpleThreadWorkPackage(rpriority,rdispatcherid) {}

	char const * getPackageName() const
	{
		return "ProcessWorkPackage";
	}
};

struct ReadWorkPackageReturnInterface
{
	virtual void put(ReadWorkPackage * P) = 0;
};

struct ReadDoneInterface
{
	virtual void readDone() = 0;
};

struct ReturnPatchVectorInterface
{
	virtual void returnPatchVectorProcessed(patch_vec_ptr P) = 0;
	// virtual void returnPatchVectorUnprocessed(patch_vec_ptr P) = 0;
};

struct EnquePatchVectorInterface
{
	virtual void enquePatchVector(std::string const & ref, patch_vec_ptr P) = 0;
};

struct OutputEntry
{
	std::string out;
	std::string err;
	std::string id;

	OutputEntry()
	{
	}
	OutputEntry(std::string const & rout, std::string const & rerr, std::string const rid)
	: out(rout), err(rerr), id(rid)
	{
	}
	OutputEntry(std::istream & in)
	:
		out(libmaus2::util::StringSerialisation::deserialiseString(in)),
		err(libmaus2::util::StringSerialisation::deserialiseString(in)),
		id(libmaus2::util::StringSerialisation::deserialiseString(in))
	{

	}

	std::istream & deserialise(std::istream & in)
	{
		*this = OutputEntry(in);
		return in;
	}

	std::ostream & serialise(std::ostream & ostr) const
	{
		libmaus2::util::StringSerialisation::serialiseString(ostr,out);
		libmaus2::util::StringSerialisation::serialiseString(ostr,err);
		libmaus2::util::StringSerialisation::serialiseString(ostr,id);
		return ostr;
	}

	bool operator<(OutputEntry const & O) const
	{
		int const r = libmaus2::bambam::StrCmpNum::strcmpnum(id.c_str(),O.id.c_str());
		return (r < 0);
	}
};

struct ReadWorkPackageDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	ReadWorkPackageReturnInterface & retintf;
	ReadDoneInterface & readDone;
	ReturnPatchVectorInterface & retPatchIntf;
	EnquePatchVectorInterface & enqPatchInf;

	libmaus2::bambam::StrCmpNum comp;

	std::string const fullfn;
	std::string const fragmentfn;

	libmaus2::fastx::FastaPeeker infull;
	libmaus2::fastx::FastaPeeker infrag;

	libmaus2::parallel::PosixSpinLock lock;

	libmaus2::parallel::LockedFreeList< patch_vec, PatchVectorAllocator, PatchVectorTypeInfo > & vecFreeList;

	ReadWorkPackageDispatcher(
		ReadWorkPackageReturnInterface & rretintf,
		ReadDoneInterface & rreadDone,
		ReturnPatchVectorInterface & rretPatchIntf,
		EnquePatchVectorInterface & renqPatchInf,
		std::string const & rfullfn,
		std::string const & rfragmentfn,
		libmaus2::parallel::LockedFreeList< patch_vec, PatchVectorAllocator, PatchVectorTypeInfo > & rvecFreeList
	) : retintf(rretintf), readDone(rreadDone), retPatchIntf(rretPatchIntf), enqPatchInf(renqPatchInf), fullfn(rfullfn), fragmentfn(rfragmentfn),
	    infull(fullfn), infrag(fragmentfn), vecFreeList(rvecFreeList)
	{

	}

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		ReadWorkPackage * rwp = dynamic_cast<ReadWorkPackage *>(P);

		bool readdone = false;

		if ( lock.trylock() )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(lock,true /* prelocked */);

			patch_vec_ptr V;

			while (
				(!readdone)
				&&
				(V = vecFreeList.getIf())
			)
			{
				V->resize(0);

				std::string ref;

				libmaus2::fastx::FastAReader::pattern_type inpat;
				libmaus2::fastx::FastAReader::pattern_type fragpat;

				while ( infull.peekNext(inpat) && infrag.peekNext(fragpat) )
				{
					std::string const lbase = getBase(fragpat.getShortStringId());

					if ( V->size() && lbase != V->front().base )
						break;

					int const r = comp.strcmpnum(inpat.sid.c_str(),lbase.c_str());

					if ( r < 0 )
					{
						infull.getNext(inpat);
					}
					else if ( r > 0 )
					{
						infrag.getNext(fragpat);
					}
					else
					{
						assert ( lbase == inpat.sid );
						infrag.getNext(fragpat);

						if ( ! V->size() )
							ref = inpat.spattern;
						V->push_back(FragmentPatch(lbase,std::pair<uint64_t,uint64_t>(),fragpat.spattern));
					}
				}

				if ( V->size() )
				{
					enqPatchInf.enquePatchVector(ref,V);
				}
				else
				{
					vecFreeList.put(V);
					readdone = true;
				}
			}
		}

		retintf.put(rwp);

		if ( readdone )
			readDone.readDone();
	}
};

struct ProcessWorkPackageReturnInterface
{
	virtual void put(ProcessWorkPackage * P) = 0;
};

struct ProcessWorkPackageDispatcher : public libmaus2::parallel::SimpleThreadWorkPackageDispatcher
{
	ProcessWorkPackageReturnInterface & retintf;
	ReturnPatchVectorInterface & retvecintf;
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<OutputEntry> & SBO;
	libmaus2::parallel::PosixSpinLock SBOlock;

	int w_err;
	int w_back;
	int64_t maxclip = 50;

	struct Kmer
	{
		uint64_t v;
		uint64_t p;

		Kmer() {}
		Kmer(uint64_t const rv, uint64_t const rp) : v(rv), p(rp) {}

		bool operator<(Kmer const & K) const
		{
			if ( v != K.v )
				return v < K.v;
			else
				return p < K.p;
		}
	};

	static uint64_t getKmers(std::string const & s, uint64_t const k, libmaus2::autoarray::AutoArray<Kmer> & K)
	{
		uint64_t const n = s.size();
		uint64_t const numk = (n >= k) ? (n-k+1) : 0;
		uint64_t const kmask = libmaus2::math::lowbits(2*(k-1));
		assert ( k );

		uint64_t v = 0;
		for ( uint64_t i = 0; i < k-1; ++i )
		{
			v <<= 2;
			v |= libmaus2::fastx::mapChar(s[i]);
		}

		uint64_t o = 0;
		for ( uint64_t i = 0; i < numk ; ++i )
		{
			v &= kmask;
			v <<= 2;
			v |= libmaus2::fastx::mapChar(s[i+k-1]);
			K.push(o,Kmer(v,i));
		}

		std::sort(K.begin(),K.begin()+o);

		return o;
	}

	struct HeapTodo
	{
		Kmer const * ka;
		Kmer const * kb_c;
		Kmer const * kb_e;

		HeapTodo()
		{

		}

		HeapTodo(
			Kmer const * rka,
			Kmer const * rkb_c,
			Kmer const * rkb_e)
		: ka(rka), kb_c(rkb_c), kb_e(rkb_e)
		{

		}

		bool hasNext() const
		{
			return kb_c+1 != kb_e;
		}

		HeapTodo getNext()
		{
			HeapTodo T = *this;
			++T.kb_c;
			return T;
		}

		int64_t getDiag() const
		{
			return static_cast<int64_t>(ka->p) - static_cast<int64_t>(kb_c->p);
		}

		bool operator<(HeapTodo const & O) const
		{
			return getDiag() < O.getDiag();
		}

		bool operator>(HeapTodo const & O) const
		{
			return getDiag() > O.getDiag();
		}
	};

	bool seedAndExtend(
		std::string const & sa, std::string const & sb,
		uint64_t const k, int64_t const tspace,
		int64_t const aid, int64_t const bid,
		uint64_t const minlen,
		libmaus2::dazzler::align::Overlap & ROVL,
		uint64_t const maxtraces
	)
	{
		std::string const ra = sa;
		std::string const rb = sb;

		libmaus2::autoarray::AutoArray<Kmer> KA;
		libmaus2::autoarray::AutoArray<Kmer> KB;

		// compute k-mers in read a and b
		uint64_t const oa = getKmers(ra,k,KA);
		uint64_t const ob = getKmers(rb,k,KB);

		uint64_t ia = 0;
		uint64_t ib = 0;

		std::priority_queue< HeapTodo,std::vector<HeapTodo>,std::greater<HeapTodo> > Q;

		// look for matching k-mers
		while ( ia < oa && ib < ob )
		{
			if ( KA[ia].v < KB[ib].v )
				++ia;
			else if ( KB[ib].v < KA[ia].v )
				++ib;
			else
			{
				assert ( KA[ia].v == KB[ib].v );

				uint64_t ja = ia+1;
				while ( ja < oa && KA[ja].v == KA[ia].v )
					++ja;
				uint64_t jb = ib+1;
				while ( jb < ob && KB[jb].v == KB[ib].v )
					++jb;

				// reverse region in B so positions on B are in decreasing order
				std::reverse(KB.begin()+ib,KB.begin()+jb);

				for ( uint64_t z = ia; z < ja; ++z )
					Q.push(HeapTodo(KA.begin()+z,KB.begin()+ib,KB.begin()+jb));

				ia = ja;
				ib = jb;
			}
		}

		libmaus2::lcs::NNP nnp(w_err,w_back);
		libmaus2::lcs::NNPTraceContainer nnptracecontainer;
		libmaus2::lcs::AlignmentTraceContainer ATC;
		libmaus2::autoarray::AutoArray<Match const *> QR;
		libmaus2::util::SimpleQueue<libmaus2::geometry::RangeSet<Match>::search_q_element> Rtodo;

		std::map< int64_t, libmaus2::geometry::RangeSet<Match>::shared_ptr_type > D;

		#if defined(D_DEBUG)
		std::map< int64_t, std::vector<libmaus2::math::IntegerInterval<int64_t> > > TT;
		#endif

		libmaus2::lcs::NNPAlignResult bestres(0,1,0,1,1);

		int64_t prevdiag = std::numeric_limits<int64_t>::min();

		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPV;
		libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPVold;
		std::set< libmaus2::dazzler::align::TracePoint > TPVseen;
		std::map< uint64_t, libmaus2::dazzler::align::Overlap > OVLseen;

		// process matching k-mers in order of increasing diagonal index
		for ( uint64_t traceid = 0; ! Q.empty(); ++traceid )
		{
			if ( traceid >= maxtraces )
				return false;

			HeapTodo H = Q.top();
			Q.pop();

			if ( H.hasNext() )
				Q.push(H.getNext());

			int64_t const nd = H.getDiag();

			if ( traceid && traceid % (1024*16) == 0 )
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] trace " << traceid << " diag " << nd << " " << sa.size() << " " << sb.size() << std::endl;
			}

			assert ( nd >= prevdiag );
			prevdiag = nd;

			// remove diagonals which we no longer need
			while ( D.begin() != D.end() && D.begin()->first < nd )
			{
				// std::cerr << "removing " << D.begin()->first << std::endl;
				D.erase(D.begin());
			}
			// std::cerr << nd << " " << H.ka->p << " " << H.kb_c->p << std::endl;

			// check whether seed was already processed (was used on a previous alignment)
			if ( D.find(nd) != D.end() && D.find(nd)->second->search(Match(std::min(H.ka->p,H.kb_c->p),k),Rtodo) )
				continue;

			// compute alignment
			libmaus2::lcs::NNPAlignResult res;
			if ( aid == bid )
			{
				res = nnp.align(
					ra.begin(),ra.end(),H.ka->p,
					ra.begin(),ra.end(),H.kb_c->p,
					nnptracecontainer,
					aid==bid,
					nnp.getDefaultMinDiag(),
					nnp.getDefaultMaxDiag(),
					true /* runsuffixpositive */
				);
			}
			else
			{
				res = nnp.align(
					ra.begin(),ra.end(),H.ka->p,
					rb.begin(),rb.end(),H.kb_c->p,
					nnptracecontainer,
					aid==bid,
					nnp.getDefaultMinDiag(),
					nnp.getDefaultMaxDiag(),
					true /* runsuffixpositive */
				);
			}

			// compute dense trace
			nnptracecontainer.computeTrace(ATC);

			libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = ATC.ta;
			libmaus2::lcs::AlignmentTraceContainer::step_type const * const te = ATC.te;

			assert ( ta <= te );

			// register matches so we can avoid starting from a seed which is already covered by an alignment
			int64_t apos = res.abpos;
			int64_t bpos = res.bbpos;

			uint64_t matchcount = 0;
			for ( ; ta != te; ++ta )
			{
				switch ( *ta )
				{
					case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
					case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
						if ( matchcount )
						{
							int64_t const d = apos - bpos;
							int64_t const off = std::min(apos,bpos)-matchcount;

							// std::cerr << "M " << matchcount << " " << d << " " << off << std::endl;

							if ( D.find(d) == D.end() )
							{
								libmaus2::geometry::RangeSet<Match>::shared_ptr_type R(
									new libmaus2::geometry::RangeSet<Match>(
										std::min(ra.size()+k,rb.size()+k)
									)
								);
								D[d] = R;
							}

							D.find(d)->second->insert(Match(off,matchcount));

							#if defined(D_DEBUG)
							TT[d].push_back(
								libmaus2::math::IntegerInterval<int64_t>(off,off+matchcount-1)
							);
							#endif

							matchcount = 0;
						}
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
						++matchcount;
						break;
					default:
						break;
				}
				switch ( *ta )
				{
					case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
						bpos++;
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
						apos++;
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
						apos++;
						bpos++;
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
						apos++;
						bpos++;
						break;
					default:
						break;
				}
			}

			if ( matchcount )
			{
				int64_t const d = apos - bpos;
				int64_t const off = std::min(apos,bpos)-matchcount;

				// std::cerr << "M " << matchcount << " " << d << " " << off << std::endl;

				if ( D.find(d) == D.end() )
				{
					libmaus2::geometry::RangeSet<Match>::shared_ptr_type R(
						new libmaus2::geometry::RangeSet<Match>(
							std::min(ra.size()+k,rb.size()+k)
						)
					);
					D[d] = R;
				}

				D.find(d)->second->insert(Match(off,matchcount));

				#if defined(D_DEBUG)
				TT[d].push_back(
					libmaus2::math::IntegerInterval<int64_t>(off,off+matchcount-1)
				);
				#endif

				matchcount = 0;
			}

			assert ( apos == static_cast<int64_t>(res.aepos) );
			assert ( bpos == static_cast<int64_t>(res.bepos) );

			// if match is sufficiently long
			if ( res.aepos - res.abpos >= minlen )
			{
				//std::cerr << res << std::endl;

				if ( res.getErrorRate() < bestres.getErrorRate() )
					bestres = res;

				// compute dazzler style overlap data structure
				libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
					0 /* flags */,
					aid,
					bid,
					res.abpos,
					res.aepos,
					res.bbpos,
					res.bepos,
					tspace,
					ATC
				);

				// get dazzler trace points
				uint64_t const tpvo = OVL.getTracePoints(tspace,traceid,TPV,0);
				libmaus2::math::IntegerInterval<int64_t> IA(res.abpos,res.aepos-1);

				bool dup = false;

				// check whether this or a previous alignment is a duplicate
				for ( uint64_t i = 0; (!dup) && i < tpvo; ++i )
				{
					// look for trace point
					std::set<libmaus2::dazzler::align::TracePoint>::const_iterator it = TPVseen.lower_bound(
						libmaus2::dazzler::align::TracePoint(TPV[i].apos,TPV[i].bpos,0)
					);
					std::vector<uint64_t> killlist;

					for (
						;
						(!dup)
						&&
						it != TPVseen.end()
						&&
						it->apos == TPV[i].apos
						&&
						it->bpos == TPV[i].bpos
						;
						++it
					)
					{
						uint64_t const oldtraceid = it->id;

						assert ( OVLseen.find(oldtraceid) != OVLseen.end() );

						libmaus2::dazzler::align::Overlap const & OVLold = OVLseen.find(oldtraceid)->second;

						libmaus2::math::IntegerInterval<int64_t> IO(OVLold.path.abpos,OVLold.path.aepos-1);
						libmaus2::math::IntegerInterval<int64_t> IC = IA.intersection(IO);

						if ( IA.diameter() <= IO.diameter() && IC.diameter() >= 0.95 * IA.diameter() )
						{
							dup = true;

							#if 0
							libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
							std::cerr << "dup ?\n";
							std::cerr << OVL << std::endl;
							std::cerr << OVLold << std::endl;
							#endif
						}
						else if ( IO.diameter() <= IA.diameter() && IC.diameter() >= 0.95 * IO.diameter() )
						{
							killlist.push_back(oldtraceid);

							#if 0
							libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
							std::cerr << "dup ?\n";
							std::cerr << OVLold << std::endl;
							std::cerr << OVL << std::endl;
							#endif
						}
					}

					for ( uint64_t i = 0; i < killlist.size(); ++i )
					{
						uint64_t const oldtraceid = killlist[i];
						libmaus2::dazzler::align::Overlap const & OVLold = OVLseen.find(oldtraceid)->second;

						uint64_t const tpvoo = OVLold.getTracePoints(tspace,oldtraceid,TPVold,0);
						for ( uint64_t j = 0; j < tpvoo; ++j )
						{
							assert ( TPVseen.find(TPVold[j]) != TPVseen.end() );
							TPVseen.erase(TPVold[j]);
						}

						OVLseen.erase(oldtraceid);
					}

					#if 0
					if ( TPVseen.find(std::pair<int64_t,int64_t>(TPV[i].apos,TPV[i].bpos)) != TPVseen.end() )
					{
						dup = true;
						break;
					}
					#endif
				}

				if ( ! dup )
				{
					for ( uint64_t i = 0; i < tpvo; ++i )
						TPVseen.insert(TPV[i]);
					OVLseen[traceid] = OVL;
				}

				// std::cerr << OVL << std::endl;

				#if 0
				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					std::cerr,
					ra.begin()+res.abpos,res.aepos-res.abpos,
					rb.begin()+res.bbpos,res.bepos-res.bbpos,
					80,
					ATC.ta,ATC.te
				);
				#endif
			}
		}

		int64_t maxspan = -1;
		libmaus2::dazzler::align::Overlap maxovl;

		for ( std::map<uint64_t,libmaus2::dazzler::align::Overlap>::const_iterator ita = OVLseen.begin();
			ita != OVLseen.end(); ++ita )
		{
			libmaus2::dazzler::align::Overlap const & OVL = ita->second;
			// std::cerr << OVL << std::endl;

			int64_t const span = OVL.path.aepos-OVL.path.abpos;
			if ( span > maxspan )
			{
				maxspan = span;
				maxovl = OVL;
			}
		}

		if ( maxspan != -1 )
		{
			ROVL = maxovl;
			return true;
		}
		else
		{
			return false;
		}
	}

	ProcessWorkPackageDispatcher(
		ProcessWorkPackageReturnInterface & rretintf,
		ReturnPatchVectorInterface & rretvecintf,
		libmaus2::sorting::SerialisingSortingBufferedOutputFile<OutputEntry> & rSBO,
		int const rw_err,
		int const rw_back,
		int const rmaxclip
	) : retintf(rretintf), retvecintf(rretvecintf), SBO(rSBO), w_err(rw_err), w_back(rw_back), maxclip(rmaxclip)
	{

	}

	virtual void dispatch(libmaus2::parallel::SimpleThreadWorkPackage * P, libmaus2::parallel::SimpleThreadPoolInterfaceEnqueTermInterface & /* tpi */)
	{
		ProcessWorkPackage * rwp = dynamic_cast<ProcessWorkPackage *>(P);

		#if 0
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "ProcessWorkPackageDispatcher enter" << std::endl;
		}
		#endif

		try
		{
			assert ( rwp->V->size() );

			std::string const base = rwp->V->at(0).base;

			std::ostringstream errstr;

			std::vector < FragmentPatch > patchout;

			// iterate over fragments
			for ( uint64_t i = 0; i < rwp->V->size(); ++i )
			{
				// get fragment
				FragmentPatch const & patch = rwp->V->at(i);

				// get fragment bases
				std::string const firstS = patch.S;

				struct TodoEntry
				{
					uint64_t from;
					uint64_t to;
					uint64_t reffrom;
					uint64_t refto;

					TodoEntry() {}
					TodoEntry(uint64_t const rfrom, uint64_t const rto, uint64_t const rreffrom, uint64_t const rrefto) : from(rfrom), to(rto), reffrom(rreffrom), refto(rrefto) {}
				};

				// push todo entry
				std::deque<TodoEntry> todo;
				todo.push_back(TodoEntry(0,firstS.size(),0,rwp->ref.size()));

				struct OutMatch
				{
					uint64_t abpos;
					uint64_t aepos;
					uint64_t bbpos;
					uint64_t bepos;

					OutMatch() {}
					OutMatch(uint64_t const rabpos, uint64_t const raepos, uint64_t const rbbpos, uint64_t const rbepos)
					: abpos(rabpos), aepos(raepos), bbpos(rbbpos), bepos(rbepos) {}

					bool operator<(OutMatch const & O) const
					{
						return abpos < O.abpos;
					}
				};

				std::vector < OutMatch > VO;

				uint64_t todohandled;
				for ( todohandled = 0 ; todo.size() && todohandled < 8192; ++todohandled )
				{
					TodoEntry T = todo.front();
					todo.pop_front();

					// std::cerr << "todo " << T.from << "," << T.to << " " << T.reffrom << "," << T.refto << std::endl;

					std::string const S = firstS.substr(T.from,T.to-T.from);
					std::string const subref = rwp->ref.substr(T.reffrom,T.refto-T.reffrom);

					libmaus2::dazzler::align::Overlap OVL;
					bool const ok = seedAndExtend(subref,S,16,100,0,1,50 /* minlen */,OVL,8192 /* max traces */);

					if ( ! ok )
					{
						libmaus2::lcs::SuffixArrayLCS::LCSResult const lcsres = libmaus2::lcs::SuffixArrayLCS::lcsmin(subref,S);
						libmaus2::lcs::NNP nnp(w_err,w_back);
						libmaus2::lcs::NNPTraceContainer nnptrace;
						libmaus2::lcs::NNPAlignResult const nnpres = nnp.align(
							subref.begin(),
							subref.end(),
							lcsres.maxpos_a,
							S.begin(),
							S.end(),
							lcsres.maxpos_b,
							nnptrace,
							false /* self */,
							nnp.getDefaultMinDiag(),
							nnp.getDefaultMaxDiag(),
							true /* runsuffixpositive */
						);

						OVL.path.abpos = nnpres.abpos;
						OVL.path.aepos = nnpres.aepos;
						OVL.path.bbpos = nnpres.bbpos;
						OVL.path.bepos = nnpres.bepos;
					}

					OVL.path.abpos += T.reffrom;
					OVL.path.aepos += T.reffrom;

					if ( OVL.path.bepos-OVL.path.bbpos >= maxclip )
					{
						VO.push_back(OutMatch(OVL.path.abpos,OVL.path.aepos,T.from+OVL.path.bbpos,T.from+OVL.path.bepos));
						//std::cerr << "push " << OVL << std::endl;

						if ( OVL.path.bbpos >= maxclip )
						{
							// std::cerr << "reentering " << OVL.path.bbpos << std::endl;
							uint64_t const unaligned = OVL.path.bbpos;
							todo.push_back(TodoEntry(T.from,T.from+unaligned,0,OVL.path.abpos));
						}
						if ( static_cast<int64_t>(S.size()) - OVL.path.bepos >= maxclip )
						{
							// std::cerr << "reentering " << S.size() - OVL.path.bepos << std::endl;
							uint64_t const unaligned = static_cast<int64_t>(S.size()) - OVL.path.bepos;
							todo.push_back(TodoEntry(T.to-unaligned,T.to,OVL.path.aepos,rwp->ref.size()));
						}
					}
				}

				if ( todo.size() )
				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << "[V] bailed out of todo loop" << std::endl;
					if ( rwp->V->size() )
					{
						std::cerr << "[V] base " << patch.base << std::endl;
						std::cerr << "[V] S " << patch.S << std::endl;
						std::cerr << "[V] ref " << rwp->ref << std::endl;
					}
				}

				std::sort(VO.begin(),VO.end());

				bool done = false;
				uint64_t mergehandled = 0;
				for ( mergehandled = 0; (! done) && mergehandled < 8192; ++ mergehandled )
				{
					done = true;

					int64_t remove = -1;

					for ( uint64_t i = 0; done && i < VO.size(); ++i )
						for ( uint64_t j = i+1; done && j < VO.size(); ++j )
						{
							if ( VO[i].aepos <= VO[j].aepos && VO[i].bepos <= VO[j].bbpos )
							{
								#if 0
								std::cerr << "merging "
									<< VO[i].abpos << "," << VO[i].aepos << "," << VO[i].bbpos << "," << VO[i].bepos << " and "
									<< VO[j].abpos << "," << VO[j].aepos << "," << VO[j].bbpos << "," << VO[j].bepos << std::endl;
								#endif

								VO[i].aepos = VO[j].aepos;
								VO[i].bepos = VO[j].bepos;

								remove = j;
								done = false;
							}
						}

					uint64_t o = 0;
					for ( uint64_t i = 0; i < VO.size(); ++i )
						if ( static_cast<int64_t>(i) != remove )
							VO[o++] = VO[i];
					VO.resize(o);
					std::sort(VO.begin(),VO.end());
				}
				if ( ! done )
				{
					std::cerr << "[V] bailed out of merge loop" << std::endl;
				}


				if ( VO.size() && VO.back().bepos != patch.S.size() )
				{
					OutMatch & O = VO.back();

					std::string const subref = rwp->ref.substr(O.abpos);
					std::string const subpatch = firstS.substr(O.bbpos);

					libmaus2::lcs::NPL np;
					np.np(subref.begin(),subref.end(),subpatch.begin(),subpatch.end());

					std::pair<uint64_t,uint64_t> P = np.getTraceContainer().getStringLengthUsed();

					#if 0
					std::cerr << np.traceToString() << std::endl;
					std::cerr << O.bbpos + P.second << " " << patch.S.size() << std::endl;
					#endif

					O.aepos = O.abpos + P.first;
					O.bepos = O.bbpos + P.second;
				}

				if ( VO.size() && VO.front().bbpos != 0 )
				{
					OutMatch & O = VO.front();

					std::string subref = rwp->ref.substr(0,O.aepos);
					std::string subpatch = firstS.substr(0,O.bepos);

					std::reverse(subref.begin(),subref.end());
					std::reverse(subpatch.begin(),subpatch.end());

					libmaus2::lcs::NPL np;
					np.np(subref.begin(),subref.end(),subpatch.begin(),subpatch.end());

					std::pair<uint64_t,uint64_t> P = np.getTraceContainer().getStringLengthUsed();

					//std::cerr << np.traceToString() << std::endl;

					O.abpos = O.aepos - P.first;
					O.bbpos = O.bepos - P.second;

					//std::cerr << O.bepos - P.second << std::endl;
				}

				#if 0
				for ( uint64_t i = 0; i < VO.size(); ++i )
				{
					std::cerr << "merged " << VO[i].abpos << "," << VO[i].aepos << " " << VO[i].bbpos << "," << VO[i].bepos << " " << patch.S.size() << std::endl;
				}
				#endif


				#if 0
				if ( VO.size() )
				{
					if ( VO[0].abpos <= 20 && VO[0].bbpos <= 20 )
					{
						VO[0].abpos = 0;
						VO[0].bbpos = 0;
					}
					if ( rwp->ref.size() - VO.back().aepos <= 20 && firstS.size() - VO.back().bepos <= 20 )
					{
						VO.back().aepos = rwp->ref.size();
						VO.back().bepos = firstS.size();
					}
				}
				#endif

				for ( uint64_t i = 0; i < VO.size(); ++i )
				{
					OutMatch const O = VO[i];

					errstr << patch.base << " [" << O.abpos << "," << O.aepos << ")" << "[" << O.bbpos << "," << O.bepos << ") " << firstS.size() << "\n";

					FragmentPatch outpatch = patch;
					outpatch.P = std::pair<uint64_t,uint64_t>(O.abpos,O.aepos);
					outpatch.S = firstS.substr(O.bbpos,O.bepos-O.bbpos);
					patchout.push_back(outpatch);
				}
			}

			std::ostringstream out;
			handle(patchout, rwp->ref, out);

			#if 0
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::coutlock);
				std::cout << out.str();
			}
			#endif

			#if 0
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << errstr.str();
			}
			#endif

			{
				libmaus2::parallel::ScopePosixSpinLock slock(SBOlock);
				SBO.put(OutputEntry(out.str(),errstr.str(),base));
			}
		}
		catch(std::exception const & ex)
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << ex.what() << std::endl;

			if ( rwp && rwp->V && rwp->V->size() )
				std::cerr << "base " << rwp->V->at(0).base << std::endl;
		}

		#if 0
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "ProcessWorkPackageDispatcher exit" << std::endl;
		}
		#endif

		retvecintf.returnPatchVectorProcessed(rwp->V);
		retintf.put(rwp);
	}
};


struct WorkControl :
	public ReadWorkPackageReturnInterface,
	public ProcessWorkPackageReturnInterface,
	public ReadDoneInterface,
	public ReturnPatchVectorInterface,
	public EnquePatchVectorInterface
{
	libmaus2::parallel::SimpleThreadPoolWorkPackageFreeList<ReadWorkPackage> readPackageFreeList;
	libmaus2::parallel::SimpleThreadPoolWorkPackageFreeList<ProcessWorkPackage> procPackageFreeList;

	uint64_t const numthreads;
	std::string const fullfn;
	std::string const fragmentfn;
	libmaus2::parallel::SimpleThreadPool STP;
	uint64_t nextdispid;
	uint64_t const readdispid;
	uint64_t const procdispid;

	libmaus2::parallel::LockedFreeList< patch_vec, PatchVectorAllocator, PatchVectorTypeInfo > vecFreeList;

	ReadWorkPackageDispatcher RWPD;
	ProcessWorkPackageDispatcher PWPD;

	libmaus2::parallel::LockedBool readDoneFlag;

	uint64_t volatile unfinished;
	libmaus2::parallel::PosixSpinLock unfinishedlock;

	uint64_t volatile finished;
	libmaus2::parallel::PosixSpinLock finishedlock;

	WorkControl(
		uint64_t const rnumthreads, std::string const & rfullfn, std::string const & rfragmentfn,
		libmaus2::sorting::SerialisingSortingBufferedOutputFile<OutputEntry> & SBO,
		int const rw_err,
		int const rw_back,
		int const rmaxclip
	)
	: numthreads(rnumthreads), fullfn(rfullfn), fragmentfn(rfragmentfn), STP(numthreads), nextdispid(0), readdispid(nextdispid++), procdispid(nextdispid++), vecFreeList(2*numthreads),
	  RWPD(*this,*this,*this,*this,fullfn,fragmentfn,vecFreeList), PWPD(*this,*this,SBO,rw_err,rw_back,rmaxclip),
	  readDoneFlag(false), unfinished(0), finished(0)
	{
		STP.registerDispatcher(readdispid,&RWPD);
		STP.registerDispatcher(procdispid,&PWPD);
	}

	void start()
	{
		ReadWorkPackage * P = dynamic_cast<ReadWorkPackage *>(readPackageFreeList.getPackage());
		*P = ReadWorkPackage(0 /* priority */,readdispid);
		STP.enque(P);
	}

	uint64_t getUnfinished()
	{
		uint64_t l;
		unfinishedlock.lock();
		l = unfinished;
		unfinishedlock.unlock();
		return l;
	}

	uint64_t getFinished()
	{
		uint64_t l;
		finishedlock.lock();
		l = finished;
		finishedlock.unlock();
		return l;
	}

	void wait()
	{
		while (
			( ! ( readDoneFlag.get() && (getUnfinished() == 0) ) )
			&&
			(!STP.isInPanicMode())
		)
		{
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] waiting " << readDoneFlag.get() << " " << getUnfinished() << " " << getFinished() << std::endl;
			}
			sleep(1);
		}

		STP.terminate();
		STP.join();
	}

	virtual void put(ReadWorkPackage * P)
	{
		readPackageFreeList.returnPackage(P);
	}
	virtual void put(ProcessWorkPackage * P)
	{
		procPackageFreeList.returnPackage(P);
	}

	virtual void readDone()
	{
		readDoneFlag.set(true);
	}

	#if 0
	virtual void returnPatchVectorUnprocessed(patch_vec_ptr V)
	{
		vecFreeList.put(V);
	}
	#endif

	virtual void returnPatchVectorProcessed(patch_vec_ptr V)
	{
		vecFreeList.put(V);

		if ( ! readDoneFlag.get() )
		{
			ReadWorkPackage * P = dynamic_cast<ReadWorkPackage *>(readPackageFreeList.getPackage());
			*P = ReadWorkPackage(0 /* priority */,readdispid);
			STP.enque(P);
		}

		unfinishedlock.lock();
		--unfinished;
		unfinishedlock.unlock();

		finishedlock.lock();
		finished++;
		finishedlock.unlock();
	}

	virtual void enquePatchVector(std::string const & ref, patch_vec_ptr V)
	{
		ProcessWorkPackage * PWP = dynamic_cast<ProcessWorkPackage *>(procPackageFreeList.getPackage());
		*PWP = ProcessWorkPackage(0 /* priority */,procdispid);
		PWP->ref = ref;
		PWP->V = V;

		unfinishedlock.lock();
		unfinished++;
		unfinishedlock.unlock();

		STP.enque(PWP);
	}
};


static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static int64_t getDefaultMaxClip()
{
	return 50;
}

int mapconstoraw(libmaus2::util::ArgParser const & arg)
{
	std::string const fullfn = arg[0];
	std::string const fragmentfn = arg[1];

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	std::string const sortfile = tmpfilebase + ".sort";
	libmaus2::util::TempFileRemovalContainer::addTempFile(sortfile);

	int const w_err = arg.uniqueArgPresent("w") ? arg.getUnsignedNumericArg<uint64_t>("w") : libmaus2::lcs::NNP::getDefaultMaxWindowError();
	int const w_back = arg.uniqueArgPresent("b") ? arg.getUnsignedNumericArg<uint64_t>("b") : libmaus2::lcs::NNP::getDefaultMaxBack();
	int64_t const maxclip = arg.uniqueArgPresent("c") ? arg.getUnsignedNumericArg<uint64_t>("c") : getDefaultMaxClip();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<OutputEntry> SBO(sortfile);

	WorkControl WC(numthreads,fullfn,fragmentfn,SBO,w_err,w_back,maxclip);
	WC.start();
	WC.wait();

	std::cerr << "[V] processing finished" << std::endl;

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<OutputEntry>::merger_ptr_type Pmerger(SBO.getMerger());
	OutputEntry OE;
	while ( Pmerger->getNext(OE) )
	{
		std::cout << OE.out;
		std::cerr << OE.err;
	}
	Pmerger.reset();
	return EXIT_SUCCESS;
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
			std::cerr << "usage: " << arg.progname << " [options] uncorrected.fasta corrected_fragments.fasta\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return mapconstoraw(arg);
		}

	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
