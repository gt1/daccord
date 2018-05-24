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
#if ! defined(HANDLECONTEXT_HPP)
#define HANDLECONTEXT_HPP

#include <DebruijnGraphContainer.hpp>
#include <ActiveElement.hpp>
#include <TraceType.hpp>
#include <DecodedReadContainer.hpp>
#include <ComputeOffsetLikely.hpp>
#include <DebruijnGraph.hpp>
#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/SuffixArrayLCS.hpp>
#include <libmaus2/dazzler/align/OverlapDataInterface.hpp>
#include <libmaus2/dazzler/align/LasRangeDecoder.hpp>
#include <libmaus2/sorting/ParallelStableSort.hpp>

// #define HANDLE_DEBUG

struct HandleContext // : public DebruijnGraphContainer
{
	typedef HandleContext this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	#if defined(HANDLE_DEBUG)
	static void printDebugInfo(
		std::ostream &
			logstr
			,
		std::ostream &
			#if defined(HANDLE_DEEP_DEBUG)
			err
			#endif
			,
		uint64_t const
			kmersize
			,
		uint64_t const
			filterfreq
			,
		uint8_t const *
			refu
			,
		std::pair<uint64_t,uint64_t> const
			slref
			,
		DebruijnGraphInterface &
			DG
		#if defined(HANDLE_DEEP_DEBUG)
			,
		libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > const & MA,
		uint64_t const MAo
		#endif
	)
	{
		#if defined(HANDLE_DEBUG)
		logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] refseq " << std::string(refu,refu+slref.first) << std::endl;
		#endif

		#if defined(HANDLE_DEEP_DEBUG)
		DG.toDot(logstr);
		#endif

		#if defined(HANDLE_DEEP_DEBUG)
		DG.print(err);
		#endif

		#if defined(HANDLE_DEBUG)
		uint64_t const numk = (slref.first >= kmersize) ? slref.first - kmersize + 1 : 0;
		for ( uint64_t i = 0; i < numk; ++i )
		{
			uint64_t w = 0;
			for ( uint64_t j = 0; j < kmersize; ++j )
			{
				w <<= 2;
				w |= libmaus2::fastx::mapChar(refu[i+j]);
			}

			if ( DG.getNodeVirtual(w) )
			{
				#if defined(HANDLE_DEEP_DEBUG)
				logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] ref " << i << " " << std::string(refu+i,refu+i+kmersize) << " node " << DG.printNode(*DG.getNodeVirtual(w)) << std::endl;
				#endif
			}
			else
			{
				logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] ref " << i << " " << std::string(refu+i,refu+i+kmersize) << " no node " << std::endl;
			}
		}
		for ( uint64_t i = 0; i+1 < numk; ++i )
		{
			uint64_t w0 = 0;
			for ( uint64_t j = 0; j < kmersize; ++j )
			{
				w0 <<= 2;
				w0 |= libmaus2::fastx::mapChar(refu[i+j]);
			}
			uint64_t w1 = 0;
			for ( uint64_t j = 0; j < kmersize; ++j )
			{
				w1 <<= 2;
				w1 |= libmaus2::fastx::mapChar(refu[i+j+1]);
			}

			if ( DG.getNodeVirtual(w0) && DG.getNodeVirtual(w1) && !DG.isEdgeActiveVirtual(w0,w1) )
			{
				logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] ref " << i << " "
					<< std::string(refu+i,refu+i+kmersize) << " -> " << std::string(refu+i+1,refu+i+1+kmersize) << " no link" << std::endl;
				Links L;
				DG.getActiveSuccessorsVirtual(w0,L);
				for ( uint64_t i = 0; i < L.size(); ++i )
					logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] ref " << i << " active sym " << L.getSym(i) << " freq " << L.getFreq(i) << std::endl;
				DG.getSuccessorsVirtual(w0,L);
				for ( uint64_t i = 0; i < L.size(); ++i )
					logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] ref " << i << " all sym " << L.getSym(i) << " freq " << L.getFreq(i) << std::endl;
			}
			// isEdgeActive
		}

		for ( uint64_t i = 0; i+1 < numk; ++i )
		{
			uint64_t w0 = 0;
			for ( uint64_t j = 0; j < kmersize; ++j )
			{
				w0 <<= 2;
				w0 |= libmaus2::fastx::mapChar(refu[i+j]);
			}
			uint64_t w1 = 0;
			for ( uint64_t j = 0; j < kmersize; ++j )
			{
				w1 <<= 2;
				w1 |= libmaus2::fastx::mapChar(refu[i+j+1]);
			}

			logstr << "[filterfreq=" << kmersize << "," << filterfreq << "," << i << "," << numk-i-1 << "] "
				<< std::string(refu+i,refu+i+kmersize) << " -> " << std::string(refu+i+1,refu+i+1+kmersize);

			if ( DG.getNodeVirtual(w0) && DG.getNodeVirtual(w1) && DG.isEdgeActiveVirtual(w0,w1) )
			{
				logstr << " edge weight " << DG.getActiveEdgeWeightVirtual(w0,w1);

				std::string sx(refu+i+1,refu+i+1+kmersize);
				for ( int64_t j = 0; j < 3; ++j )
				{
					if ( libmaus2::fastx::mapChar(sx.back()) != j )
					{
						std::string sy = sx;
						sy.back() = libmaus2::fastx::remapChar(j);

						uint64_t w2 = 0;
						for ( uint64_t j = 0; j < kmersize; ++j )
						{
							w2 <<= 2;
							w2 |= libmaus2::fastx::mapChar(sy[j]);
						}

						if ( DG.getNodeVirtual(w0) && DG.getNodeVirtual(w2) && DG.isEdgeActiveVirtual(w0,w2) )
						{
							logstr << " | " << std::string(refu+i,refu+i+kmersize) << " -> " << sy << " weight " << DG.getActiveEdgeWeightVirtual(w0,w2);
						}
					}
				}
			}
			else
			{
				logstr << " no edge";
				if ( ! DG.getNodeVirtual(w0) )
					logstr << " no source";
				else
					logstr << " source";
				if ( ! DG.getNodeVirtual(w1) )
					logstr << " no target";
				else
					logstr << " target";

				if ( DG.getNodeVirtual(w0) )
				{
					Links L;
					DG.getSuccessorsVirtual(w0,L);
					for ( uint64_t i = 0; i < L.size(); ++i )
						logstr << " (" << libmaus2::fastx::remapChar(L.getSym(i)) << "," << L.getFreq(i) << ")";
				}
			}

			logstr << std::endl;

			DG.printStretches(w0,logstr);
		}
		#endif

		#if defined(HANDLE_DEEP_DEBUG)
		for ( uint64_t i = 0; i < MAo; ++i )
		{
			logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] S " << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
		}
		#endif
	}
	#endif

	struct PileElement
	{
		int64_t apos;
		int64_t apre;
		char sym;

		PileElement() {}
		PileElement(
			int64_t const rapos,
			int64_t const rapre,
			char const rsym
		) : apos(rapos), apre(rapre), sym(rsym) {}

		bool operator<(PileElement const & P) const
		{
			if ( apos != P.apos )
				return apos < P.apos;
			else if ( apre != P.apre )
				return apre < P.apre;
			else
				return sym < P.sym;
		}

		std::string toString() const
		{
			std::ostringstream ostr;
			ostr << "PileElement(apos=" << apos << ",apre=" << apre << ",sym=" << sym << ")";
			return ostr.str();
		}
	};

	uint64_t const maxalign;
	uint64_t const windowsize;
	uint64_t const advancesize;
	int64_t const tspace;
	OffsetLikely const & offsetLikely;

	// produce full sequence (even if uncorrected)
	bool const producefull;
	int const verbose;
	uint64_t const numthreads;
	uint64_t const minwindowcov;
	uint64_t const eminrate;
	uint64_t const minlen;
	int64_t const minfilterfreq;
	int64_t const maxfilterfreq;
	std::string const eq20;

	libmaus2::parallel::SynchronousCounter<uint64_t> & wellcounter;

	libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & readDataFreeList;
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList;
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList2;

	// trace free list
	libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> & traceFreeList;

	struct ThreadContext : public DebruijnGraphContainer
	{
		typedef ThreadContext this_type;
		typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
		typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

		// aligner
		libmaus2::lcs::Aligner::unique_ptr_type PNP;
		libmaus2::lcs::Aligner & NP;

		libmaus2::autoarray::AutoArray < PileElement > PV;
		libmaus2::autoarray::AutoArray < PileElement > NPV;
		libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > MA;

		// end of interval heap
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > E;

		libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type> Mtraces;

		libmaus2::autoarray::AutoArray<double> VVVV;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > PVI;

		libmaus2::autoarray::AutoArray<uint64_t> Vdist;
		libmaus2::autoarray::AutoArray<std::pair<uint64_t,uint64_t> > VPdist;

		libmaus2::autoarray::AutoArray<char> CO;

		// set of active traces/reads
		std::map < uint64_t, ActiveElement > activeset;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > Vend;

		libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > Atrace;

		ThreadContext(
			double const est_cor,
			uint64_t const kmersizelow,
			uint64_t const kmersizehigh,
			std::map < uint64_t, KmerLimit::shared_ptr_type > const & MKL
		)
		: DebruijnGraphContainer(est_cor,kmersizelow,kmersizehigh,MKL), PNP(DebruijnGraphBase::getAligner()), NP(*PNP), E(1024)
		{

		}
	};

	libmaus2::autoarray::AutoArray < ThreadContext::unique_ptr_type > Pthreadcontext;

	HandleContext(
		uint64_t const rmaxalign,
		libmaus2::parallel::SynchronousCounter<uint64_t> & rwellcounter,
		uint64_t const rwindowsize,
		uint64_t const radvancesize,
		// reads
		libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & rreadDataFreeList,
		libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & rreadDecoderFreeList,
		libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & rreadDecoderFreeList2,
		// trace free list
		libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> & rtraceFreeList,
		// trace point space for overlaps
		int64_t const rtspace,
		// kmer length
		OffsetLikely const & roffsetLikely,
		// produce full sequence (even if uncorrected)
		bool const rproducefull,
		double const est_cor,
		uint64_t const kmersizelow,
		uint64_t const kmersizehigh,
		int const rverbose,
		uint64_t const rminwindowcov,
		uint64_t const reminrate,
		uint64_t const rminlen,
		int64_t const rminfilterfreq,
		int64_t const rmaxfilterfreq,
		std::map < uint64_t, KmerLimit::shared_ptr_type > const & MKL,
		uint64_t const rnumthreads
	) : // DebruijnGraphContainer(est_cor,kmersizelow,kmersizehigh,MKL),
	    maxalign(rmaxalign),
	    windowsize(rwindowsize), advancesize(radvancesize),
	    tspace(rtspace), offsetLikely(roffsetLikely),
	    producefull(rproducefull),
	    verbose(rverbose),
	    numthreads(rnumthreads),
	    minwindowcov(rminwindowcov),
	    eminrate(reminrate),
	    minlen(rminlen),
	    minfilterfreq(rminfilterfreq),
	    maxfilterfreq(rmaxfilterfreq),
	    eq20(20,'='),
	    wellcounter(rwellcounter),
	    readDataFreeList(rreadDataFreeList),
	    readDecoderFreeList(rreadDecoderFreeList),
	    readDecoderFreeList2(rreadDecoderFreeList2),
	    traceFreeList(rtraceFreeList),
	    Pthreadcontext(numthreads)
	    #if 0
	    PNP(DebruijnGraphBase::getAligner()), NP(*PNP),
	    E(1024)
	    #endif
	{
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			ThreadContext::unique_ptr_type tptr(new ThreadContext(est_cor,kmersizelow,kmersizehigh,MKL));
			Pthreadcontext[i] = UNIQUE_PTR_MOVE(tptr);
		}
	}

	struct Windows
	{
		uint64_t const l;
		uint64_t const a;
		uint64_t const w;
		uint64_t const n;

		// compute number of windows
		static uint64_t computeN(uint64_t const l, uint64_t const a, uint64_t const w)
		{
			uint64_t const npre = (l + a >= w) ? ((l + a - w) / a) : 0;

			if ( npre )
			{
				if ( (npre-1)*a+w == l )
					return npre;
				else
					return npre+1;
			}
			else
			{
				if ( l >= w )
					return 1;
				else
					return 0;
			}
		}

		Windows(uint64_t const rl, uint64_t const ra, uint64_t const rw)
		: l(rl), a(ra), w(rw), n( computeN(l,a,w))
		{

		}

		uint64_t size() const
		{
			return n;
		}

		std::pair<uint64_t,uint64_t> operator[](uint64_t const i) const
		{
			assert ( i < size() );

			if ( i*a+w <= l )
				return std::pair<uint64_t,uint64_t>(i*a,i*a+w);
			else
				return std::pair<uint64_t,uint64_t>(l-w,l);
		}

		uint64_t offset(uint64_t const i) const
		{
			if ( i+1 < size() )
			{
				std::pair<uint64_t,uint64_t> const Wi = operator[](i);
				std::pair<uint64_t,uint64_t> const Wi1 = operator[](i+1);

				assert ( Wi1.first > Wi.first );

				return Wi1.first - Wi.first;
			}
			else
			{
				return 0;
			}
		}
	};

	void operator()(
		std::ostream & out,
		std::ostream & err,
		// overlaps
		libmaus2::dazzler::align::Overlap const * const ita,
		libmaus2::dazzler::align::Overlap const * const ite
		#if defined(HANDLE_DEBUG)
			,
		// BAM alignment of OVL-A read to reference
		libmaus2::bambam::BamAlignment const & algn,
		//
		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> const & /* Vbam */,
		// reference text
		std::string const & text,
		std::map<double,uint64_t> & ED
		#endif
	)
	{
		ThreadContext & threadcontext = *(Pthreadcontext[0]);
		libmaus2::lcs::Aligner & NP = threadcontext.NP;
		libmaus2::autoarray::AutoArray < PileElement > & PV = threadcontext.PV;
		libmaus2::autoarray::AutoArray < PileElement > & NPV = threadcontext.NPV;
		libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > & MA = threadcontext.MA;
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & E = threadcontext.E;
		libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type> & Mtraces = threadcontext.Mtraces;
		libmaus2::autoarray::AutoArray<double> & VVVV = threadcontext.VVVV;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > & PVI = threadcontext.PVI;
		libmaus2::autoarray::AutoArray<uint64_t> & Vdist = threadcontext.Vdist;
		libmaus2::autoarray::AutoArray<std::pair<uint64_t,uint64_t> > & VPdist = threadcontext.VPdist;
		libmaus2::autoarray::AutoArray<char> & CO = threadcontext.CO;
		std::map < uint64_t, ActiveElement > & activeset = threadcontext.activeset;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > & Vend = threadcontext.Vend;
		// libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > & Atrace = threadcontext.Atrace;

		assert ( E.empty() );
		assert ( activeset.empty() );

		DecodedReadContainer RC(readDataFreeList,readDecoderFreeList);
		DecodedReadContainer RC2(readDataFreeList,readDecoderFreeList2);

		// number of overlaps with A read
		uint64_t const nintv = ite-ita;

		// resize Vend if necessary
		if ( Vend.size() < nintv )
		{
			Vend = libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > ();
			Vend = libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > (nintv,false);
		}
		// collect (bread,aepos) pairs
		for ( uint64_t z = 0; z < nintv; ++z )
			Vend[z] = std::pair<uint64_t,uint64_t>(ita[z].bread,ita[z].path.aepos);
		// sort
		std::sort(Vend.begin(),Vend.begin()+nintv);
		// keep last end for each bread
		uint64_t ovend = 0;
		{
			uint64_t low = 0;
			while ( low < nintv )
			{
				uint64_t high = low+1;
				while ( high < nintv && Vend[low].first == Vend[high].first )
					++high;

				Vend[ovend] = Vend[high-1];
				std::swap(Vend[ovend].first,Vend[ovend].second);
				ovend++;

				low = high;
			}
		}
		std::sort(Vend.begin(),Vend.begin()+ovend);
		uint64_t ivend = 0;

		uint64_t maxaepos = 0;
		for ( uint64_t z = 0; z < nintv; ++z )
			if ( ita[z].path.aepos > static_cast<int64_t>(maxaepos) )
				maxaepos = ita[z].path.aepos;

		#if ! defined(WINDOWALIGNAVOID)
		bool const windowstracepointaligned = ((windowsize % tspace) == 0) && ((advancesize % tspace) == 0) && ((maxaepos % tspace) == 0);
		#else
		bool const windowstracepointaligned = false;
		#endif

		// std::cerr << "nintv=" << nintv << std::endl;

		double maxerate = 0.0;
		double minerate = 1.0;
		for ( uint64_t i = 0; i < nintv; ++i )
		{
			double const erate = ita[i].getErrorRate();
			if ( erate > maxerate )
				maxerate = erate;
			if ( erate < minerate )
				minerate = erate;
		}
		double const ediv = (maxerate > minerate) ? (maxerate - minerate) : 1.0;

		#if 0
		for ( uint64_t i = 0; i < nintv; ++i )
			err << (ita[i].getErrorRate() - minerate) / ediv << std::endl;
		#endif

		int64_t const aread = nintv ? ita[0].aread : -1;
		for ( uint64_t i = 1; i < nintv; ++i )
			assert ( ita[i].aread == aread );

		libmaus2::timing::RealTimeClock handlertc; handlertc.start();

		typedef std::pair<uint64_t,uint64_t> upair;

		int64_t const aid = (ita != ite) ? ita->aread : -1;

		#if defined(HANDLE_DEBUG)
		// cig op vector
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
		// number of cigar operations in ground truth alignment
		uint32_t const ncigar = algn.getCigarOperations(cigop);
		// get trace from free list
		trace_type::shared_ptr_type Pbamtrace = traceFreeList.get();
		// turn cigar of reference/ground truth to trace
		libmaus2::bambam::CigarStringParser::cigarToTrace(cigop.begin(),cigop.begin()+ncigar,*Pbamtrace,true /* ignore unknown */);
		// part of text covered by read
		bool const bamrev = algn.isReverse();
		std::string const psreftext = text.substr(algn.getPos()-algn.getFrontDel(),algn.getReferenceLength());
		// reference text matched by A read
		std::string const sreftext = bamrev ? libmaus2::fastx::reverseComplementUnmapped(psreftext) : psreftext;
		// get start and end pointer
		libmaus2::lcs::AlignmentTraceContainer::step_type * rta = Pbamtrace->ta;
		libmaus2::lcs::AlignmentTraceContainer::step_type * rte = Pbamtrace->te;
		// reverse vector if we are on the reverse complement
		if ( bamrev )
			std::reverse(rta,rte);
		// pointer to reference string
		uint8_t const * refu = reinterpret_cast<uint8_t const *>(sreftext.c_str());
		// front soft clipping
		uint64_t frontsoftclip = algn.getFrontSoftClipping();
		#endif

		// traces
		if ( Mtraces.size() < nintv )
		{
			Mtraces = libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type>();
			Mtraces = libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type>(nintv,false);

			for ( uint64_t z = 0; z < nintv; ++z )
				Mtraces[z] = trace_type::shared_ptr_type();
		}


		#if 0
		uint64_t const ylimit = (maxaepos + advancesize >= windowsize) ? ((maxaepos + advancesize - windowsize) / advancesize) : 0;
		assert ( ylimit * advancesize + windowsize > maxaepos );
		#endif

		uint64_t PVo = 0;

		int failcount = 0;
		int insufcount = 0;

		#if defined(HANDLE_DEBUG)
		libmaus2::lcs::NP HDNP;
		#endif

		Windows W(maxaepos,advancesize,windowsize);

		assert (
			(
				maxaepos < windowsize
				&&
				W.size() == 0
			)
			||
			(
				W.size()
				&&
				W[W.size()-1].second == maxaepos
			)
		);

		char const * emptystr = "";
		std::ostringstream logstr;

		uint64_t z = 0;
		for ( uint64_t y = 0; y < W.size(); ++y )
		{
			bool printlogstr = false;

			logstr.str(emptystr);
			logstr.clear();

			// start of window on A
			uint64_t const astart = W[y].first; // y * advancesize;
			uint64_t const aend   = W[y].second; // astart + windowsize;
			assert ( aend <= maxaepos );

			// remove reads no longer needed from read cache
			while ( ivend < ovend && Vend[ivend].first < astart )
			{
				uint64_t const id = Vend[ivend++].second;
				RC2.erase(id);
			}

			logstr << eq20 << " aread=" << aid << " window y=" << y << "/" << W.size() << " [" << astart << ',' << aend << ')' << '\n';

			assert (
				(!windowstracepointaligned)
				||
				((astart % tspace) == 0)
			);
			assert (
				(!windowstracepointaligned)
				||
				((aend % tspace) == 0)
			);


			#if defined(HANDLE_DEBUG)
			// bases used on reference
			uint64_t const refcov = frontsoftclip >= windowsize ? 0 : windowsize - frontsoftclip;

			// part of reference not handled by front soft clipping
			std::pair<uint64_t,uint64_t> const advref = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refcov);
			assert ( advref.first == refcov || rta + advref.second == rte );
			std::pair<uint64_t,uint64_t> const slref = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advref.second);
			#endif

			// add new active intervals
			while ( z < nintv && static_cast<int64_t>(astart) >= ita[z].path.abpos )
			{
				if ( ita[z].path.aepos >= static_cast<int64_t>(astart) )
				{
					if ( ! windowstracepointaligned )
					{
						uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread));
						uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread));

						trace_type::shared_ptr_type Ptrace = traceFreeList.get();

						ita[z].computeTrace(ua,ub,tspace,*Ptrace,NP);

						Mtraces[z] = Ptrace;
					}

					// offset from start of window
					uint64_t const aoff = astart - ita[z].path.abpos;

					uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)) + astart;
					assert ( (!windowstracepointaligned) || ((((ua - reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread))) % tspace) == 0)));
					assert (((ua - reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)))) == static_cast< ::std::ptrdiff_t > (astart));

					libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = 0;
					libmaus2::lcs::AlignmentTraceContainer::step_type const * te = 0;
					uint64_t uboff;

					if ( windowstracepointaligned )
					{
						uboff = ita[z].path.bbpos + ita[z].getBBlockOffset(((astart - ita[z].path.abpos)+tspace-1)/tspace);
					}
					else
					{
						// get trace
						trace_type::shared_ptr_type Ptrace = Mtraces[z];
						ta = Ptrace->ta;
						te = Ptrace->te;

						// see how many operations there are up to start of region we are interested in
						std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(ta,te,aoff);

						bool const ok = adv.first == aoff;
						if ( ! ok )
						{
							std::cerr << "adv.first=" << adv.first << " aoff=" << aoff << " astart=" << astart << std::endl;
							std::cerr << ita[z] << std::endl;
							assert ( adv.first == aoff );
						}

						uboff = ita[z].path.bbpos + libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ta,ta+adv.second).second;

						// advance in trace
						ta += adv.second;
					}

					uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread)) + uboff;

					// add to active set
					uint64_t const escore = static_cast<uint64_t>(((ita[z].getErrorRate() - minerate) / ediv) * std::numeric_limits<uint32_t>::max());
					assert ( escore <= std::numeric_limits<uint32_t>::max() );
					uint64_t const eindex = (escore<<32) | z;

					activeset[eindex] = ActiveElement(ua,ub,ta,te,uboff,ita[z].getErrorRate());

					// add to erase list
					E.pushBump(upair(ita[z].path.aepos,eindex));
				}

				// next
				z += 1;
			}
			// cleanup
			while ( (! (E.empty())) && E.top().first < aend )
			{
				upair UP = E.pop();
				uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
				trace_type::shared_ptr_type Ptrace = Mtraces[z];
				traceFreeList.put(Ptrace);
				Mtraces[z] = trace_type::shared_ptr_type();
				activeset.erase(UP.second);
			}

			uint64_t MAo = 0;

			uint8_t const * w_ua = activeset.size() ? activeset.begin()->second.ua : 0;

			// iterate over active read list
			for ( std::map<uint64_t,ActiveElement>::iterator s_ita = activeset.begin(); s_ita != activeset.end(); ++s_ita )
			{
				// active element
				ActiveElement & AE = s_ita->second;

				// read id
				uint64_t const key = s_ita->first & std::numeric_limits<uint32_t>::max();
				// get overlap
				libmaus2::dazzler::align::Overlap const & OVL = ita[key];

				#if ! defined(NDEBUG)
				// sanity checks
				assert ( OVL.path.abpos <= static_cast<int64_t>(astart ) );
				assert ( OVL.path.aepos >= static_cast<int64_t>(aend   ) );
				#endif

				uint64_t bwindowsize;
				uint64_t badvancesize;

				if ( windowstracepointaligned )
				{
					uint64_t const aoff = astart - ita[key].path.abpos;

					bwindowsize = 0;

					for ( uint64_t j = 0; j < windowsize / tspace; ++j )
						bwindowsize += OVL.path.path[ (aoff + tspace - 1) / tspace + j ] . second;

					badvancesize = 0;

					for ( uint64_t j = 0; j < advancesize / tspace; ++j )
						badvancesize += OVL.path.path[ (aoff + tspace - 1) / tspace + j ] . second;
				}
				else
				{
					// get end of region in trace
					std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,windowsize);

					bool const ok = adv.first == windowsize;

					if ( ! ok )
					{
						std::cerr << OVL << std::endl;
						std::cerr << "astart=" << astart << std::endl;
						std::cerr << "aend=" << aend << std::endl;
					}

					assert ( adv.first == windowsize );
					// string length
					std::pair<uint64_t,uint64_t> const sl = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+adv.second);

					bwindowsize = sl.second;

					// compute how much we advance to set up for next window
					std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,W.offset(y)); // advancesize
					// in string length used on A and B read
					std::pair<uint64_t,uint64_t> const sladv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+advadv.second);
					// update trace pointer
					AE.ta += advadv.second;

					badvancesize = sladv.second;
				}

				// push A read if this is the first instance of this loop
				if (
					! MAo
					&&
					(&RC != &RC2)
				)
					MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ua,windowsize));
				if ( MAo < maxalign )
				{
					// push B read
					MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ub,bwindowsize));
				}

				// update active element by advancing to next window
				AE.ua += W.offset(y); // advancesize
				AE.ub += badvancesize;
				AE.uboff += badvancesize;
			}

			int64_t maxvprodindex = -1;
			if ( MAo )
			{
				// maximum position on first read as minimum of support
				int64_t minSupLen = static_cast<int64_t>(MA[0].second)-1;
				// support maximum
				int64_t maxSupLen = minSupLen;
				for ( uint64_t j = 1; j < MAo; ++j )
				{
					int64_t const lastpos = static_cast<int64_t>(MA[j].second)-1;
					minSupLen = std::min(minSupLen,lastpos);
					maxSupLen = std::max(maxSupLen,lastpos);
				}
				if ( minSupLen < 0 )
					minSupLen = 0;
				if ( maxSupLen < 0 )
					maxSupLen = 0;

				uint64_t supStart = offsetLikely.getSupportLow(minSupLen);
				uint64_t supEnd   = offsetLikely.getSupportHigh(maxSupLen);

				assert ( (supStart == 0) || (minSupLen >= static_cast<int64_t>(offsetLikely.DPnorm[supStart-1].size())) );
				assert ( (supEnd == offsetLikely.DPnorm.size()) ||  (maxSupLen < static_cast<int64_t>(offsetLikely.DPnorm[supEnd].firstsign)) );

				double maxval = std::numeric_limits<double>::min();
				// iterate over ref positions in support
				for ( uint64_t i = supStart; i < supEnd; ++i )
				{
					DotProduct const & DP = offsetLikely.DPnorm[i];

					double vprod = 1.0;

					for ( uint64_t j = 0; j < MAo; ++j )
					{
						uint64_t const len = MA[j].second;
						if ( len )
						{
							uint64_t const lastpos = len-1;

							vprod *= DP[lastpos];
						}
					}

					if ( vprod > maxval )
					{
						maxval = vprod;
						maxvprodindex = i;
					}
				}
			}

			// length estimation via conditional probability failed, try comparing density distributions
			if ( maxvprodindex == -1 )
			{
				int64_t maxoff = -1;
				double maxoffv = std::numeric_limits<double>::min();


				uint64_t Vdisto = 0;
				uint64_t VPdisto = 0;
				for ( uint64_t i = 0; i < MAo; ++i )
					Vdist.push(Vdisto,MA[i].second);
				std::sort(Vdist.begin(),Vdist.begin()+Vdisto);
				{
					uint64_t low = 0;
					while ( low < Vdisto )
					{
						uint64_t high = low+1;
						while ( high < Vdisto && Vdist[high] == Vdist[low] )
							++high;

						VPdist.push(VPdisto,std::pair<uint64_t,uint64_t>(Vdist[low],high-low));

						low = high;
					}
				}

				uint64_t VVVVo = 0;
				for ( uint64_t i = 0; i < VPdisto; ++i )
				{
					while ( ! ( VPdist[i].first < VVVVo ) )
						VVVV.push(VVVVo,0);
					VVVV[VPdist[i].first] = VPdist[i].second-1;
				}

				for ( uint64_t i = 0; i < offsetLikely.DPnormSquare.size(); ++i )
				{
					double const v = offsetLikely.DPnormSquare[i].dotproduct(VVVV.begin(),VVVVo);
					if ( v > maxoffv )
					{
						maxoff = i;
						maxoffv = v;
					}
					#if 0
					if ( v > 1e-3 )
						logstr << "prod " << (i+1) << " = " << v << " slref=" << slref.first << std::endl;
					#endif
				}

				if ( maxoff != -1 && maxoffv >= 1e-3 )
				{
					maxvprodindex = maxoff;
					// logstr << "prod " << (maxvprodindex+1) << " = " << maxoffv << " slref=" << slref.first << std::endl;
				}
			}

			#if 0
			if ( y+1 == W.size() )
			{
				std::cerr << "window " << astart << "," << aend << " maxaepos=" << maxaepos << " MAo=" << MAo << std::endl;
			}
			#endif

			if (
				(MAo >= minwindowcov)
				#if 0
				&&
				y == 1624
				#endif
			)
			{
				int64_t const elength = maxvprodindex+1;

				logstr << "elength " << elength
				#if defined(HANDLE_DEBUG)
					<< " slref=" << slref.first
				#endif
					<< std::endl;

				#if 0
				for ( uint64_t i = 0; i < MAo; ++i )
				{
					err << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
				}
				#endif

				bool pathfailed = true;
				int64_t filterfreq = -1;

				uint64_t minindex = 0;
				uint64_t minrate = eminrate;
				DebruijnGraphInterface * minDG = 0;

				for ( uint64_t adgi = 0; adgi < threadcontext.ADG.size(); ++adgi )
				{
					DebruijnGraphInterface & DG = *(threadcontext.ADG[adgi]);

					// logstr << "refseq " << std::string(refu,refu+slref.first) << std::endl;
					filterfreq = maxfilterfreq;

					for ( ; filterfreq >= minfilterfreq ; --filterfreq )
					{
						#if defined(HANDLE_TIME)
						libmaus2::timing::RealTimeClock rtc;
						#endif

						// set up debruijn graph
						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.setup(MA.begin(), MAo);
						#if defined(HANDLE_TIME)
						err << "setup " << rtc.getElapsedSeconds() << std::endl;
						#endif

						// filter by frequency (need at least 2)
						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.filterFreq(std::max(filterfreq,static_cast<int64_t>(1)),MAo);
						#if defined(HANDLE_TIME)
						err << "filterfreq " << rtc.getElapsedSeconds() << std::endl;
						#endif

						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
						#if defined(HANDLE_TIME)
						err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
						#endif

						if ( filterfreq == 0 )
						{
							logstr << "*** trying to fill in gaps" << std::endl;

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.getLevelSuccessors(2);
							#if defined(HANDLE_TIME)
							err << "getLevelSuccessors " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.setupNodes();
							#if defined(HANDLE_TIME)
							err << "setupNodes " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.setupAddHeap(MAo);
							#if defined(HANDLE_TIME)
							err << "setupAddHeap " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
							#if defined(HANDLE_TIME)
							err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
							#endif
						}

						uint64_t mintry = 0;
						uint64_t const maxtries = 3;
						bool lconsok = false;

						do
						{
							#if defined(HANDLE_TIME)
							rtc.start();
							#endif

							#if defined(HANDLE_DEEP_DEBUG)
							DG.toDot(err);
							#endif

							bool const consok = DG.traverse(elength - 4,elength + 4,MA.begin(),MAo,16 /* maxfronpath */,16 /* maxfullpath */);
							#if defined(HANDLE_TIME)
							err << "traverse " << rtc.getElapsedSeconds() << std::endl;
							#endif

							if ( consok )
							{
								std::pair<uint64_t,uint64_t> const MR = DG.checkCandidatesU(MA.begin(),MAo);

								if ( MR.second < minrate )
								{
									lconsok = true;

									#if defined(HANDLE_DEBUG)
									if ( minDG )
									{
										logstr << "replacing result for " << minDG->getKmerSize() << "," << minrate << " by " << DG.getKmerSize() << "," << MR.second << std::endl;
										printlogstr = true;
									}
									#endif

									minrate = MR.second;
									minindex = MR.first;
									minDG = &DG;
								}
								else if ( minDG )
								{
									lconsok = true;
								}
								break;
							}
							else
							{
								if ( ++mintry >= maxtries )
									break;
							}
						} while (
							DG.addNextFromHeap(&logstr)
						);

						bool printdebuginfo = false;

						if ( !lconsok )
						{
							logstr << "[filterfreq=" << DG.getKmerSize() << "," << filterfreq << "] failed " << MAo << std::endl;
							printdebuginfo = true;

							if ( verbose >= 2 )
								printlogstr = true;
						}
						else
						{
							pathfailed = false;
							break;
						}

						if ( printdebuginfo )
						{
							#if defined(HANDLE_DEBUG)
							printDebugInfo(logstr,err,DG.getKmerSize(),filterfreq,refu,slref,DG
								#if defined(HANDLE_DEEP_DEBUG)
								,MA,MAo
								#endif
							);
							#endif
						}
					}
				}

				if ( ! pathfailed )
				{
					assert ( minDG );

					std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(minindex);

					#if defined(HANDLE_DEBUG)
					logstr << eq20 << " aread=" << aid << " window y=" << y << "/" << W.size() << " MAo=" << MAo << " maxvprodindex=" << maxvprodindex
						<< " " << std::string(consensus.first,consensus.second)
						<< std::endl;
					#endif

					#if defined(HANDLE_DEBUG)
					libmaus2::lcs::AlignmentStatistics AS;
					libmaus2::lcs::AlignmentStatistics ASD;
					for ( uint64_t i = 0; i < MAo; ++i )
					{
						#if 0
						std::cerr << "slref.first=" << slref.first << " MA[i].second=" << MA[i].second << std::endl;
						std::cerr << std::string(refu,refu+slref.first) << std::endl;
						std::cerr << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
						#endif
						NP.align(refu,slref.first,MA[i].first,MA[i].second);
						AS += NP.getTraceContainer().getAlignmentStatistics();

						NP.align(consensus.first,consensus.second-consensus.first,MA[i].first,MA[i].second);
						ASD += NP.getTraceContainer().getAlignmentStatistics();
					}

					if ( ASD.getEditDistance() != minrate )
						err << "[!] minrate=" << minrate << " ASD.getEditDistance()=" << ASD.getEditDistance() << " AS.getEditDistance()=" << AS.getEditDistance() << " AS.getErrorRate()=" << AS.getErrorRate() << std::endl;

					if ( minrate > AS.getEditDistance() )
					{
						err <<
							"[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] minrate=" << minrate << " ASD.getEditDistance()=" << ASD.getEditDistance() << " AS.getEditDistance()=" << AS.getEditDistance() << " AS.getErrorRate()=" << AS.getErrorRate() << std::endl;

						logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] consensus error " << minrate << " ref error " << AS.getEditDistance() << std::endl;

						for ( uint64_t i = 0; i < minDG->getNumCandidates(); ++i )
						{
							std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(i);
							logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] candidate "
								<< std::string(consensus.first,consensus.second) << " "
								<< minDG->getCandidateErrorU(MA.begin(),MAo,i) << " weight " << minDG->getCandidateWeight(i) << std::endl;
						}

						#if 1
						libmaus2::lcs::AlignmentTraceContainer ATC;
						HDNP.splitAlign(refu,refu+slref.first,consensus.first,consensus.second,ATC);

						// minDG->NP.align(refu,slref.first,consensus.first,consensus.second-consensus.first);

						libmaus2::lcs::AlignmentPrint::printAlignmentLines(
							logstr,
							refu,slref.first,consensus.first,consensus.second-consensus.first,
							160,
							ATC.ta,ATC.te
						);
						#endif

						printDebugInfo(logstr,err,minDG->getKmerSize(),filterfreq,refu,slref,*minDG
							#if defined(HANDLE_DEEP_DEBUG)
							,MA,MAo
							#endif
						);

						for ( uint64_t i = 0; i < MAo; ++i )
						{
							logstr << "MA[" << i << "]=" << std::string(MA[i].first,MA[i].first + MA[i].second) << std::endl;
						}

						printlogstr = true;
					}
					#endif // HANDLE_DEBUG

					uint8_t const * aread = w_ua; // MA[0].first;
					uint64_t const areadlen = windowsize; // MA[0].second;
					uint8_t const * cdata = consensus.first;
					uint64_t const clen = consensus.second - consensus.first;

					NP.align(aread,areadlen,cdata,clen);

					#if 0 // defined(HANDLE_DEBUG)
					if ( refcov == windowsize )
					{
						double const conserror = NP.getTraceContainer().getAlignmentStatistics().getErrorRate();
						// std::cerr << "referror=" << referror << " conserror=" << conserror << std::endl;
						double const errdif = conserror - referror;
						ED [ errdif ] ++;
					}
					#endif

					uint64_t apos = astart;

					libmaus2::lcs::BaseConstants::step_type const * ta = NP.getTraceContainer().ta;
					libmaus2::lcs::BaseConstants::step_type const * te = NP.getTraceContainer().te;

					while ( ta != te )
					{
						uint64_t numins = 0;
						while ( ta != te && *ta == libmaus2::lcs::BaseConstants::STEP_INS )
						{
							++numins;
							++ta;
						}
						for ( uint64_t i = 0; i < numins; ++i )
						{
							assert ( cdata != consensus.second );
							PileElement PE(apos,(-static_cast<int64_t>(numins)) + static_cast<int64_t>(i),*(cdata++));
							PV.push(PVo,PE);
						}

						if ( ta != te )
						{
							assert ( *ta != libmaus2::lcs::BaseConstants::STEP_INS );

							switch ( *(ta++) )
							{
								case libmaus2::lcs::BaseConstants::STEP_MATCH:
								case libmaus2::lcs::BaseConstants::STEP_MISMATCH:
								{
									assert ( cdata != consensus.second );

									PileElement PE(apos++,0,*(cdata++));
									PV.push(PVo,PE);
									break;
								}
								case libmaus2::lcs::BaseConstants::STEP_DEL:
								{
									PileElement PE(apos++,0,'D');
									PV.push(PVo,PE);
									break;
								}
								default:
									break;
							}
						}
					}

					assert ( apos == aend );
				}

				if ( pathfailed )
					++failcount;
			}
			else
			{
				logstr << "insufficient depth " << MAo << std::endl;
				insufcount += 1;
			}

			#if defined(HANDLE_DEBUG)
			// how much to advance on reference or to reduce front soft clipping
			uint64_t const refadv = frontsoftclip >= W.offset(y) ? 0 : W.offset(y) - frontsoftclip;
			// advancement in trace operations
			std::pair<uint64_t,uint64_t> const advrefadv = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refadv);
			// advancement in string length
			std::pair<uint64_t,uint64_t> const slrefadv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advrefadv.second);
			// update front soft clipping
			frontsoftclip -= std::min(frontsoftclip,W.offset(y));
			// update trace start
			rta += advrefadv.second;
			// update reference pointer
			refu += slrefadv.first;
			#endif

			if ( printlogstr )
				err << logstr.str();
		}

		while ( !E.empty() )
		{
			upair UP = E.pop();
			activeset.erase(UP.second);

			uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
			trace_type::shared_ptr_type Ptrace = Mtraces[z];
			traceFreeList.put(Ptrace);
			Mtraces[z] = trace_type::shared_ptr_type();
		}

		while ( ivend < ovend )
		{
			uint64_t const id = Vend[ivend++].second;
			RC2.erase(id);
		}

		std::sort(PV.begin(),PV.begin()+PVo);

		if ( producefull && (ita != ite) )
		{
			uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita->aread));

			uint64_t next = 0;
			uint64_t low = 0;
			uint64_t NPVo = 0;

			while ( low < PVo )
			{
				// look for consecutive pile elements
				uint64_t high = low+1;
				while ( high < PVo && PV[low].apos == PV[high].apos )
					++high;

				// fill up before elements
				for ( ; static_cast<int64_t>(next) < PV[low].apos; ++next )
					NPV.push(NPVo,PileElement(next,0,::tolower(ua[next])));

				// copy pile elements
				for ( uint64_t i = low; i < high; ++i )
					NPV.push(NPVo,PV[i]);

				next = PV[low].apos+1;
				low = high;
			}

			// fill up until end
			uint64_t const rl = RC.getReadLength(ita->aread);
			for ( ; next < rl; ++next )
				NPV.push(NPVo,PileElement(next,0,::tolower(ua[next])));

			PV.swap(NPV);
			PVo = NPVo;

			for ( uint64_t i = 1; i < PVo; ++i )
				assert ( PV[i-1].apos <= PV[i].apos );
		}

		uint64_t PVIo = 0;
		// std::vector< std::pair<uint64_t,uint64_t> > PVI;

		#if 0
		for ( uint64_t i = 0; i < PVo; ++i )
			std::cerr << PV[i].apos << "," << PV[i].apre << std::endl;
		#endif

		uint64_t il = 0;
		while ( il < PVo )
		{
			uint64_t ih = il+1;

			while (
				ih != PVo
				&&
				(PV[ih].apos - PV[ih-1].apos) <= 1
			)
				++ih;

			uint64_t const first = PV[il].apos;
			uint64_t const last = PV[ih-1].apos;

			if ( last-first >= 100 )
			{
				// std::cerr << "first=" << first << " last=" << last << std::endl;
				PVI.push(PVIo,std::pair<uint64_t,uint64_t>(il,ih));
			}

			il = ih;
		}

		for ( uint64_t z = 0; z < PVIo; ++z )
		{
			std::pair<uint64_t,uint64_t> const P = PVI[z];

			uint64_t const first = PV[P.first].apos;
			uint64_t const last = PV[P.second-1].apos;

			// std::vector<char> CO;
			uint64_t COo = 0;

			int64_t l = P.second;
			int64_t depth = -1;

			while ( l > static_cast<int64_t>(P.first) )
			{
				int64_t h = --l;

				while (
					l >= 0
					&&
					PV[l].apos == PV[h].apos
					&&
					PV[l].apre == PV[h].apre
				)
					--l;

				l += 1;

				for ( int64_t i = l; i <= h; ++i )
				{
					assert ( PV[i].apos == PV[l].apos );
					assert ( PV[i].apre == PV[l].apre );
				}

				uint64_t const ld = (h-l)+1;

				if ( PV[l].apre == 0 )
					depth = ld;

				std::pair<uint64_t,uint64_t> C[] =
				{
					std::pair<uint64_t,uint64_t>(0,'A'),
					std::pair<uint64_t,uint64_t>(0,'C'),
					std::pair<uint64_t,uint64_t>(0,'G'),
					std::pair<uint64_t,uint64_t>(0,'T'),
					std::pair<uint64_t,uint64_t>(0,'D'),
					std::pair<uint64_t,uint64_t>(0,'a'),
					std::pair<uint64_t,uint64_t>(0,'c'),
					std::pair<uint64_t,uint64_t>(0,'g'),
					std::pair<uint64_t,uint64_t>(0,'t'),
					std::pair<uint64_t,uint64_t>(0,0)
				};

				#if defined(PV_DEBUG)
				logstr << "PV " << PV[l].apos << "," << PV[l].apre << " " << depth << " ";
				#endif

				for ( int64_t i = l; i <= h; ++i )
				{
					#if defined(PV_DEBUG)
					logstr << PV[i].sym;
					#endif

					switch ( PV[i].sym )
					{
						case 'A': C[0].first++; break;
						case 'C': C[1].first++; break;
						case 'G': C[2].first++; break;
						case 'T': C[3].first++; break;
						case 'D': C[4].first++; break;
						case 'a': C[5].first++; break;
						case 'c': C[6].first++; break;
						case 'g': C[7].first++; break;
						case 't': C[8].first++; break;
						default: assert(0); break;
					}
				}

				for ( int64_t i = ld; i < depth; ++i )
					C[4].first++;

				std::sort(&C[0],&C[sizeof(C)/sizeof(C[0])],std::greater<std::pair<uint64_t,uint64_t> >());

				#if defined(PV_DEBUG)
				err.put(' ');
				for ( uint64_t i = 0; C[i].first; ++i )
					logstr << "(" << static_cast<char>(C[i].second) << "," << C[i].first << ")";
				err << std::endl;
				#endif

				if ( C[0].first && C[0].second != 'D' )
					CO.push(COo,C[0].second);
			}

			std::reverse(CO.begin(),CO.begin()+COo);

			if ( producefull || COo >= minlen )
			{
				out << '>' << (aid+1) << '/' << wellcounter++ << '/' << first << '_' << first + COo << " A=[" << first << "," << last << "]" << "\n";
				char const * zp = CO.begin();
				char const * ze = zp + COo;
				uint64_t const linewidth = 80;
				while ( zp != ze )
				{
					uint64_t const r = ze-zp;
					uint64_t const tp = std::min(r,linewidth);
					out.write(zp,tp);
					out.put('\n');
					zp += tp;
				}
			}

			#if 0
			{
				std::string const SCO(CO.begin(),CO.begin()+COo);

				// consensus text
				uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(SCO.c_str());
				uint8_t const * uSCOe = uSCO + SCO.size();
				// last read position for consensus
				uint64_t const last = PV[P.second-1].apos;

				std::string const conspart(uSCO, uSCOe);
				std::string const readpart(RC.getForwardRead(ita[0].aread)+first,RC.getForwardRead(ita[0].aread)+last+1);

				libmaus2::lcs::NP NPcons;

				NPcons.np(readpart.begin(),readpart.end(),conspart.begin(),conspart.end());

				int64_t rpos = first;

				libmaus2::lcs::AlignmentTraceContainer::step_type const * consta = NPcons.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * conste = NPcons.te;

				if ( rpos % advancesize != 0 )
				{
					int64_t const skip = advancesize - (rpos % advancesize);

					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,skip);
					assert ( Pcons.first == skip );
					consta += Pcons.second;

					rpos += skip;
				}

				assert ( consta == conste || (rpos % advancesize == 0) );

				while ( libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize).first == windowsize )
				{
					// std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize);
					//double const conserate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(consta,consta+Pcons.second).getErrorRate();
					consta += libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,advancesize).second;
				}
			}
			#endif

			#if defined(HANDLE_DEBUG)
			{
				uint8_t const * uref = reinterpret_cast<uint8_t const *>(sreftext.c_str());
				uint8_t const * urefe = uref + sreftext.size();

				// err << SCO << std::endl;

				uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(SCO.c_str());
				uint8_t const * uSCOe = uSCO + SCO.size();

				libmaus2::lcs::SuffixArrayLCS::LCSResult const lcsres = libmaus2::lcs::SuffixArrayLCS::lcsmin(std::string(uref,urefe),std::string(uSCO,uSCOe));
				libmaus2::lcs::NNP nnp;
				libmaus2::lcs::NNPTraceContainer nnptrace;
				// align ref and consensus
				libmaus2::lcs::NNPAlignResult const nnpres = nnp.align(uref,urefe,lcsres.maxpos_a,uSCO,uSCOe,lcsres.maxpos_b,nnptrace);
				//nnptrace.printTraceLines(err,uref + nnpres.abpos,uSCO + nnpres.bbpos);

				libmaus2::lcs::AlignmentTraceContainer ATC;
				nnptrace.computeTrace(ATC);

				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					uref+nnpres.abpos,
					nnpres.aepos-nnpres.abpos,
					uSCO+nnpres.bbpos,
					nnpres.bepos-nnpres.bbpos,
					80,
					ATC.ta,
					ATC.te
				);

				libmaus2::lcs::AlignmentStatistics const AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(ATC.ta,ATC.te);

				// compare alignments ref/read to cons/read
				std::string const refpart(uref + nnpres.abpos,uref + nnpres.aepos);
				std::string const conspart(uSCO + nnpres.bbpos, uSCO + nnpres.bepos);
				std::string const readpart(RC.getForwardRead(ita[0].aread)+first,RC.getForwardRead(ita[0].aread)+last+1);

				libmaus2::lcs::NP NPref;
				libmaus2::lcs::NP NPcons;

				NPref.np(readpart.begin(),readpart.end(),refpart.begin(),refpart.end());
				NPcons.np(readpart.begin(),readpart.end(),conspart.begin(),conspart.end());

				#if 0
				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					readpart.begin(),
					readpart.size(),
					refpart.begin(),
					refpart.size(),
					80,
					NPref.ta,
					NPref.te
				);

				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					readpart.begin(),
					readpart.size(),
					conspart.begin(),
					conspart.size(),
					80,
					NPcons.ta,
					NPcons.te
				);
				#endif

				int64_t rpos = first;

				libmaus2::lcs::AlignmentTraceContainer::step_type const * refta = NPref.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * refte = NPref.te;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * consta = NPcons.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * conste = NPcons.te;

				if ( rpos % advancesize != 0 )
				{
					int64_t const skip = advancesize - (rpos % advancesize);

					std::pair<uint64_t,uint64_t> const Pref = libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,skip);
					assert ( static_cast<int64_t>(Pref.first) == skip );
					refta += Pref.second;

					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,skip);
					assert ( static_cast<int64_t>(Pcons.first) == skip );
					consta += Pcons.second;

					rpos += skip;
				}

				while ( libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,windowsize).first == windowsize )
				{
					std::pair<uint64_t,uint64_t> const Pref = libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,windowsize);
					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize);

					double const referate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(refta,refta+Pref.second).getErrorRate();
					double const conserate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(consta,consta+Pcons.second).getErrorRate();

					if ( std::abs(conserate-referate) >= 1e-6 )
						err << rpos << " (" << referate << "," << conserate << "," << referate-conserate << ")\n";

					ED[referate - conserate]++;

					refta += libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,advancesize).second;
					consta += libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,advancesize).second;
					rpos += advancesize;
				}

				err << "read id " << ita[0].aread << " " << AS << " " << nnpres << std::endl;

				// err << SCO << std::endl;
			}
			#endif
		}

		// err << "read id " << ita[0].aread << " " << AS << " failcount=" << failcount << " insufcount=" << insufcount << " unbroken=" << Vsufsub << " " << nnpres << " pvdense=" << pvdense << std::endl;

		// return traces
		for ( uint64_t i = 0; i < nintv; ++i )
			assert ( ! Mtraces[i] );
			// if ( Mtraces[i] )
				//traceFreeList.put(Mtraces[i]);

		#if defined(HANDLE_DEBUG)
		// return BAM trace for ground truth
		traceFreeList.put(Pbamtrace);
		#endif

		// DG.printSize(err);

		err << "[V] read id " << ita[0].aread << " time " << handlertc << std::endl;
	}

	void operator()(
		std::ostream & out,
		std::ostream & err,
		// overlaps
		libmaus2::dazzler::align::OverlapDataInterface const * ita,
		libmaus2::dazzler::align::OverlapDataInterface const * ite
		#if defined(HANDLE_DEBUG)
			,
		// BAM alignment of OVL-A read to reference
		libmaus2::bambam::BamAlignment const & algn,
		//
		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> const & /* Vbam */,
		// reference text
		std::string const & text,
		std::map<double,uint64_t> & ED
		#endif
	)
	{
		ThreadContext & threadcontext = *(Pthreadcontext[0]);
		libmaus2::lcs::Aligner & NP = threadcontext.NP;
		libmaus2::autoarray::AutoArray < PileElement > & PV = threadcontext.PV;
		libmaus2::autoarray::AutoArray < PileElement > & NPV = threadcontext.NPV;
		libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > & MA = threadcontext.MA;
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & E = threadcontext.E;
		libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type> & Mtraces = threadcontext.Mtraces;
		libmaus2::autoarray::AutoArray<double> & VVVV = threadcontext.VVVV;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > & PVI = threadcontext.PVI;
		libmaus2::autoarray::AutoArray<uint64_t> & Vdist = threadcontext.Vdist;
		libmaus2::autoarray::AutoArray<std::pair<uint64_t,uint64_t> > & VPdist = threadcontext.VPdist;
		libmaus2::autoarray::AutoArray<char> & CO = threadcontext.CO;
		std::map < uint64_t, ActiveElement > & activeset = threadcontext.activeset;
		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > & Vend = threadcontext.Vend;
		libmaus2::autoarray::AutoArray<std::pair<uint16_t,uint16_t> > & Atrace = threadcontext.Atrace;

		assert ( E.empty() );
		assert ( activeset.empty() );

		DecodedReadContainer RC(readDataFreeList,readDecoderFreeList);
		DecodedReadContainer RC2(readDataFreeList,readDecoderFreeList2);

		// number of overlaps with A read
		uint64_t const nintv = ite-ita;

		// resize Vend if necessary
		if ( Vend.size() < nintv )
		{
			Vend = libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > ();
			Vend = libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > (nintv,false);
		}
		// collect (bread,aepos) pairs
		for ( uint64_t z = 0; z < nintv; ++z )
			Vend[z] = std::pair<uint64_t,uint64_t>(ita[z].bread(),ita[z].aepos());
		// sort
		std::sort(Vend.begin(),Vend.begin()+nintv);
		// keep last end for each bread
		uint64_t ovend = 0;
		{
			uint64_t low = 0;
			while ( low < nintv )
			{
				uint64_t high = low+1;
				while ( high < nintv && Vend[low].first == Vend[high].first )
					++high;

				Vend[ovend] = Vend[high-1];
				std::swap(Vend[ovend].first,Vend[ovend].second);
				ovend++;

				low = high;
			}
		}
		std::sort(Vend.begin(),Vend.begin()+ovend);
		uint64_t ivend = 0;

		uint64_t maxaepos = 0;
		for ( uint64_t z = 0; z < nintv; ++z )
			if ( ita[z].aepos() > static_cast<int64_t>(maxaepos) )
				maxaepos = ita[z].aepos();

		// std::cerr << "nintv=" << nintv << std::endl;

		double maxerate = 0.0;
		double minerate = 1.0;
		for ( uint64_t i = 0; i < nintv; ++i )
		{
			double const erate = ita[i].getErrorRate();
			if ( erate > maxerate )
				maxerate = erate;
			if ( erate < minerate )
				minerate = erate;
		}
		double const ediv = (maxerate > minerate) ? (maxerate - minerate) : 1.0;

		#if 0
		for ( uint64_t i = 0; i < nintv; ++i )
			err << (ita[i].getErrorRate() - minerate) / ediv << std::endl;
		#endif

		int64_t const aread = nintv ? ita[0].aread() : -1;
		for ( uint64_t i = 1; i < nintv; ++i )
			assert ( ita[i].aread() == aread );

		libmaus2::timing::RealTimeClock handlertc; handlertc.start();

		typedef std::pair<uint64_t,uint64_t> upair;

		int64_t const aid = (ita != ite) ? ita->aread() : -1;

		#if defined(HANDLE_DEBUG)
		// cig op vector
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
		// number of cigar operations in ground truth alignment
		uint32_t const ncigar = algn.getCigarOperations(cigop);
		// get trace from free list
		trace_type::shared_ptr_type Pbamtrace = traceFreeList.get();
		// turn cigar of reference/ground truth to trace
		libmaus2::bambam::CigarStringParser::cigarToTrace(cigop.begin(),cigop.begin()+ncigar,*Pbamtrace,true /* ignore unknown */);
		// part of text covered by read
		bool const bamrev = algn.isReverse();
		std::string const psreftext = text.substr(algn.getPos()-algn.getFrontDel(),algn.getReferenceLength());
		// reference text matched by A read
		std::string const sreftext = bamrev ? libmaus2::fastx::reverseComplementUnmapped(psreftext) : psreftext;
		// get start and end pointer
		libmaus2::lcs::AlignmentTraceContainer::step_type * rta = Pbamtrace->ta;
		libmaus2::lcs::AlignmentTraceContainer::step_type * rte = Pbamtrace->te;
		// reverse vector if we are on the reverse complement
		if ( bamrev )
			std::reverse(rta,rte);
		// pointer to reference string
		uint8_t const * refu = reinterpret_cast<uint8_t const *>(sreftext.c_str());
		// front soft clipping
		uint64_t frontsoftclip = algn.getFrontSoftClipping();
		#endif

		// traces
		if ( Mtraces.size() < nintv )
		{
			Mtraces = libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type>();
			Mtraces = libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type>(nintv,false);

			for ( uint64_t z = 0; z < nintv; ++z )
				Mtraces[z] = trace_type::shared_ptr_type();
		}

		uint64_t PVo = 0;

		int failcount = 0;
		int insufcount = 0;

		#if defined(HANDLE_DEBUG)
		libmaus2::lcs::NP HDNP;
		#endif

		Windows W(maxaepos,advancesize,windowsize);

		assert (
			(
				maxaepos < windowsize
				&&
				W.size() == 0
			)
			||
			(
				W.size()
				&&
				W[W.size()-1].second == maxaepos
			)
		);

		char const * emptystr = "";
		std::ostringstream logstr;

		uint64_t z = 0;
		for ( uint64_t y = 0; y < W.size(); ++y )
		{
			bool printlogstr = false;

			logstr.str(emptystr);
			logstr.clear();

			// start of window on A
			uint64_t const astart = W[y].first; // y * advancesize;
			uint64_t const aend   = W[y].second; // astart + windowsize;
			assert ( aend <= maxaepos );

			// remove reads no longer needed from read cache
			while ( ivend < ovend && Vend[ivend].first < astart )
			{
				uint64_t const id = Vend[ivend++].second;
				RC2.erase(id);
			}

			logstr << eq20 << " aread=" << aid << " window y=" << y << "/" << W.size() << " [" << astart << ',' << aend << ')' << '\n';

			#if defined(HANDLE_DEBUG)
			// bases used on reference
			uint64_t const refcov = frontsoftclip >= windowsize ? 0 : windowsize - frontsoftclip;

			// part of reference not handled by front soft clipping
			std::pair<uint64_t,uint64_t> const advref = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refcov);
			assert ( advref.first == refcov || rta + advref.second == rte );
			std::pair<uint64_t,uint64_t> const slref = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advref.second);
			#endif

			// add new active intervals
			while ( z < nintv && static_cast<int64_t>(astart) >= ita[z].abpos() )
			{
				if ( ita[z].aepos() >= static_cast<int64_t>(astart) )
				{
					{
						uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread()));
						uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread()) : RC2.getForwardRead(ita[z].bread()));

						trace_type::shared_ptr_type Ptrace = traceFreeList.get();

						ita[z].computeTrace(Atrace,tspace,ua,ub,*Ptrace,NP);

						Mtraces[z] = Ptrace;
					}

					// offset from start of window
					uint64_t const aoff = astart - ita[z].abpos();

					uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread())) + astart;
					assert (((ua - reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread())))) == static_cast< ::std::ptrdiff_t > (astart));

					libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = 0;
					libmaus2::lcs::AlignmentTraceContainer::step_type const * te = 0;
					uint64_t uboff;

					{
						// get trace
						trace_type::shared_ptr_type Ptrace = Mtraces[z];
						ta = Ptrace->ta;
						te = Ptrace->te;

						// see how many operations there are up to start of region we are interested in
						std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(ta,te,aoff);

						bool const ok = adv.first == aoff;
						if ( ! ok )
						{
							std::cerr << "adv.first=" << adv.first << " aoff=" << aoff << " astart=" << astart << std::endl;
							std::cerr << ita[z] << std::endl;
							assert ( adv.first == aoff );
						}

						uboff = ita[z].bbpos() + libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ta,ta+adv.second).second;

						// advance in trace
						ta += adv.second;
					}

					uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread()) : RC2.getForwardRead(ita[z].bread())) + uboff;

					// add to active set
					uint64_t const escore = static_cast<uint64_t>(((ita[z].getErrorRate() - minerate) / ediv) * std::numeric_limits<uint32_t>::max());
					assert ( escore <= std::numeric_limits<uint32_t>::max() );
					uint64_t const eindex = (escore<<32) | z;

					activeset[eindex] = ActiveElement(ua,ub,ta,te,uboff,ita[z].getErrorRate());

					// add to erase list
					E.pushBump(upair(ita[z].aepos(),eindex));
				}

				// next
				z += 1;
			}
			// cleanup
			while ( (! (E.empty())) && E.top().first < aend )
			{
				upair UP = E.pop();
				uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
				trace_type::shared_ptr_type Ptrace = Mtraces[z];
				traceFreeList.put(Ptrace);
				Mtraces[z] = trace_type::shared_ptr_type();
				activeset.erase(UP.second);
			}

			uint64_t MAo = 0;

			uint8_t const * w_ua = activeset.size() ? activeset.begin()->second.ua : 0;

			// iterate over active read list
			for ( std::map<uint64_t,ActiveElement>::iterator s_ita = activeset.begin(); s_ita != activeset.end(); ++s_ita )
			{
				// active element
				ActiveElement & AE = s_ita->second;

				// read id
				uint64_t const key = s_ita->first & std::numeric_limits<uint32_t>::max();
				// get overlap
				libmaus2::dazzler::align::OverlapDataInterface const & OVL = ita[key];

				#if ! defined(NDEBUG)
				// sanity checks
				assert ( OVL.abpos() <= static_cast<int64_t>(astart ) );
				assert ( OVL.aepos() >= static_cast<int64_t>(aend   ) );
				#endif

				uint64_t bwindowsize;
				uint64_t badvancesize;

				{
					// get end of region in trace
					std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,windowsize);

					bool const ok = adv.first == windowsize;

					if ( ! ok )
					{
						std::cerr << OVL << std::endl;
						std::cerr << "astart=" << astart << std::endl;
						std::cerr << "aend=" << aend << std::endl;
					}

					assert ( adv.first == windowsize );
					// string length
					std::pair<uint64_t,uint64_t> const sl = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+adv.second);

					bwindowsize = sl.second;

					// compute how much we advance to set up for next window
					std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,W.offset(y)); // advancesize
					// in string length used on A and B read
					std::pair<uint64_t,uint64_t> const sladv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+advadv.second);
					// update trace pointer
					AE.ta += advadv.second;

					badvancesize = sladv.second;
				}

				// push A read if this is the first instance of this loop
				if (
					! MAo
					&&
					(&RC != &RC2)
				)
					MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ua,windowsize));
				if ( MAo < maxalign )
				{
					// push B read
					MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ub,bwindowsize));
				}

				// update active element by advancing to next window
				AE.ua += W.offset(y); // advancesize
				AE.ub += badvancesize;
				AE.uboff += badvancesize;
			}

			int64_t maxvprodindex = -1;
			if ( MAo )
			{
				// maximum position on first read as minimum of support
				int64_t minSupLen = static_cast<int64_t>(MA[0].second)-1;
				// support maximum
				int64_t maxSupLen = minSupLen;
				for ( uint64_t j = 1; j < MAo; ++j )
				{
					int64_t const lastpos = static_cast<int64_t>(MA[j].second)-1;
					minSupLen = std::min(minSupLen,lastpos);
					maxSupLen = std::max(maxSupLen,lastpos);
				}
				if ( minSupLen < 0 )
					minSupLen = 0;
				if ( maxSupLen < 0 )
					maxSupLen = 0;

				uint64_t supStart = offsetLikely.getSupportLow(minSupLen);
				uint64_t supEnd   = offsetLikely.getSupportHigh(maxSupLen);

				assert ( (supStart == 0) || (minSupLen >= static_cast<int64_t>(offsetLikely.DPnorm[supStart-1].size())) );
				assert ( (supEnd == offsetLikely.DPnorm.size()) ||  (maxSupLen < static_cast<int64_t>(offsetLikely.DPnorm[supEnd].firstsign)) );

				double maxval = std::numeric_limits<double>::min();
				// iterate over ref positions in support
				for ( uint64_t i = supStart; i < supEnd; ++i )
				{
					DotProduct const & DP = offsetLikely.DPnorm[i];

					double vprod = 1.0;

					for ( uint64_t j = 0; j < MAo; ++j )
					{
						uint64_t const len = MA[j].second;
						if ( len )
						{
							uint64_t const lastpos = len-1;

							vprod *= DP[lastpos];
						}
					}

					if ( vprod > maxval )
					{
						maxval = vprod;
						maxvprodindex = i;
					}
				}
			}

			// length estimation via conditional probability failed, try comparing density distributions
			if ( maxvprodindex == -1 )
			{
				int64_t maxoff = -1;
				double maxoffv = std::numeric_limits<double>::min();


				uint64_t Vdisto = 0;
				uint64_t VPdisto = 0;
				for ( uint64_t i = 0; i < MAo; ++i )
					Vdist.push(Vdisto,MA[i].second);
				std::sort(Vdist.begin(),Vdist.begin()+Vdisto);
				{
					uint64_t low = 0;
					while ( low < Vdisto )
					{
						uint64_t high = low+1;
						while ( high < Vdisto && Vdist[high] == Vdist[low] )
							++high;

						VPdist.push(VPdisto,std::pair<uint64_t,uint64_t>(Vdist[low],high-low));

						low = high;
					}
				}

				uint64_t VVVVo = 0;
				for ( uint64_t i = 0; i < VPdisto; ++i )
				{
					while ( ! ( VPdist[i].first < VVVVo ) )
						VVVV.push(VVVVo,0);
					VVVV[VPdist[i].first] = VPdist[i].second-1;
				}

				for ( uint64_t i = 0; i < offsetLikely.DPnormSquare.size(); ++i )
				{
					double const v = offsetLikely.DPnormSquare[i].dotproduct(VVVV.begin(),VVVVo);
					if ( v > maxoffv )
					{
						maxoff = i;
						maxoffv = v;
					}
					#if 0
					if ( v > 1e-3 )
						logstr << "prod " << (i+1) << " = " << v << " slref=" << slref.first << std::endl;
					#endif
				}

				if ( maxoff != -1 && maxoffv >= 1e-3 )
				{
					maxvprodindex = maxoff;
					// logstr << "prod " << (maxvprodindex+1) << " = " << maxoffv << " slref=" << slref.first << std::endl;
				}
			}

			#if 0
			if ( y+1 == W.size() )
			{
				std::cerr << "window " << astart << "," << aend << " maxaepos=" << maxaepos << " MAo=" << MAo << std::endl;
			}
			#endif

			if (
				(MAo >= minwindowcov)
				#if 0
				&&
				y == 1624
				#endif
			)
			{
				int64_t const elength = maxvprodindex+1;

				logstr << "elength " << elength
				#if defined(HANDLE_DEBUG)
					<< " slref=" << slref.first
				#endif
					<< std::endl;

				#if 0
				for ( uint64_t i = 0; i < MAo; ++i )
				{
					err << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
				}
				#endif

				bool pathfailed = true;
				int64_t filterfreq = -1;

				uint64_t minindex = 0;
				uint64_t minrate = eminrate;
				DebruijnGraphInterface * minDG = 0;

				for ( uint64_t adgi = 0; adgi < threadcontext.ADG.size(); ++adgi )
				{
					DebruijnGraphInterface & DG = *(threadcontext.ADG[adgi]);

					// logstr << "refseq " << std::string(refu,refu+slref.first) << std::endl;
					filterfreq = maxfilterfreq;

					for ( ; filterfreq >= minfilterfreq ; --filterfreq )
					{
						#if defined(HANDLE_TIME)
						libmaus2::timing::RealTimeClock rtc;
						#endif

						// set up debruijn graph
						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.setup(MA.begin(), MAo);
						#if defined(HANDLE_TIME)
						err << "setup " << rtc.getElapsedSeconds() << std::endl;
						#endif

						// filter by frequency (need at least 2)
						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.filterFreq(std::max(filterfreq,static_cast<int64_t>(1)),MAo);
						#if defined(HANDLE_TIME)
						err << "filterfreq " << rtc.getElapsedSeconds() << std::endl;
						#endif

						#if defined(HANDLE_TIME)
						rtc.start();
						#endif
						DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
						#if defined(HANDLE_TIME)
						err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
						#endif

						if ( filterfreq == 0 )
						{
							logstr << "*** trying to fill in gaps" << std::endl;

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.getLevelSuccessors(2);
							#if defined(HANDLE_TIME)
							err << "getLevelSuccessors " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.setupNodes();
							#if defined(HANDLE_TIME)
							err << "setupNodes " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.setupAddHeap(MAo);
							#if defined(HANDLE_TIME)
							err << "setupAddHeap " << rtc.getElapsedSeconds() << std::endl;
							#endif

							#if defined(HANDLE_TIME)
							rtc.start();
							#endif
							DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
							#if defined(HANDLE_TIME)
							err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
							#endif
						}

						uint64_t mintry = 0;
						uint64_t const maxtries = 3;
						bool lconsok = false;

						do
						{
							#if defined(HANDLE_TIME)
							rtc.start();
							#endif

							#if defined(HANDLE_DEEP_DEBUG)
							DG.toDot(err);
							#endif

							bool const consok = DG.traverse(elength - 4,elength + 4,MA.begin(),MAo,16 /* maxfronpath */,16 /* maxfullpath */);
							#if defined(HANDLE_TIME)
							err << "traverse " << rtc.getElapsedSeconds() << std::endl;
							#endif

							if ( consok )
							{
								std::pair<uint64_t,uint64_t> const MR = DG.checkCandidatesU(MA.begin(),MAo);

								if ( MR.second < minrate )
								{
									lconsok = true;

									#if defined(HANDLE_DEBUG)
									if ( minDG )
									{
										logstr << "replacing result for " << minDG->getKmerSize() << "," << minrate << " by " << DG.getKmerSize() << "," << MR.second << std::endl;
										printlogstr = true;
									}
									#endif

									minrate = MR.second;
									minindex = MR.first;
									minDG = &DG;
								}
								else if ( minDG )
								{
									lconsok = true;
								}
								break;
							}
							else
							{
								if ( ++mintry >= maxtries )
									break;
							}
						} while (
							DG.addNextFromHeap(&logstr)
						);

						bool printdebuginfo = false;

						if ( !lconsok )
						{
							logstr << "[filterfreq=" << DG.getKmerSize() << "," << filterfreq << "] failed " << MAo << std::endl;
							printdebuginfo = true;

							if ( verbose >= 2 )
								printlogstr = true;
						}
						else
						{
							pathfailed = false;
							break;
						}

						if ( printdebuginfo )
						{
							#if defined(HANDLE_DEBUG)
							printDebugInfo(logstr,err,DG.getKmerSize(),filterfreq,refu,slref,DG
								#if defined(HANDLE_DEEP_DEBUG)
								,MA,MAo
								#endif
							);
							#endif
						}
					}
				}

				if ( ! pathfailed )
				{
					assert ( minDG );

					std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(minindex);

					#if defined(HANDLE_DEBUG)
					logstr << eq20 << " aread=" << aid << " window y=" << y << "/" << W.size() << " MAo=" << MAo << " maxvprodindex=" << maxvprodindex
						<< " " << std::string(consensus.first,consensus.second)
						<< std::endl;
					#endif

					#if defined(HANDLE_DEBUG)
					libmaus2::lcs::AlignmentStatistics AS;
					libmaus2::lcs::AlignmentStatistics ASD;
					for ( uint64_t i = 0; i < MAo; ++i )
					{
						#if 0
						std::cerr << "slref.first=" << slref.first << " MA[i].second=" << MA[i].second << std::endl;
						std::cerr << std::string(refu,refu+slref.first) << std::endl;
						std::cerr << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
						#endif
						NP.align(refu,slref.first,MA[i].first,MA[i].second);
						AS += NP.getTraceContainer().getAlignmentStatistics();

						NP.align(consensus.first,consensus.second-consensus.first,MA[i].first,MA[i].second);
						ASD += NP.getTraceContainer().getAlignmentStatistics();
					}

					if ( ASD.getEditDistance() != minrate )
						err << "[!] minrate=" << minrate << " ASD.getEditDistance()=" << ASD.getEditDistance() << " AS.getEditDistance()=" << AS.getEditDistance() << " AS.getErrorRate()=" << AS.getErrorRate() << std::endl;

					if ( minrate > AS.getEditDistance() )
					{
						err <<
							"[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] minrate=" << minrate << " ASD.getEditDistance()=" << ASD.getEditDistance() << " AS.getEditDistance()=" << AS.getEditDistance() << " AS.getErrorRate()=" << AS.getErrorRate() << std::endl;

						logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] consensus error " << minrate << " ref error " << AS.getEditDistance() << std::endl;

						for ( uint64_t i = 0; i < minDG->getNumCandidates(); ++i )
						{
							std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(i);
							logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] candidate "
								<< std::string(consensus.first,consensus.second) << " "
								<< minDG->getCandidateErrorU(MA.begin(),MAo,i) << " weight " << minDG->getCandidateWeight(i) << std::endl;
						}

						#if 1
						libmaus2::lcs::AlignmentTraceContainer ATC;
						HDNP.splitAlign(refu,refu+slref.first,consensus.first,consensus.second,ATC);

						// minDG->NP.align(refu,slref.first,consensus.first,consensus.second-consensus.first);

						libmaus2::lcs::AlignmentPrint::printAlignmentLines(
							logstr,
							refu,slref.first,consensus.first,consensus.second-consensus.first,
							160,
							ATC.ta,ATC.te
						);
						#endif

						printDebugInfo(logstr,err,minDG->getKmerSize(),filterfreq,refu,slref,*minDG
							#if defined(HANDLE_DEEP_DEBUG)
							,MA,MAo
							#endif
						);

						for ( uint64_t i = 0; i < MAo; ++i )
						{
							logstr << "MA[" << i << "]=" << std::string(MA[i].first,MA[i].first + MA[i].second) << std::endl;
						}

						printlogstr = true;
					}
					#endif // HANDLE_DEBUG

					uint8_t const * aread = w_ua; // MA[0].first;
					uint64_t const areadlen = windowsize; // MA[0].second;
					uint8_t const * cdata = consensus.first;
					uint64_t const clen = consensus.second - consensus.first;

					NP.align(aread,areadlen,cdata,clen);

					#if 0 // defined(HANDLE_DEBUG)
					if ( refcov == windowsize )
					{
						double const conserror = NP.getTraceContainer().getAlignmentStatistics().getErrorRate();
						// std::cerr << "referror=" << referror << " conserror=" << conserror << std::endl;
						double const errdif = conserror - referror;
						ED [ errdif ] ++;
					}
					#endif

					uint64_t apos = astart;

					libmaus2::lcs::BaseConstants::step_type const * ta = NP.getTraceContainer().ta;
					libmaus2::lcs::BaseConstants::step_type const * te = NP.getTraceContainer().te;

					while ( ta != te )
					{
						uint64_t numins = 0;
						while ( ta != te && *ta == libmaus2::lcs::BaseConstants::STEP_INS )
						{
							++numins;
							++ta;
						}
						for ( uint64_t i = 0; i < numins; ++i )
						{
							assert ( cdata != consensus.second );
							PileElement PE(apos,(-static_cast<int64_t>(numins)) + static_cast<int64_t>(i),*(cdata++));
							PV.push(PVo,PE);
						}

						if ( ta != te )
						{
							assert ( *ta != libmaus2::lcs::BaseConstants::STEP_INS );

							switch ( *(ta++) )
							{
								case libmaus2::lcs::BaseConstants::STEP_MATCH:
								case libmaus2::lcs::BaseConstants::STEP_MISMATCH:
								{
									assert ( cdata != consensus.second );

									PileElement PE(apos++,0,*(cdata++));
									PV.push(PVo,PE);
									break;
								}
								case libmaus2::lcs::BaseConstants::STEP_DEL:
								{
									PileElement PE(apos++,0,'D');
									PV.push(PVo,PE);
									break;
								}
								default:
									break;
							}
						}
					}

					assert ( apos == aend );
				}

				if ( pathfailed )
					++failcount;
			}
			else
			{
				logstr << "insufficient depth " << MAo << std::endl;
				insufcount += 1;
			}

			#if defined(HANDLE_DEBUG)
			// how much to advance on reference or to reduce front soft clipping
			uint64_t const refadv = frontsoftclip >= W.offset(y) ? 0 : W.offset(y) - frontsoftclip;
			// advancement in trace operations
			std::pair<uint64_t,uint64_t> const advrefadv = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refadv);
			// advancement in string length
			std::pair<uint64_t,uint64_t> const slrefadv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advrefadv.second);
			// update front soft clipping
			frontsoftclip -= std::min(frontsoftclip,W.offset(y));
			// update trace start
			rta += advrefadv.second;
			// update reference pointer
			refu += slrefadv.first;
			#endif

			if ( printlogstr )
				err << logstr.str();
		}

		while ( !E.empty() )
		{
			upair UP = E.pop();
			activeset.erase(UP.second);

			uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
			trace_type::shared_ptr_type Ptrace = Mtraces[z];
			traceFreeList.put(Ptrace);
			Mtraces[z] = trace_type::shared_ptr_type();
		}

		while ( ivend < ovend )
		{
			uint64_t const id = Vend[ivend++].second;
			RC2.erase(id);
		}

		std::sort(PV.begin(),PV.begin()+PVo);

		if ( producefull && (ita != ite) )
		{
			uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita->aread()));

			uint64_t next = 0;
			uint64_t low = 0;
			uint64_t NPVo = 0;

			while ( low < PVo )
			{
				// look for consecutive pile elements
				uint64_t high = low+1;
				while ( high < PVo && PV[low].apos == PV[high].apos )
					++high;

				// fill up before elements
				for ( ; static_cast<int64_t>(next) < PV[low].apos; ++next )
					NPV.push(NPVo,PileElement(next,0,::tolower(ua[next])));

				// copy pile elements
				for ( uint64_t i = low; i < high; ++i )
					NPV.push(NPVo,PV[i]);

				next = PV[low].apos+1;
				low = high;
			}

			// fill up until end
			uint64_t const rl = RC.getReadLength(ita->aread());
			for ( ; next < rl; ++next )
				NPV.push(NPVo,PileElement(next,0,::tolower(ua[next])));

			PV.swap(NPV);
			PVo = NPVo;

			for ( uint64_t i = 1; i < PVo; ++i )
				assert ( PV[i-1].apos <= PV[i].apos );
		}

		uint64_t PVIo = 0;
		// std::vector< std::pair<uint64_t,uint64_t> > PVI;

		#if 0
		for ( uint64_t i = 0; i < PVo; ++i )
			std::cerr << PV[i].apos << "," << PV[i].apre << std::endl;
		#endif

		uint64_t il = 0;
		while ( il < PVo )
		{
			uint64_t ih = il+1;

			while (
				ih != PVo
				&&
				(PV[ih].apos - PV[ih-1].apos) <= 1
			)
				++ih;

			uint64_t const first = PV[il].apos;
			uint64_t const last = PV[ih-1].apos;

			if ( last-first >= 100 )
			{
				// std::cerr << "first=" << first << " last=" << last << std::endl;
				PVI.push(PVIo,std::pair<uint64_t,uint64_t>(il,ih));
			}

			il = ih;
		}

		for ( uint64_t z = 0; z < PVIo; ++z )
		{
			std::pair<uint64_t,uint64_t> const P = PVI[z];

			uint64_t const first = PV[P.first].apos;
			uint64_t const last = PV[P.second-1].apos;

			// std::vector<char> CO;
			uint64_t COo = 0;

			int64_t l = P.second;
			int64_t depth = -1;

			while ( l > static_cast<int64_t>(P.first) )
			{
				int64_t h = --l;

				while (
					l >= 0
					&&
					PV[l].apos == PV[h].apos
					&&
					PV[l].apre == PV[h].apre
				)
					--l;

				l += 1;

				for ( int64_t i = l; i <= h; ++i )
				{
					assert ( PV[i].apos == PV[l].apos );
					assert ( PV[i].apre == PV[l].apre );
				}

				uint64_t const ld = (h-l)+1;

				if ( PV[l].apre == 0 )
					depth = ld;

				std::pair<uint64_t,uint64_t> C[] =
				{
					std::pair<uint64_t,uint64_t>(0,'A'),
					std::pair<uint64_t,uint64_t>(0,'C'),
					std::pair<uint64_t,uint64_t>(0,'G'),
					std::pair<uint64_t,uint64_t>(0,'T'),
					std::pair<uint64_t,uint64_t>(0,'D'),
					std::pair<uint64_t,uint64_t>(0,'a'),
					std::pair<uint64_t,uint64_t>(0,'c'),
					std::pair<uint64_t,uint64_t>(0,'g'),
					std::pair<uint64_t,uint64_t>(0,'t'),
					std::pair<uint64_t,uint64_t>(0,0)
				};

				#if defined(PV_DEBUG)
				logstr << "PV " << PV[l].apos << "," << PV[l].apre << " " << depth << " ";
				#endif

				for ( int64_t i = l; i <= h; ++i )
				{
					#if defined(PV_DEBUG)
					logstr << PV[i].sym;
					#endif

					switch ( PV[i].sym )
					{
						case 'A': C[0].first++; break;
						case 'C': C[1].first++; break;
						case 'G': C[2].first++; break;
						case 'T': C[3].first++; break;
						case 'D': C[4].first++; break;
						case 'a': C[5].first++; break;
						case 'c': C[6].first++; break;
						case 'g': C[7].first++; break;
						case 't': C[8].first++; break;
						default: assert(0); break;
					}
				}

				for ( int64_t i = ld; i < depth; ++i )
					C[4].first++;

				std::sort(&C[0],&C[sizeof(C)/sizeof(C[0])],std::greater<std::pair<uint64_t,uint64_t> >());

				#if defined(PV_DEBUG)
				err.put(' ');
				for ( uint64_t i = 0; C[i].first; ++i )
					logstr << "(" << static_cast<char>(C[i].second) << "," << C[i].first << ")";
				err << std::endl;
				#endif

				if ( C[0].first && C[0].second != 'D' )
					CO.push(COo,C[0].second);
			}

			std::reverse(CO.begin(),CO.begin()+COo);

			if ( producefull || COo >= minlen )
			{
				out << '>' << (aid+1) << '/' << wellcounter++ << '/' << first << '_' << first + COo << " A=[" << first << "," << last << "]" << "\n";
				char const * zp = CO.begin();
				char const * ze = zp + COo;
				uint64_t const linewidth = 80;
				while ( zp != ze )
				{
					uint64_t const r = ze-zp;
					uint64_t const tp = std::min(r,linewidth);
					out.write(zp,tp);
					out.put('\n');
					zp += tp;
				}
			}

			#if 0
			{
				std::string const SCO(CO.begin(),CO.begin()+COo);

				// consensus text
				uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(SCO.c_str());
				uint8_t const * uSCOe = uSCO + SCO.size();
				// last read position for consensus
				uint64_t const last = PV[P.second-1].apos;

				std::string const conspart(uSCO, uSCOe);
				std::string const readpart(RC.getForwardRead(ita[0].aread())+first,RC.getForwardRead(ita[0].aread())+last+1);

				libmaus2::lcs::NP NPcons;

				NPcons.np(readpart.begin(),readpart.end(),conspart.begin(),conspart.end());

				int64_t rpos = first;

				libmaus2::lcs::AlignmentTraceContainer::step_type const * consta = NPcons.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * conste = NPcons.te;

				if ( rpos % advancesize != 0 )
				{
					int64_t const skip = advancesize - (rpos % advancesize);

					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,skip);
					assert ( Pcons.first == skip );
					consta += Pcons.second;

					rpos += skip;
				}

				assert ( consta == conste || (rpos % advancesize == 0) );

				while ( libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize).first == windowsize )
				{
					// std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize);
					//double const conserate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(consta,consta+Pcons.second).getErrorRate();
					consta += libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,advancesize).second;
				}
			}
			#endif

			#if defined(HANDLE_DEBUG)
			{
				uint8_t const * uref = reinterpret_cast<uint8_t const *>(sreftext.c_str());
				uint8_t const * urefe = uref + sreftext.size();

				// err << SCO << std::endl;

				uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(SCO.c_str());
				uint8_t const * uSCOe = uSCO + SCO.size();

				libmaus2::lcs::SuffixArrayLCS::LCSResult const lcsres = libmaus2::lcs::SuffixArrayLCS::lcsmin(std::string(uref,urefe),std::string(uSCO,uSCOe));
				libmaus2::lcs::NNP nnp;
				libmaus2::lcs::NNPTraceContainer nnptrace;
				// align ref and consensus
				libmaus2::lcs::NNPAlignResult const nnpres = nnp.align(uref,urefe,lcsres.maxpos_a,uSCO,uSCOe,lcsres.maxpos_b,nnptrace);
				//nnptrace.printTraceLines(err,uref + nnpres.abpos,uSCO + nnpres.bbpos);

				libmaus2::lcs::AlignmentTraceContainer ATC;
				nnptrace.computeTrace(ATC);

				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					uref+nnpres.abpos,
					nnpres.aepos-nnpres.abpos,
					uSCO+nnpres.bbpos,
					nnpres.bepos-nnpres.bbpos,
					80,
					ATC.ta,
					ATC.te
				);

				libmaus2::lcs::AlignmentStatistics const AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(ATC.ta,ATC.te);

				// compare alignments ref/read to cons/read
				std::string const refpart(uref + nnpres.abpos,uref + nnpres.aepos);
				std::string const conspart(uSCO + nnpres.bbpos, uSCO + nnpres.bepos);
				std::string const readpart(RC.getForwardRead(ita[0].aread())+first,RC.getForwardRead(ita[0].aread())+last+1);

				libmaus2::lcs::NP NPref;
				libmaus2::lcs::NP NPcons;

				NPref.np(readpart.begin(),readpart.end(),refpart.begin(),refpart.end());
				NPcons.np(readpart.begin(),readpart.end(),conspart.begin(),conspart.end());

				#if 0
				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					readpart.begin(),
					readpart.size(),
					refpart.begin(),
					refpart.size(),
					80,
					NPref.ta,
					NPref.te
				);

				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					err,
					readpart.begin(),
					readpart.size(),
					conspart.begin(),
					conspart.size(),
					80,
					NPcons.ta,
					NPcons.te
				);
				#endif

				int64_t rpos = first;

				libmaus2::lcs::AlignmentTraceContainer::step_type const * refta = NPref.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * refte = NPref.te;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * consta = NPcons.ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * conste = NPcons.te;

				if ( rpos % advancesize != 0 )
				{
					int64_t const skip = advancesize - (rpos % advancesize);

					std::pair<uint64_t,uint64_t> const Pref = libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,skip);
					assert ( static_cast<int64_t>(Pref.first) == skip );
					refta += Pref.second;

					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,skip);
					assert ( static_cast<int64_t>(Pcons.first) == skip );
					consta += Pcons.second;

					rpos += skip;
				}

				while ( libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,windowsize).first == windowsize )
				{
					std::pair<uint64_t,uint64_t> const Pref = libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,windowsize);
					std::pair<uint64_t,uint64_t> const Pcons = libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,windowsize);

					double const referate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(refta,refta+Pref.second).getErrorRate();
					double const conserate = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(consta,consta+Pcons.second).getErrorRate();

					if ( std::abs(conserate-referate) >= 1e-6 )
						err << rpos << " (" << referate << "," << conserate << "," << referate-conserate << ")\n";

					ED[referate - conserate]++;

					refta += libmaus2::lcs::AlignmentTraceContainer::advanceA(refta,refte,advancesize).second;
					consta += libmaus2::lcs::AlignmentTraceContainer::advanceA(consta,conste,advancesize).second;
					rpos += advancesize;
				}

				err << "read id " << ita[0].aread() << " " << AS << " " << nnpres << std::endl;

				// err << SCO << std::endl;
			}
			#endif
		}

		// err << "read id " << ita[0].aread << " " << AS << " failcount=" << failcount << " insufcount=" << insufcount << " unbroken=" << Vsufsub << " " << nnpres << " pvdense=" << pvdense << std::endl;

		// return traces
		for ( uint64_t i = 0; i < nintv; ++i )
			assert ( ! Mtraces[i] );
			// if ( Mtraces[i] )
				//traceFreeList.put(Mtraces[i]);

		#if defined(HANDLE_DEBUG)
		// return BAM trace for ground truth
		traceFreeList.put(Pbamtrace);
		#endif

		// DG.printSize(err);

		err << "[V] read id " << ita[0].aread() << " time " << handlertc << std::endl;
	}

	void operator()(
		std::ostream & out,
		std::ostream & err,
		// overlaps
		std::string const & lasfn,
		libmaus2::dazzler::align::BinIndexDecoder const & binindex,
		uint64_t const aread
	)
	{
		uint64_t const rl = DecodedReadContainer(readDataFreeList,readDecoderFreeList).getReadLength(aread);

		Windows W(rl,advancesize,windowsize);

		uint64_t const numwindows = W.size();

		#if 0
		uint64_t const packsperthread = 8;
		uint64_t const tnumpacks = numthreads * packsperthread;
		uint64_t const windowsperpack = (numwindows + tnumpacks - 1)/tnumpacks;
		uint64_t const numpacks = windowsperpack ? (numwindows + windowsperpack - 1)/windowsperpack : 0;
		#endif

		std::string const emptystr;

		std::vector < uint64_t > VPVo(numthreads,0);
		int volatile loopfailed = 0;

		uint64_t volatile windowsfinished = 0;
		uint64_t const windowmod = 256;

		DecodedReadContainer RC (readDataFreeList,readDecoderFreeList );
		uint8_t const * g_ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(aread));

		struct VolatilePair
		{
			uint64_t volatile first;
			uint64_t volatile second;

			VolatilePair() {}
			VolatilePair(uint64_t const ra, uint64_t const rb) : first(ra), second(rb) {}
		};

		struct TodoList
		{
			uint64_t const numthreads;
			uint64_t const numwindows;

			libmaus2::autoarray::AutoArray < VolatilePair > A;
			libmaus2::autoarray::AutoArray < libmaus2::parallel::PosixSpinLock::unique_ptr_type > L;

			TodoList(uint64_t const rnumthreads, uint64_t const rnumwindows)
			: numthreads(rnumthreads), numwindows(rnumwindows), A(numthreads), L(numthreads)
			{
				uint64_t const windowsperthread = (numwindows + numthreads - 1)/numthreads;

				for ( uint64_t i = 0; i < numthreads; ++i )
				{
					uint64_t const low = std::min(i * windowsperthread,numwindows);
					uint64_t const high = std::min(low + windowsperthread, numwindows);
					A[i] = VolatilePair(low,high);

					libmaus2::parallel::PosixSpinLock::unique_ptr_type TL(new libmaus2::parallel::PosixSpinLock);
					L[i] = UNIQUE_PTR_MOVE(TL);
				}
			}

			bool getNext(uint64_t const tid, uint64_t & idlow, uint64_t & idhigh)
			{
				enum state { normal, steal };
				state s = normal;

				while ( true )
				{
					switch ( s )
					{
						case normal:
						{
							bool lok = false;

							L[tid]->lock();

							if ( A[tid].second != A[tid].first )
							{
								assert ( A[tid].first < A[tid].second );

								idlow  = A[tid].first++;
								idhigh = A[tid].second;
								lok = true;
							}

							L[tid]->unlock();

							if ( lok )
							{
								return true;
							}
							else
							{
								s = steal;
							}

							break;
						}
						case steal:
						{
							for ( uint64_t t = 0; t < numthreads; ++t )
							{
								L[t]->lock();
							}

							int64_t m = -1;
							uint64_t s = 0;
							for ( uint64_t t = 0; t < numthreads; ++t )
								if ( A[t].second != A[t].first )
								{
									uint64_t const l = A[t].second - A[t].first;

									if ( l > s )
									{
										s = l;
										m = t;
									}
								}

							bool lok = false;

							if ( s )
							{
								assert ( m >= 0 );
								assert ( static_cast<uint64_t>(m) != tid );

								uint64_t const s2 = (s+1)/2;
								uint64_t const low   = A[m].first;
								uint64_t const split = A[m].second - s2;
								uint64_t const high  = A[m].second;

								assert ( s2 > 0 );
								assert ( split >= low );
								assert ( high > split );

								A[m].first    = low;
								A[m].second   = split;
								A[tid].first  = split;
								A[tid].second = high;

								idlow  = A[tid].first++;
								idhigh = A[tid].second;
								lok = true;
							}

							for ( uint64_t t = 0; t < numthreads; ++t )
							{
								L[t]->unlock();
							}

							if ( lok )
							{
								return true;
							}
							else
							{
								return false;
							}
						}
					}
				}
			}
		};

		TodoList TL(numthreads,numwindows);

		#if defined(_OPENMP)
		#pragma omp parallel num_threads(numthreads)
		#endif
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			try
			{

				uint64_t PVo = VPVo[tid];

				ThreadContext & threadcontext = *(Pthreadcontext.at(tid));
				libmaus2::lcs::Aligner & NP = threadcontext.NP;
				libmaus2::autoarray::AutoArray < PileElement > & PV = threadcontext.PV;
				libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > & MA = threadcontext.MA;
				libmaus2::autoarray::AutoArray<trace_type::shared_ptr_type> & Mtraces = threadcontext.Mtraces;
				libmaus2::autoarray::AutoArray<double> & VVVV = threadcontext.VVVV;
				libmaus2::autoarray::AutoArray<uint64_t> & Vdist = threadcontext.Vdist;
				libmaus2::autoarray::AutoArray<std::pair<uint64_t,uint64_t> > & VPdist = threadcontext.VPdist;

				libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & E = threadcontext.E;
				std::map < uint64_t, ActiveElement > & activeset = threadcontext.activeset;

				DecodedReadContainer RC2(readDataFreeList,readDecoderFreeList2);

				uint64_t nextwindow = std::numeric_limits<uint64_t>::max();
				uint64_t packstart = std::numeric_limits<uint64_t>::max();
				uint64_t packend   = std::numeric_limits<uint64_t>::max();

				#if 0
				uint64_t const wlow = pack * windowsperpack;
				uint64_t const whigh = std::min(wlow + windowsperpack, numwindows);
				#endif

				#if 0
				uint64_t const packstart = W[wlow   ].first;
				uint64_t const packend   = W[whigh-1].second;
				#endif


				#if 0
				{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] processing pack=" << pack << "/" << numpacks << " wlow=" << wlow << " whigh=" << whigh << " packstart=" << packstart << " packend=" << packend << std::endl;
				}
				#endif

				libmaus2::dazzler::align::LasRangeDecoder::unique_ptr_type PLRD;

				#if 0
				libmaus2::dazzler::align::LasRangeDecoder LRD(lasfn,binindex);
				LRD.setup(aread,rl,packstart,packend);
				#endif

				assert ( E.empty() );
				assert ( activeset.empty() );

				std::ostringstream logstr;
				uint64_t z = 0;
				uint64_t failcount = 0;
				uint64_t insufcount = 0;
				uint64_t tllow;
				uint64_t tlhigh;

				// for ( uint64_t y = wlow; y < whigh; ++y )
				while ( TL.getNext(tid,tllow,tlhigh) )
				{
					#if 0
					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "[V] thread " << tid << " tllow=" << tllow << std::endl;
					}
					#endif

					if ( tllow != nextwindow )
					{
						packstart = W[tllow].first;
						packend = W[tlhigh-1].second;

						#if 0
						{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "[V] thread " << tid << " switching to window interval [" << tllow << "," << tlhigh << ") base interval [" << packstart << "," << packend << ")" << std::endl;
						}
						#endif

						PLRD.reset();
						libmaus2::dazzler::align::LasRangeDecoder::unique_ptr_type TLRD(new libmaus2::dazzler::align::LasRangeDecoder(lasfn,binindex));
						PLRD = UNIQUE_PTR_MOVE(TLRD);
						PLRD->setup(aread,rl,packstart,packend);

						while ( !E.empty() )
						{
							std::pair<uint64_t,uint64_t> UP = E.pop();
							activeset.erase(UP.second);

							uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
							trace_type::shared_ptr_type Ptrace = Mtraces[z];
							traceFreeList.put(Ptrace);
							Mtraces[z] = trace_type::shared_ptr_type();
						}

						assert ( E.empty() );
						assert ( activeset.empty() );

						z = 0;

						nextwindow = tllow;
					}

					uint64_t const y = tllow;

					nextwindow += 1;

					bool printlogstr = false;

					logstr.str(emptystr);
					logstr.clear();

					// start of window on A
					uint64_t const astart = W[y].first;
					uint64_t const aend   = W[y].second;

					libmaus2::dazzler::align::Overlap OVL;
					while ( PLRD->peekNext(OVL) && static_cast<int64_t>(astart) >= OVL.path.abpos )
					{
						#if 0
						if ( y == 11896 )
						{
							libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
							std::cerr << "[V] thread " << tid << " tllow=" << tllow << " " << OVL << std::endl;
						}
						#endif

						PLRD->getNext(OVL);

						if ( OVL.path.aepos >= static_cast<int64_t>(astart) )
						{
							uint64_t const lz = z;

							uint8_t const * l_ub = reinterpret_cast<uint8_t const *>(OVL.isInverse() ? RC2.getReverseComplementRead(OVL.bread) : RC2.getForwardRead(OVL.bread));
							trace_type::shared_ptr_type Ptrace = traceFreeList.get();
							OVL.computeTrace(g_ua,l_ub,tspace,*Ptrace,NP);
							Mtraces.push(z,Ptrace);

							// offset from start of window
							uint64_t const aoff = astart - OVL.path.abpos;
							uint8_t const * ua = g_ua + astart;

							#if 0
							std::cerr << "[V] adding z=" << lz << " " << OVL << std::endl;
							std::cerr << "[V] aoff=" << aoff << std::endl;
							#endif

							libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = 0;
							libmaus2::lcs::AlignmentTraceContainer::step_type const * te = 0;
							uint64_t uboff;

							{
								// get trace
								trace_type::shared_ptr_type Ptrace = Mtraces[lz];
								ta = Ptrace->ta;
								te = Ptrace->te;

								// see how many operations there are up to start of region we are interested in
								std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(ta,te,aoff);

								bool const ok = adv.first == aoff;
								if ( ! ok )
								{
									std::cerr << "adv.first=" << adv.first << " aoff=" << aoff << " astart=" << astart << std::endl;
									std::cerr << OVL << std::endl;
									assert ( adv.first == aoff );
								}

								uboff = OVL.path.bbpos + libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ta,ta+adv.second).second;

								// advance in trace
								ta += adv.second;
							}

							uint8_t const * ub = l_ub + uboff;

							// add to active set
							uint64_t const escore = static_cast<uint64_t>(OVL.getErrorRate() * std::numeric_limits<uint32_t>::max());
							assert ( escore <= std::numeric_limits<uint32_t>::max() );
							uint64_t const eindex = (escore<<32) | lz;

							activeset[eindex] = ActiveElement(ua,ub,ta,te,uboff,OVL.getErrorRate());

							// add to erase list
							E.pushBump(std::pair<uint64_t,uint64_t>(OVL.path.aepos,eindex));
						}
					}

					// cleanup
					while ( (! (E.empty())) && E.top().first < aend )
					{
						std::pair<uint64_t,uint64_t> UP = E.pop();
						uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
						trace_type::shared_ptr_type Ptrace = Mtraces[z];
						traceFreeList.put(Ptrace);
						Mtraces[z] = trace_type::shared_ptr_type();
						activeset.erase(UP.second);
					}

					uint64_t MAo = 0;

					uint8_t const * w_ua = activeset.size() ? activeset.begin()->second.ua : 0;

					// iterate over active read list
					for ( std::map<uint64_t,ActiveElement>::iterator s_ita = activeset.begin(); s_ita != activeset.end(); ++s_ita )
					{
						// active element
						ActiveElement & AE = s_ita->second;

						// read id
						// uint64_t const key = s_ita->first & std::numeric_limits<uint32_t>::max();

						uint64_t bwindowsize;
						uint64_t badvancesize;

						{
							// get end of region in trace
							std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,windowsize);

							// bool const ok = adv.first == windowsize;

							assert ( adv.first == windowsize );
							// string length
							std::pair<uint64_t,uint64_t> const sl = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+adv.second);

							bwindowsize = sl.second;

							// compute how much we advance to set up for next window
							std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,W.offset(y)); // advancesize
							// in string length used on A and B read
							std::pair<uint64_t,uint64_t> const sladv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+advadv.second);
							// update trace pointer
							AE.ta += advadv.second;

							badvancesize = sladv.second;
						}

						// push A read if this is the first instance of this loop
						if ( MAo < maxalign )
						{
							// push B read
							MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ub,bwindowsize));
						}

						// update active element by advancing to next window
						AE.ua += W.offset(y); // advancesize
						AE.ub += badvancesize;
						AE.uboff += badvancesize;
					}

					int64_t maxvprodindex = -1;
					if ( MAo )
					{
						// maximum position on first read as minimum of support
						int64_t minSupLen = static_cast<int64_t>(MA[0].second)-1;
						// support maximum
						int64_t maxSupLen = minSupLen;
						for ( uint64_t j = 1; j < MAo; ++j )
						{
							int64_t const lastpos = static_cast<int64_t>(MA[j].second)-1;
							minSupLen = std::min(minSupLen,lastpos);
							maxSupLen = std::max(maxSupLen,lastpos);
						}
						if ( minSupLen < 0 )
							minSupLen = 0;
						if ( maxSupLen < 0 )
							maxSupLen = 0;

						uint64_t supStart = offsetLikely.getSupportLow(minSupLen);
						uint64_t supEnd   = offsetLikely.getSupportHigh(maxSupLen);

						assert ( (supStart == 0) || (minSupLen >= static_cast<int64_t>(offsetLikely.DPnorm[supStart-1].size())) );
						assert ( (supEnd == offsetLikely.DPnorm.size()) ||  (maxSupLen < static_cast<int64_t>(offsetLikely.DPnorm[supEnd].firstsign)) );

						double maxval = std::numeric_limits<double>::min();
						// iterate over ref positions in support
						for ( uint64_t i = supStart; i < supEnd; ++i )
						{
							DotProduct const & DP = offsetLikely.DPnorm[i];

							double vprod = 1.0;

							for ( uint64_t j = 0; j < MAo; ++j )
							{
								uint64_t const len = MA[j].second;
								if ( len )
								{
									uint64_t const lastpos = len-1;

									vprod *= DP[lastpos];
								}
							}

							if ( vprod > maxval )
							{
								maxval = vprod;
								maxvprodindex = i;
							}
						}
					}

					// length estimation via conditional probability failed, try comparing density distributions
					if ( maxvprodindex == -1 )
					{
						int64_t maxoff = -1;
						double maxoffv = std::numeric_limits<double>::min();


						uint64_t Vdisto = 0;
						uint64_t VPdisto = 0;
						for ( uint64_t i = 0; i < MAo; ++i )
							Vdist.push(Vdisto,MA[i].second);
						std::sort(Vdist.begin(),Vdist.begin()+Vdisto);
						{
							uint64_t low = 0;
							while ( low < Vdisto )
							{
								uint64_t high = low+1;
								while ( high < Vdisto && Vdist[high] == Vdist[low] )
									++high;

								VPdist.push(VPdisto,std::pair<uint64_t,uint64_t>(Vdist[low],high-low));

								low = high;
							}
						}

						uint64_t VVVVo = 0;
						for ( uint64_t i = 0; i < VPdisto; ++i )
						{
							while ( ! ( VPdist[i].first < VVVVo ) )
								VVVV.push(VVVVo,0);
							VVVV[VPdist[i].first] = VPdist[i].second-1;
						}

						for ( uint64_t i = 0; i < offsetLikely.DPnormSquare.size(); ++i )
						{
							double const v = offsetLikely.DPnormSquare[i].dotproduct(VVVV.begin(),VVVVo);
							if ( v > maxoffv )
							{
								maxoff = i;
								maxoffv = v;
							}
							#if 0
							if ( v > 1e-3 )
								logstr << "prod " << (i+1) << " = " << v << " slref=" << slref.first << std::endl;
							#endif
						}

						if ( maxoff != -1 && maxoffv >= 1e-3 )
						{
							maxvprodindex = maxoff;
							// logstr << "prod " << (maxvprodindex+1) << " = " << maxoffv << " slref=" << slref.first << std::endl;
						}
					}

					#if 0
					if ( y+1 == W.size() )
					{
						std::cerr << "window " << astart << "," << aend << " maxaepos=" << maxaepos << " MAo=" << MAo << std::endl;
					}
					#endif

					if (
						(MAo >= minwindowcov)
						#if 0
						&&
						y == 1624
						#endif
					)
					{
						int64_t const elength = maxvprodindex+1;

						logstr << "elength " << elength << std::endl;

						#if 0
						for ( uint64_t i = 0; i < MAo; ++i )
						{
							err << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
						}
						#endif

						bool pathfailed = true;
						int64_t filterfreq = -1;

						uint64_t minindex = 0;
						uint64_t minrate = eminrate;
						DebruijnGraphInterface * minDG = 0;

						for ( uint64_t adgi = 0; adgi < threadcontext.ADG.size(); ++adgi )
						{
							DebruijnGraphInterface & DG = *(threadcontext.ADG[adgi]);

							// logstr << "refseq " << std::string(refu,refu+slref.first) << std::endl;
							filterfreq = maxfilterfreq;

							for ( ; filterfreq >= minfilterfreq ; --filterfreq )
							{
								#if defined(HANDLE_TIME)
								libmaus2::timing::RealTimeClock rtc;
								#endif

								// set up debruijn graph
								#if defined(HANDLE_TIME)
								rtc.start();
								#endif
								DG.setup(MA.begin(), MAo);
								#if defined(HANDLE_TIME)
								err << "setup " << rtc.getElapsedSeconds() << std::endl;
								#endif

								// filter by frequency (need at least 2)
								#if defined(HANDLE_TIME)
								rtc.start();
								#endif
								DG.filterFreq(std::max(filterfreq,static_cast<int64_t>(1)),MAo);
								#if defined(HANDLE_TIME)
								err << "filterfreq " << rtc.getElapsedSeconds() << std::endl;
								#endif

								#if defined(HANDLE_TIME)
								rtc.start();
								#endif
								DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
								#if defined(HANDLE_TIME)
								err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
								#endif

								if ( filterfreq == 0 )
								{
									logstr << "*** trying to fill in gaps" << std::endl;

									#if defined(HANDLE_TIME)
									rtc.start();
									#endif
									DG.getLevelSuccessors(2);
									#if defined(HANDLE_TIME)
									err << "getLevelSuccessors " << rtc.getElapsedSeconds() << std::endl;
									#endif

									#if defined(HANDLE_TIME)
									rtc.start();
									#endif
									DG.setupNodes();
									#if defined(HANDLE_TIME)
									err << "setupNodes " << rtc.getElapsedSeconds() << std::endl;
									#endif

									#if defined(HANDLE_TIME)
									rtc.start();
									#endif
									DG.setupAddHeap(MAo);
									#if defined(HANDLE_TIME)
									err << "setupAddHeap " << rtc.getElapsedSeconds() << std::endl;
									#endif

									#if defined(HANDLE_TIME)
									rtc.start();
									#endif
									DG.computeFeasibleKmerPositions(offsetLikely, 1e-3 /* thres */);
									#if defined(HANDLE_TIME)
									err << "computeFeasibleKmerPositions " << rtc.getElapsedSeconds() << std::endl;
									#endif
								}

								uint64_t mintry = 0;
								uint64_t const maxtries = 3;
								bool lconsok = false;

								do
								{
									#if defined(HANDLE_TIME)
									rtc.start();
									#endif

									#if defined(HANDLE_DEEP_DEBUG)
									DG.toDot(err);
									#endif

									bool const consok = DG.traverse(elength - 4,elength + 4,MA.begin(),MAo,16 /* maxfronpath */,16 /* maxfullpath */);
									#if defined(HANDLE_TIME)
									err << "traverse " << rtc.getElapsedSeconds() << std::endl;
									#endif

									if ( consok )
									{
										std::pair<uint64_t,uint64_t> const MR = DG.checkCandidatesU(MA.begin(),MAo);

										if ( MR.second < minrate )
										{
											lconsok = true;

											minrate = MR.second;
											minindex = MR.first;
											minDG = &DG;
										}
										else if ( minDG )
										{
											lconsok = true;
										}
										break;
									}
									else
									{
										if ( ++mintry >= maxtries )
											break;
									}
								} while (
									DG.addNextFromHeap(&logstr)
								);

								//bool printdebuginfo = false;

								if ( !lconsok )
								{
									logstr << "[filterfreq=" << DG.getKmerSize() << "," << filterfreq << "] failed " << MAo << std::endl;
									// printdebuginfo = true;

									if ( verbose >= 2 )
										printlogstr = true;
								}
								else
								{
									pathfailed = false;
									break;
								}
							}
						}

						if ( ! pathfailed )
						{
							assert ( minDG );

							std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(minindex);

							uint8_t const * aread = w_ua; // MA[0].first;
							uint64_t const areadlen = windowsize; // MA[0].second;
							uint8_t const * cdata = consensus.first;
							uint64_t const clen = consensus.second - consensus.first;

							NP.align(aread,areadlen,cdata,clen);

							uint64_t apos = astart;

							libmaus2::lcs::BaseConstants::step_type const * ta = NP.getTraceContainer().ta;
							libmaus2::lcs::BaseConstants::step_type const * te = NP.getTraceContainer().te;

							while ( ta != te )
							{
								uint64_t numins = 0;
								while ( ta != te && *ta == libmaus2::lcs::BaseConstants::STEP_INS )
								{
									++numins;
									++ta;
								}
								for ( uint64_t i = 0; i < numins; ++i )
								{
									assert ( cdata != consensus.second );
									PileElement PE(apos,(-static_cast<int64_t>(numins)) + static_cast<int64_t>(i),*(cdata++));
									PV.push(PVo,PE);
								}

								if ( ta != te )
								{
									assert ( *ta != libmaus2::lcs::BaseConstants::STEP_INS );

									switch ( *(ta++) )
									{
										case libmaus2::lcs::BaseConstants::STEP_MATCH:
										case libmaus2::lcs::BaseConstants::STEP_MISMATCH:
										{
											assert ( cdata != consensus.second );

											PileElement PE(apos++,0,*(cdata++));
											PV.push(PVo,PE);
											break;
										}
										case libmaus2::lcs::BaseConstants::STEP_DEL:
										{
											PileElement PE(apos++,0,'D');
											PV.push(PVo,PE);
											break;
										}
										default:
											break;
									}
								}
							}

							assert ( apos == aend );
						}

						if ( pathfailed )
							++failcount;
					}
					else
					{
						logstr << "insufficient depth " << MAo << std::endl;
						insufcount += 1;
					}

					if ( printlogstr )
					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						err << logstr.str();
					}


					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						windowsfinished += 1;

						if ( windowsfinished % windowmod == 0 )
							std::cerr << "[V] aread " << aread << " windows finished " << windowsfinished << "/" << numwindows << std::endl;
					}
				}

				while ( !E.empty() )
				{
					std::pair<uint64_t,uint64_t> UP = E.pop();
					activeset.erase(UP.second);

					uint64_t const z = UP.second & ((static_cast<uint64_t>(1ull)<<32)-1);
					trace_type::shared_ptr_type Ptrace = Mtraces[z];
					traceFreeList.put(Ptrace);
					Mtraces[z] = trace_type::shared_ptr_type();
				}

				assert ( E.empty() );
				assert ( activeset.empty() );

				#if 0
				for ( uint64_t i = 0; i < PVo; ++i )
				{
					std::cerr << "PV[" << i << "]=" << PV[i].toString() << std::endl;
				}
				#endif

				VPVo.at(tid) = PVo;
			}
			catch(std::exception const & ex)
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[E] tid=" << tid << " " << ex.what() << std::endl;
				loopfailed = 1;
			}
		}

		if ( windowsfinished % windowmod != 0 )
			std::cerr << "[V] aread " << aread << " windows finished " << windowsfinished << "/" << numwindows << std::endl;

		if ( loopfailed )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] parallel loop failed in HandleContext" << std::endl;
			lme.finish();
			throw lme;
		}

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
		#endif
		for ( uint64_t t = 0; t < numthreads; ++t )
		{
			ThreadContext & threadcontext = *(Pthreadcontext.at(t));
			libmaus2::autoarray::AutoArray < PileElement > & PV = threadcontext.PV;
			uint64_t const PVo = VPVo.at(t);
			std::sort(PV.begin(),PV.begin()+PVo);
		}

		std::vector<uint64_t> VM;

		while ( VPVo.size() > 1 )
		{
			// std::cerr << "[V] PV merge round with " << VPVo.size() << " vectors" << std::endl;

			uint64_t o = 0;

			for ( uint64_t j = 0; j < VPVo.size(); j += 2, ++o )
				if ( j+1 < VPVo.size() )
				{
					libmaus2::autoarray::AutoArray < PileElement > & PV0 = Pthreadcontext[j+0]->PV;
					libmaus2::autoarray::AutoArray < PileElement > & PV1 = Pthreadcontext[j+1]->PV;
					uint64_t const s = VPVo.at(j+0)+VPVo.at(j+1);
					libmaus2::autoarray::AutoArray < PileElement > PV2(s,false);
					libmaus2::sorting::ParallelStableSort::parallelMerge<PileElement *,std::less<PileElement> >(
						PV0.begin(),PV0.begin()+VPVo.at(j+0),
						PV1.begin(),PV1.begin()+VPVo.at(j+1),
						PV2.begin()
					);

					VPVo.at(o) = s;
					Pthreadcontext[o]->PV.swap(PV2);
				}
				else
				{
					libmaus2::autoarray::AutoArray < PileElement > & PV0 = Pthreadcontext[j+0]->PV;
					uint64_t const s = VPVo.at(j+0);
					VPVo.at(o) = s;
					Pthreadcontext[o]->PV.swap(PV0);
				}

			VPVo.resize(o);
		}

		assert ( VPVo.size() == 1 );
		uint64_t PVo = VPVo.at(0);
		libmaus2::autoarray::AutoArray < PileElement > & PV = Pthreadcontext[0]->PV;
		libmaus2::autoarray::AutoArray < PileElement > & NPV = Pthreadcontext[0]->NPV;

		if ( producefull )
		{
			uint64_t next = 0;
			uint64_t low = 0;
			uint64_t NPVo = 0;

			while ( low < PVo )
			{
				// look for consecutive pile elements
				uint64_t high = low+1;
				while ( high < PVo && PV[low].apos == PV[high].apos )
					++high;

				// fill up before elements
				for ( ; static_cast<int64_t>(next) < PV[low].apos; ++next )
					NPV.push(NPVo,PileElement(next,0,::tolower(g_ua[next])));

				// copy pile elements
				for ( uint64_t i = low; i < high; ++i )
					NPV.push(NPVo,PV[i]);

				next = PV[low].apos+1;
				low = high;
			}

			// fill up until end
			uint64_t const rl = RC.getReadLength(aread);
			for ( ; next < rl; ++next )
				NPV.push(NPVo,PileElement(next,0,::tolower(g_ua[next])));

			PV.swap(NPV);
			PVo = NPVo;

			for ( uint64_t i = 1; i < PVo; ++i )
				assert ( PV[i-1].apos <= PV[i].apos );
		}

		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > & PVI = Pthreadcontext[0]->PVI;
		uint64_t PVIo = 0;
		libmaus2::autoarray::AutoArray<char> & CO = Pthreadcontext[0]->CO;
		// std::vector< std::pair<uint64_t,uint64_t> > PVI;

		#if 0
		for ( uint64_t i = 0; i < PVo; ++i )
			std::cerr << PV[i].apos << "," << PV[i].apre << std::endl;
		#endif

		uint64_t il = 0;
		while ( il < PVo )
		{
			uint64_t ih = il+1;

			while (
				ih != PVo
				&&
				(PV[ih].apos - PV[ih-1].apos) <= 1
			)
				++ih;

			uint64_t const first = PV[il].apos;
			uint64_t const last = PV[ih-1].apos;

			if ( last-first >= 100 )
			{
				// std::cerr << "first=" << first << " last=" << last << std::endl;
				PVI.push(PVIo,std::pair<uint64_t,uint64_t>(il,ih));
			}

			il = ih;
		}

		for ( uint64_t z = 0; z < PVIo; ++z )
		{
			std::pair<uint64_t,uint64_t> const P = PVI[z];

			uint64_t const first = PV[P.first].apos;
			uint64_t const last = PV[P.second-1].apos;

			// std::vector<char> CO;
			uint64_t COo = 0;

			int64_t l = P.second;
			int64_t depth = -1;

			while ( l > static_cast<int64_t>(P.first) )
			{
				int64_t h = --l;

				while (
					l >= 0
					&&
					PV[l].apos == PV[h].apos
					&&
					PV[l].apre == PV[h].apre
				)
					--l;

				l += 1;

				for ( int64_t i = l; i <= h; ++i )
				{
					assert ( PV[i].apos == PV[l].apos );
					assert ( PV[i].apre == PV[l].apre );
				}

				uint64_t const ld = (h-l)+1;

				if ( PV[l].apre == 0 )
					depth = ld;

				std::pair<uint64_t,uint64_t> C[] =
				{
					std::pair<uint64_t,uint64_t>(0,'A'),
					std::pair<uint64_t,uint64_t>(0,'C'),
					std::pair<uint64_t,uint64_t>(0,'G'),
					std::pair<uint64_t,uint64_t>(0,'T'),
					std::pair<uint64_t,uint64_t>(0,'D'),
					std::pair<uint64_t,uint64_t>(0,'a'),
					std::pair<uint64_t,uint64_t>(0,'c'),
					std::pair<uint64_t,uint64_t>(0,'g'),
					std::pair<uint64_t,uint64_t>(0,'t'),
					std::pair<uint64_t,uint64_t>(0,0)
				};

				#if defined(PV_DEBUG)
				logstr << "PV " << PV[l].apos << "," << PV[l].apre << " " << depth << " ";
				#endif

				for ( int64_t i = l; i <= h; ++i )
				{
					#if defined(PV_DEBUG)
					logstr << PV[i].sym;
					#endif

					switch ( PV[i].sym )
					{
						case 'A': C[0].first++; break;
						case 'C': C[1].first++; break;
						case 'G': C[2].first++; break;
						case 'T': C[3].first++; break;
						case 'D': C[4].first++; break;
						case 'a': C[5].first++; break;
						case 'c': C[6].first++; break;
						case 'g': C[7].first++; break;
						case 't': C[8].first++; break;
						default: assert(0); break;
					}
				}

				for ( int64_t i = ld; i < depth; ++i )
					C[4].first++;

				std::sort(&C[0],&C[sizeof(C)/sizeof(C[0])],std::greater<std::pair<uint64_t,uint64_t> >());

				#if defined(PV_DEBUG)
				err.put(' ');
				for ( uint64_t i = 0; C[i].first; ++i )
					logstr << "(" << static_cast<char>(C[i].second) << "," << C[i].first << ")";
				err << std::endl;
				#endif

				if ( C[0].first && C[0].second != 'D' )
					CO.push(COo,C[0].second);
			}

			std::reverse(CO.begin(),CO.begin()+COo);

			if ( producefull || COo >= minlen )
			{
				out << '>' << (aread+1) << '/' << wellcounter++ << '/' << first << '_' << first + COo << " A=[" << first << "," << last << "]" << "\n";
				char const * zp = CO.begin();
				char const * ze = zp + COo;
				uint64_t const linewidth = 80;
				while ( zp != ze )
				{
					uint64_t const r = ze-zp;
					uint64_t const tp = std::min(r,linewidth);
					out.write(zp,tp);
					out.put('\n');
					zp += tp;
				}
			}
		}
	}
};
#endif
