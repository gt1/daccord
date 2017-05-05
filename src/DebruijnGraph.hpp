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
#if ! defined(DEBRUIJNGRAPH_HPP)
#define DEBRUIJNGRAPH_HPP

#include <DebruijnGraphBase.hpp>
#include <Node.hpp>
#include <Links.hpp>
#include <OffsetLikely.hpp>
#include <DebruijnGraphInterface.hpp>

struct NodeComparator
{
	 bool operator()(Node const & A, Node const & B) const
	 {
		return A.v < B.v;
	 }
};

struct Stretch
{
	uint64_t first;
	uint64_t ext;
	uint64_t last;
	uint64_t len;
	#if defined(SWEIGHT)
	uint64_t weight;
	#endif
	#if defined(AWEIGHT)
	uint64_t aweight;
	#endif
	uint64_t stretchO;
	uint64_t feasposO;
	uint64_t feasposL;
	uint64_t cfeasposO;
	uint64_t cfeasposL;

	bool isLoop() const
	{
		return first == last;
	}

	Stretch()
	{

	}

	Stretch(
		uint64_t const rfirst,
		uint64_t const rext,
		uint64_t const rlast,
		uint64_t const rlen,
		#if defined(SWEIGHT)
		uint64_t const rweight,
		#endif
		#if defined(AWEIGHT)
		uint64_t const raweight,
		#endif
		uint64_t const rstretchO
	) : first(rfirst), ext(rext), last(rlast), len(rlen),
	    #if defined(SWEIGHT)
	    weight(rweight),
	    #endif
	    #if defined(AWEIGHT)
	    aweight(raweight),
	    #endif
	    stretchO(rstretchO), feasposL(0), cfeasposL(0)
	{

	}

	#if 0
	uint32_t hash() const
	{
		uint64_t A[2] = { first, ext };
		return libmaus2::hashing::EvaHash::hash2(reinterpret_cast<uint32_t const *>(&A[0]),4);
	}
	#endif

	bool operator<(Stretch const & O) const
	{
		if ( first != O.first )
			return first < O.first;
		else if ( ext != O.ext )
			return ext < O.ext;
		else if ( len != O.len )
			return len > O.len;
		else // if ( last != O.last )
			return last < O.last;
		#if defined(SWEIGHT)
		else if ( weight != O.weight )
			return weight < O.weight;
		#endif
		#if defined(AWEIGHT)
		else
			return aweight < O.aweight;
		#endif
	}

	bool operator==(Stretch const & O) const
	{
		return
			first == O.first
			&&
			ext == O.ext
			&&
			last == O.last
			&&
			len == O.len
			#if defined(SWEIGHT)
			&&
			weight == O.weight
			#endif
			#if defined(AWEIGHT)
			&&
			aweight == O.aweight
			#endif
			;
	}
};

std::ostream & operator<<(std::ostream & out, Stretch const & S)
{
	out << "Stretch("
		<< "first=" << S.first
		<< ",ext=" << S.ext
		<< ",last=" << S.last
		<< ",len=" << S.len
		#if defined(SWEIGHT)
		<< ",weight=" << S.weight
		#endif
		<< ",stretchO=" << S.stretchO
		<< ",feasposO=" << S.feasposO
		<< ",feasposL=" << S.feasposL
		<< ")";
	return out;
}

struct Path
{
	uint64_t len;
	uint64_t off;
	// position of last k-mer on path
	uint64_t pos;
	double weight;
	uint64_t baselen;
	#if defined(AWEIGHT)
	uint64_t aweight;
	#endif

	Path() : len(0), off(0), pos(0), weight(0.0), baselen(0)
		#if defined(AWEIGHT)
		, aweight(0)
		#endif
	{
	}

	Path(uint64_t const rlen, uint64_t const roff, uint64_t rpos, double const rweight, uint64_t const rbaselen
		#if defined(AWEIGHT)
		, uint64_t const raweight
		#endif
	)
	: len(rlen), off(roff), pos(rpos), weight(rweight), baselen(rbaselen)
	  #if defined(AWEIGHT)
	  , aweight(raweight)
	  #endif
	{

	}

	bool operator<(Path const & P) const
	{
		return weight > P.weight;
	}
};

struct ReversePath
{
	typedef uint32_t offset_type;
	typedef uint32_t kmer_type;
	typedef double weight_type;
	typedef uint16_t pos_type;
	typedef uint16_t len_type;

	static weight_type getDefaultWeight()
	{
		return 0.0;
	}

	// offset in link list
	offset_type linkoff;
	// first kmer on path
	kmer_type front;

	// path weight
	weight_type weight;
	#if defined(AWEIGHT)
	// freq weight
	uint64_t aweight;
	#endif

	// position
	pos_type pos;

	// length in number of stretches
	len_type len;
	// length in number of bases
	len_type baselen;

	ReversePath() : linkoff(0), front(0), weight(getDefaultWeight()),
	  #if defined(AWEIGHT)
	  aweight(0),
	  #endif
	  pos(0), len(0), baselen(0)
	{

	}

	ReversePath(uint64_t const rlen, uint64_t const rlinkoff, int64_t const rpos, weight_type const rweight, uint64_t rfront, uint64_t rbaselen
		#if defined(AWEIGHT)
		, uint64_t const raweight
		#endif
	)
	: linkoff(rlinkoff), front(rfront), weight(rweight),
	  #if defined(AWEIGHT)
	  aweight(raweight),
	  #endif
	  pos(rpos), len(rlen), baselen(rbaselen)
	{

	}

	bool operator<(ReversePath const & P) const
	{
		if ( front != P.front )
			return front < P.front;
		else
			return baselen < P.baselen;
	}
};

struct PathOffComparator
{
	bool operator()(Path const & A, Path const & B) const
	{
		return A.off < B.off;
	}
};

struct ReversePathFrontComparator
{
	bool operator()(ReversePath const & A, ReversePath const & B) const
	{
		return A.front < B.front;
	}
};

struct ReversePathBaseLenComparator
{
	bool operator()(ReversePath const & A, ReversePath const & B) const
	{
		return A.baselen < B.baselen;
	}
};

struct ReversePathWeightHeapComparator
{
	bool operator()(ReversePath const & A, ReversePath const & B) const
	{
		return A.weight < B.weight;
	}
};

struct PathWeightHeapComparator
{
	bool operator()(Path const & A, Path const & B) const
	{
		return A.weight < B.weight;
	}
};

struct ReversePathWeightQueueHeapComparator
{
	bool operator()(ReversePath const & A, ReversePath const & B) const
	{
		return A.weight > B.weight;
	}
};

std::ostream & operator<<(std::ostream & out, Path const & P)
{
	return out << "Path("
		<< "len=" << P.len
		<< ",off=" << P.off
		<< ",pos=" << P.pos
		<< ",weight=" << P.weight
		<< ",baselen=" << P.baselen
		#if defined(AWEIGHT)
		<< ",aweight=" << P.aweight
		#endif
		<< ")";
}

std::ostream & operator<<(std::ostream & out, ReversePath const & P)
{
	return out << "ReversePath("
		<< "len=" << P.len
		<< ",linkoff=" << P.linkoff
		<< ",pos=" << P.pos
		<< ",weight=" << P.weight
		<< ",front=" << P.front
		<< ",baselen=" << P.baselen
		#if defined(AWEIGHT)
		<< ",aweight=" << P.aweight
		#endif
		<< ")";
}

struct EdgeActivationElement
{
	uint64_t freq;
	uint64_t nodeid;
	uint64_t edgeid;

	EdgeActivationElement() {}
	EdgeActivationElement(uint64_t const rfreq, uint64_t const rnodeid, uint64_t const redgeid)
	: freq(rfreq), nodeid(rnodeid), edgeid(redgeid)
	{}

	bool operator<(EdgeActivationElement const & E) const
	{
		if ( E.freq != freq )
			return freq > E.freq;
		else if ( nodeid != E.nodeid )
			return nodeid < E.nodeid;
		else
			return edgeid < E.edgeid;
	}
};

std::ostream & operator<<(std::ostream & out, EdgeActivationElement const & E)
{
	return out << "EdgeActivationElement("
		<< "freq=" << E.freq
		<< ",nodeid=" << E.nodeid
		<< ",edgeid=" << E.edgeid
		<< ")";
}

struct ConsensusCandidate
{
	uint64_t o;
	uint64_t l;
	double weight;
	double error;

	ConsensusCandidate()
	{

	}

	ConsensusCandidate(uint64_t const ro, uint64_t const rl, double const rweight, double const rerror)
	: o(ro), l(rl), weight(rweight), error(rerror)
	{

	}

	#if 0
	bool operator<(ConsensusCandidate const & C) const
	{
		return weight > C.weight;
	}
	#endif
};

struct ConsensusCandidateDumpHeapComparator
{
	bool operator()(ConsensusCandidate const & A, ConsensusCandidate const & B) const
	{
		return A.weight < B.weight;
	}
};

struct ConsensusCandidateHeapComparator
{
	bool operator()(ConsensusCandidate const & A, ConsensusCandidate const & B) const
	{
		return A.weight > B.weight;
	}
};

struct ConsensusCandidateErrorComparator
{
	bool operator()(ConsensusCandidate const & A, ConsensusCandidate const & B) const
	{
		return A.error < B.error;
	}
};

std::ostream & operator<<(std::ostream & out, ConsensusCandidate const & E)
{
	return out << "ConsensusCandidate("
		<< "o=" << E.o
		<< ",l=" << E.l
		<< ",weight=" << E.weight
		<< ",error=" << E.error
		<< ")";
}

struct ScoreInterval
{
	uint64_t left;
	uint64_t right;
	uint64_t current;
	double weight;
	Path P;

	ScoreInterval() {}
	ScoreInterval(uint64_t const rleft, uint64_t const rright, uint64_t const rcurrent, double const rweight, Path const & rP)
	: left(rleft), right(rright), current(rcurrent), weight(rweight), P(rP) {}

	#if 0
	bool operator<(ScoreInterval const & S) const
	{
		return weight < S.weight;
	}
	#endif

	bool operator<(ScoreInterval const & S) const
	{
		return weight > S.weight;
	}
};

std::ostream & operator<<(std::ostream & out, ScoreInterval const & S)
{
	return
		out << "ScoreInterval(" << S.left << "," << S.right << "," << S.current << "," << S.weight << ")";

}

struct SeqPos
{
	uint32_t seq;
	uint32_t pos;

	SeqPos() {}
	SeqPos(uint32_t const rseq, uint32_t const rpos) : seq(rseq), pos(rpos) {}

	bool operator<(SeqPos const & O) const
	{
		if ( pos != O.pos )
			return pos < O.pos;
		else
			return seq < O.seq;
	}

	bool operator==(SeqPos const & O) const
	{
		return
			pos == O.pos &&
			seq == O.seq;
	}
};

struct PosFreq
{
	uint32_t pos;
	uint32_t freq;

	PosFreq() {}
	PosFreq(uint32_t const rpos, uint32_t const rfreq) : pos(rpos), freq(rfreq) {}
};

std::ostream & operator<<(std::ostream & out, PosFreq const & PF)
{
	return out << "PosFreq(pos=" << PF.pos << ",freq=" << PF.freq << ")";
}

struct LevelAddElement
{
	uint64_t from;
	uint64_t to;
	uint64_t v;
	uint64_t off;

	LevelAddElement() {}
	LevelAddElement(
		uint64_t const rfrom,
		uint64_t const rto,
		uint64_t const rv,
		uint64_t const roff
	) : from(rfrom), to(rto), v(rv), off(roff)
	{

	}
};

struct NodeAddElement
{
	uint64_t v;
	uint64_t pos;

	NodeAddElement() {}
	NodeAddElement(uint64_t const rv, uint64_t const rpos)
	: v(rv), pos(rpos) {}

	bool operator<(NodeAddElement const & O) const
	{
		if ( v != O.v )
			return v < O.v;
		else
			return pos < O.pos;
	}
};


struct StretchesFirstComparator
{
	bool operator()(Stretch const & A, Stretch const & B) const
	{
		return A.first < B.first;
	}
};

struct StretchesLastComparator
{
	Stretch const * P;
	uint64_t stretcho;
	uint64_t reflast;

	StretchesLastComparator() : P(0), stretcho(0) {}
	StretchesLastComparator(Stretch const * rP, uint64_t const rstretcho) : P(rP), stretcho(rstretcho) {}

	uint64_t getLast(uint64_t const i) const
	{
		if ( i < stretcho )
			return P[i].last;
		else
			return reflast;
	}

	bool operator()(uint64_t const i, uint64_t const j) const
	{
		return getLast(i) < getLast(j);
	}
};

struct PathReversePathPair
{
	Path P;
	ReversePath RP;
	double weight;

	PathReversePathPair() : P(), RP(), weight(0.0) {}
	PathReversePathPair(
		Path const & rP,
		ReversePath const & rRP,
		double const rweight
	) : P(rP), RP(rRP), weight(rweight) {}

	bool operator<(PathReversePathPair const & O) const
	{
		return weight < O.weight;
	}
};

struct PathReversePathPairDescComp
{
	bool operator()(PathReversePathPair const & A, PathReversePathPair const & B) const
	{
		return A.weight > B.weight;
	}
};

#include <libmaus2/bitio/BitVector.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/rmq/QuickDynamicRMQ.hpp>
#include <libmaus2/wavelet/WaveletTree.hpp>
#include <libmaus2/lcs/AlignmentOneAgainstManyInterface.hpp>
#include <libmaus2/lcs/AlignmentOneAgainstManyFactory.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/util/PrefixSums.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>

template<unsigned int _kmersize>
struct DebruijnGraph :
	public DebruijnGraphInterface,
	public DebruijnGraphBase
{
	// k-mer size
	static unsigned int const kmersize = _kmersize;
	//
	typedef DebruijnGraph<kmersize> this_type;
	typedef typename libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef typename libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	// compare by kmer value
	struct KmerComp
	{
		bool operator()(uint64_t const l, uint64_t const r) const
		{
			return (l >> getKmerShift()) < (r >> getKmerShift());
		}
	};

	static unsigned int const kbits = 2*kmersize;
	static unsigned int const kbucketsbits = 8;
	static uint64_t const kbuckets = 1ull << kbucketsbits;
	static uint64_t const kbucketsmask = kbuckets - 1;
	static unsigned int const kbucketrounds = (kbits + kbucketsbits - 1)/kbucketsbits;

	static unsigned int const seqbucketsbits = 8;
	static uint64_t const seqbuckets = 1ull << seqbucketsbits;
	static uint64_t const seqbucketsmask = seqbuckets-1;

	static unsigned int const posbucketsbits = 8;
	static uint64_t const posbuckets = 1ull << posbucketsbits;
	static uint64_t const posbucketsmask = posbuckets-1;

	libmaus2::autoarray::AutoArray< uint64_t > kbuckethist;
	libmaus2::autoarray::AutoArray< uint64_t > seqbuckethist;
	libmaus2::autoarray::AutoArray< uint64_t > posbuckethist;
	libmaus2::autoarray::AutoArray< uint64_t > seqlenhist;
	libmaus2::bitio::BitVector seqlenhistset;
	libmaus2::autoarray::AutoArray< uint64_t > seqlenhistlong;

	// prenodes (kmer,seq,pos)
	libmaus2::autoarray::AutoArray<uint64_t> prenodes;
	libmaus2::autoarray::AutoArray<uint64_t> prenodestmp;
	// last node for each input sequence
	libmaus2::autoarray::AutoArray<uint64_t> last;
	// length for each sequence
	libmaus2::autoarray::AutoArray<uint32_t> seqlen;
	// mask
	uint64_t const m;
	// number of prenodes
	uint64_t numprenodes;
	// number of last nodes
	uint64_t lastn;
	uint64_t seqlenn;
	// maximum number of kmers for any input string contained
	uint64_t maxk;

	// sequence,position pairs
	libmaus2::autoarray::AutoArray<SeqPos> SP;
	libmaus2::autoarray::AutoArray<SeqPos> RSP;
	// pos,freq pairs
	libmaus2::autoarray::AutoArray<PosFreq> PF;
	// pos,freq pairs
	libmaus2::autoarray::AutoArray<PosFreq> RPF;
	// node objects
	libmaus2::autoarray::AutoArray<Node> nodes;
	// number of nodes
	uint64_t numnodes;

	// non primary edge activation heap
	libmaus2::util::FiniteSizeHeap<EdgeActivationElement> EAH;

	// stretch links (k-mer words by value on stretch)
	libmaus2::autoarray::AutoArray<uint64_t> stretchLinks;
	// current add pointer for stretchLinks
	uint64_t stretchLinksO;

	std::ostream & printSize(std::ostream & out) const
	{
		uint64_t Afeasbucks = 0;
		for ( uint64_t i = 0; i < Afeasbuck.size(); ++i )
			if ( Afeasbuck[i].second )
			{
				Afeasbucks += Afeasbuck[i].second->byteSize();
			}
		Afeasbucks += Afeasbuck.byteSize();

		out
			<< "[S] k=" << kmersize << '\n'
			<< "[S] prenodes=" << prenodes.byteSize() << '\n'
			<< "[S] last=" << last.byteSize() << '\n'
			<< "[S] m=" << sizeof(m) << '\n'
			<< "[S] numprenodes=" << sizeof(numprenodes) << '\n'
			<< "[S] lastn=" << sizeof(lastn) << '\n'
			<< "[S] maxk=" << sizeof(maxk) << '\n'
			<< "[S] SP=" << SP.byteSize() << '\n'
			<< "[S] PF=" << PF.byteSize() << '\n'
			<< "[S] nodes=" << nodes.byteSize() << '\n'
			<< "[S] numnodes=" << sizeof(numnodes) << '\n'
			<< "[S] EAH=" << EAH.byteSize() << '\n'
			<< "[S] stretchLinks=" << stretchLinks.byteSize() << '\n'
			<< "[S] stretchLinksO=" << sizeof(stretchLinksO) << '\n'
			<< "[S] stretcho=" << sizeof(stretcho) << '\n'
			<< "[S] stretches=" << stretches.byteSize() << '\n'
			<< "[S] stretchBV=" << stretchBV.byteSize() << '\n'
			<< "[S] markBV=" << markBV.byteSize() << '\n'
			<< "[S] Acons=" << Acons.byteSize() << '\n'
			<< "[S] conso=" << sizeof(conso) << '\n'
			<< "[S] splitA=" << splitA.byteSize() << '\n'
			<< "[S] AP=" << AP.byteSize() << '\n'
			<< "[S] APR=" << APR.byteSize() << '\n'
			<< "[S] APRo=" << sizeof(APRo) << '\n'
			<< "[S] ACC=" << ACC.byteSize() << '\n'
			<< "[S] ACCo=" << sizeof(ACCo) << '\n'
			<< "[S] Afeaspos=" << Afeaspos.byteSize() << '\n'
			<< "[S] Afeasposo=" << sizeof(Afeasposo) << '\n'
			<< "[S] nodecache=" << nodecache.byteSize() << '\n'
			<< "[S] LS=" << LS.byteSize() << '\n'
			<< "[S] Atmpp=" << Atmpp.byteSize() << '\n'
			<< "[S] ANE=" << ANE.byteSize() << '\n'
			<< "[S] PQ=" << PQ.byteSize() << '\n'
			<< "[S] SIQ=" << SIQ.byteSize() << '\n'
			<< "[S] Astretchfeas=" << Astretchfeas.byteSize() << '\n'
			<< "[S] Afeasbuck=" << Afeasbucks << '\n'
			<< "[S] reverseStretchLinks=" << reverseStretchLinks.byteSize() << '\n'
			<< "[S] reverseStretchLinksO=" << sizeof(reverseStretchLinksO) << '\n'
			<< "[S] ARP=" << ARP.byteSize() << '\n'
			<< "[S] ARPo=" << sizeof(ARPo) << '\n'
			<< "[S] ARWT=" << ARWT.byteSize() << '\n'
			<< "[S] ARW=" << ARW.byteSize() << '\n'
			<< "[S] ARWR=" << ARWR.byteSize() << '\n'
			<< "[S] ARWR_RMQ=" << (ARWR_RMQ ? ARWR_RMQ->byteSize() : 0) << '\n'
			<< "[S] ARWWT=" << (ARWWT ? ARWWT->byteSize() : 0) << '\n'
			<< "[S] maxkmerpos=" << sizeof(maxkmerpos) << '\n'
			<< "[S] maxstretchlength=" << sizeof(maxstretchlength) << '\n'
			;

		return out;
	}

	// current add pointer for stretches
	uint64_t stretcho;
	// stretch objects
	libmaus2::autoarray::AutoArray<Stretch> stretches;
	// bit vector for loop detection while computing stretches
	libmaus2::bitio::BitVector stretchBV;
	// bit vector for marking used positions
	libmaus2::bitio::BitVector markBV;

	// consensus
	libmaus2::autoarray::AutoArray<uint8_t> Acons;
	// output pointer for consensus
	uint64_t conso;

	// stretch splitting support
	libmaus2::autoarray::AutoArray<uint64_t> splitA;
	// path elements support
	libmaus2::autoarray::AutoArray<uint64_t> AP;
	uint64_t APo;

	// reverse path elements support
	libmaus2::autoarray::AutoArray<uint64_t> APR;
	uint64_t APRo;

	// consensus candidate vector
	libmaus2::autoarray::AutoArray<ConsensusCandidate> ACC;
	// output pointer for ACC
	uint64_t ACCo;

	// feasible positions vector for k-mers
	libmaus2::autoarray::AutoArray< std::pair<uint64_t,double> > Afeaspos;
	// output pointer for Afeaspos
	uint64_t Afeasposo;

	// feasible positions vector for k-mers
	libmaus2::autoarray::AutoArray< std::pair<uint64_t,double> > Acfeaspos;
	// output pointer for Afeaspos
	uint64_t Acfeasposo;

	// node cache (pointers k-mer value -> node object id)
	typedef int32_t nodecache_value_type;
	libmaus2::autoarray::AutoArray < nodecache_value_type > nodecache;

	// support data structure for getLevelSuccessors (for filling in supposedly missing k-mers)
	libmaus2::autoarray::AutoArray<LevelAddElement> LS;
	// support data structure for getLevelSuccessors (for computing consistent positons between front and end k-mer)
	libmaus2::autoarray::AutoArray< std::pair<uint64_t,double> > Atmpp;
	// node add elements for adding supposedly missing k-mers
	libmaus2::autoarray::AutoArray< NodeAddElement > ANE;

	// path todo list for traverse
	libmaus2::util::FiniteSizeHeap<Path> PQ;
	// score interval todo list for traverse
	libmaus2::util::FiniteSizeHeap<ScoreInterval> SIQ;

	// possible start(!) positions of stretches
	struct StretchFeasObject
	{
		uint64_t p;
		double w;
		double wf;
		double wl;

		StretchFeasObject() {}
		StretchFeasObject(
			uint64_t const rp,
			double const rw,
			double const rwf,
			double const rwl
		) : p(rp), w(rw), wf(rwf), wl(rwl) {}
	};
	libmaus2::autoarray::AutoArray < StretchFeasObject > Astretchfeas;
	libmaus2::autoarray::AutoArray < StretchFeasObject > Acstretchfeas;

	// bucket sorting support for computing coherent positions
	libmaus2::autoarray::AutoArray <
		std::pair<uint64_t, libmaus2::autoarray::AutoArray < double >::shared_ptr_type >
	> Afeasbuck;

	// reverse stretch link pairs (to,from)
	libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > reverseStretchLinks;
	// output pointer for reverseStretchLinks
	uint64_t reverseStretchLinksO;

	// reverse path storage
	libmaus2::autoarray::AutoArray < ReversePath > ARP;
	// output pointer for ARP
	uint64_t ARPo;

	// reverse path average weight
	libmaus2::autoarray::AutoArray< std::pair<double,uint64_t> > ARWT;
	// reverse permutation for positions in ARWT after sorting by weight
	// essentially the weight sequence of ARWT before sorting turned into an increasing integer sequence by weight
	libmaus2::autoarray::AutoArray< uint64_t > ARW;
	// similar to ARW but heighest weight gets lowest integer (for turning the range max queries we want into range min queries)
	libmaus2::autoarray::AutoArray< uint64_t > ARWR;
	// range minimum support for ARWR / range maximum support for ARW and thus ARWT
	libmaus2::rmq::QuickDynamicRMQ < uint64_t const * >::unique_ptr_type ARWR_RMQ;
	// wavelet tree for range previous queries on ARW
	typedef libmaus2::rank::ERank222B rank_type;
	typedef libmaus2::wavelet::WaveletTree<rank_type,uint64_t> wt_type;
	typedef wt_type::unique_ptr_type wt_ptr_type;
	wt_ptr_type ARWWT;

	uint64_t maxkmerpos;
	uint64_t maxstretchlength;

	typedef libmaus2::util::FiniteSizeHeap<ReversePath,ReversePathWeightHeapComparator> rph_type;
	typedef rph_type::shared_ptr_type rph_ptr_type;
	libmaus2::autoarray::AutoArray < rph_ptr_type > ARPH;
	uint64_t ARPHo;

	libmaus2::autoarray::AutoArray < libmaus2::util::FiniteSizeHeap<Path,PathWeightHeapComparator>::shared_ptr_type > Apathheap;
	uint64_t Apathheapo;

	libmaus2::util::FiniteSizeHeap<ReversePath,ReversePathWeightQueueHeapComparator> RPST;

	libmaus2::util::FiniteSizeHeap<ConsensusCandidate,ConsensusCandidateDumpHeapComparator> CDH;
	libmaus2::util::FiniteSizeHeap<ConsensusCandidate,ConsensusCandidateHeapComparator> CH;

	uint64_t maxsupto;

	libmaus2::lcs::Aligner::unique_ptr_type SNP;

	libmaus2::lcs::AlignmentOneAgainstManyInterface::unique_ptr_type algn;

	// (freq,kmer) pairs for first kmer (used by traverse)
	libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > maxFirst;
	// (freq,kmer) pairs for last kmer (used by traverse)
	libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > maxLast;

	libmaus2::autoarray::AutoArray<uint8_t> Aconstmp;


	// kmer stored in bits 32,...,63
	static unsigned int getKmerShift() { return 32; }
	// sequence stored in bits 0,...,15
	static unsigned int getSeqShift() { return 0; }
	// position stored in bits 16,...,31
	static unsigned int getPosShift() { return 16; }

	uint64_t getKmerSize() const
	{
		return kmersize;
	}

	Node const * getNode(uint64_t const v) const
	{
		nodecache_value_type const j = nodecache[v];

		if ( j < 0 )
		{
			static Node const * p = 0;
			return p;
		}
		else
		{
			#if 0
			assert ( j < numnodes );
			assert ( nodes[j].v == v );
			#endif
			return nodes.begin() + j;
		}
	}

	int64_t getNodeId(uint64_t const v) const
	{
		return nodecache[v];
	}

	Node * getNode(uint64_t const v)
	{
		nodecache_value_type const j = nodecache[v];

		if ( j < 0 )
		{
			static Node * p = 0;
			return p;
		}
		else
		{
			#if 0
			assert ( j < numnodes );
			assert ( nodes[j].v == v );
			#endif
			return nodes.begin() + j;
		}
	}

	Node const * getNodeVirtual(uint64_t const v) const
	{
		return getNode(v);
	}

	uint64_t getLevelSuccessors(
		// source node
		uint64_t const v,
		unsigned int const s,
		libmaus2::autoarray::AutoArray<LevelAddElement> & A,
		uint64_t o
	) const
	{
		uint64_t const low = (v << (2*s))&m;
		uint64_t const high = low | ::libmaus2::math::lowbits(2*s);

		Node ref;
		ref.v = low;
		Node const * nlow = ::std::lower_bound(nodes.begin(),nodes.begin()+numnodes,ref,NodeComparator());
		ref.v = high;
		Node const * nhigh = ::std::upper_bound(nodes.begin(),nodes.begin()+numnodes,ref,NodeComparator());

		// iterate over target nodes
		for ( Node const * np = nlow; np != nhigh; ++np )
		{
			// target node kmer word
			uint64_t const nv = np->v;

			// check words which should exist in between
			for ( unsigned int i = 1; i < s; ++i )
			{
				uint64_t const vhigh = (v << (2*i)) & m;
				uint64_t const vlow = nv >> ((s-i)*2);
				uint64_t const cv = vlow | vhigh;

				// if node does not exist
				if ( ! getNode(cv) )
				{
					A.push(o,LevelAddElement(v,nv,cv,i));
				}
				// std::cerr << "i=" << i << " " << decode(v) << " " << decode(cv) << " " << decode(nv) << " " << getNode(cv) << std::endl;
			}
		}

		return o;
	}

	uint64_t getLevelSuccessors(
		unsigned int const s,
		libmaus2::autoarray::AutoArray<LevelAddElement> & A
	) const
	{
		uint64_t o = 0;
		for ( uint64_t i = 0; i < numnodes; ++i )
			o = getLevelSuccessors(nodes[i].v,s,A,o);
		return o;
	}

	void getLevelSuccessors(unsigned int const s)
	{
		uint64_t const o = getLevelSuccessors(s,LS);

		uint64_t aneo = 0;
		for ( uint64_t i = 0; i < o; ++i )
		{
			LevelAddElement const & L = LS[i];

			// std::cerr << decode(L.from) << " " << decode(L.v) << " " << decode(L.to) << " " << getNode(L.v) << " " << L.off << std::endl;

			assert ( getNode(L.from) );
			assert ( getNode(L.to) );

			Node const & from = *getNode(L.from);
			Node const & to = *getNode(L.to);

			// intersect feasible positions of from and to node
			uint64_t tmpo = 0;
			for ( uint64_t j = 0; j < from.numfeaspos; ++j )
				Atmpp.push(tmpo,
					std::pair<uint64_t,double>(
						Afeaspos[from.feaspos + j].first + s,
						Afeaspos[from.feaspos + j].second
					)
				);
			for ( uint64_t j = 0; j < to.numfeaspos; ++j )
				Atmpp.push(tmpo,Afeaspos[to.feaspos + j]);

			std::sort(Atmpp.begin(),Atmpp.begin()+tmpo);

			// check for feasible common position
			uint64_t l = 0;
			uint64_t mp = 0;
			double mweight = std::numeric_limits<double>::min();
			while ( l < tmpo )
			{
				uint64_t h = l+1;
				while ( h < tmpo && Atmpp[l].first == Atmpp[h].first )
					++h;

				assert ( h-l <= 2 );

				if ( h-l > 1 )
				{
					// position on to node
					uint64_t const pp = Atmpp[l].first;

					// sanity check
					assert ( pp >= s );

					uint64_t const p = pp - s + L.off;
					double const weight = Atmpp[l].second + Atmpp[h-1].second;

					if ( weight > mweight )
					{
						mweight = weight;
						mp = p;
					}
				}

				l = h;
			}

			if ( mweight != std::numeric_limits<double>::min() )
			{
				// std::cerr << "here " << mp << " " << mweight << std::endl;
				ANE.push(aneo,NodeAddElement(L.v,mp));
			}
		}

		std::sort(ANE.begin(),ANE.begin() + aneo);

		for ( uint64_t i = 0; i < aneo; ++i )
		{
			int64_t seqid = -1;
			for ( uint64_t j = 0; j < seqlenn && seqid < 0; ++j )
				if ( ANE[i].pos+kmersize <= seqlen[j] )
					seqid = j;

			if ( seqid != -1 )
				prenodes.push(numprenodes,combine(ANE[i].v,seqid,ANE[i].pos));
		}

		// sort nodes by kmer
		std::sort(prenodes.begin(),prenodes.begin()+numprenodes);
	}

	void clearNodeCache()
	{
		for ( uint64_t i = 0; i < numnodes; ++i )
			nodecache [ nodes[i].v ] = -1;

		#if 0
		for ( uint64_t i = 0; i < nodecache.size(); ++i )
			assert ( nodecache[i] == -1 );
		#endif
	}

	void setupNodeCache()
	{
		for ( uint64_t i = 0; i < numnodes; ++i )
			nodecache [ nodes[i].v ] = i;
	}

	// filter out all kmers occuring less than f times
	void filterFreq(uint64_t const f)
	{
		clearNodeCache();

		uint64_t o = 0;
		for ( uint64_t i = 0; i < numnodes; ++i )
			if ( nodes[i].freq >= f )
				nodes[o++] = nodes[i];

		numnodes = o;

		setupNodeCache();

		setupAddHeap();

		check();
	}

	// print bits in a word (MSB to LSB)
	static std::string printBits(uint64_t v)
	{
		std::ostringstream ostr;
		for ( int i = 63; i >= 0; --i )
			ostr << ((v >> i) & 1);
		return ostr.str();
	}

	// combine kmer,seq,pos into a single 64 bit word
	static uint64_t combine(uint64_t const v, uint64_t const seq, uint64_t const p) // const
	{
		return
			(v << getKmerShift())
			|
			(seq<<getSeqShift())
			|
			(p  << getPosShift());
	}

	// get kmer from combined word
	static uint64_t kmerMask(uint64_t const v)
	{
		return v >> getKmerShift();
	}

	// get sequence id from combined word
	static uint64_t seqMask(uint64_t const v)
	{
		return (v >> getSeqShift()) & 0xFFFFull;
	}

	// get position from combined mask
	static uint64_t posMask(uint64_t const v)
	{
		return (v >> getPosShift()) & 0xFFFFull;
	}

	uint64_t kmer(uint64_t const nodeid) const
	{
		Node const & node = nodes[nodeid];
		return node.v;
	}

	uint64_t seq(uint64_t const nodeid, uint64_t const instid) const
	{
		Node const & node = nodes[nodeid];
		return SP[node.spo + instid].seq;
	}

	uint64_t pos(uint64_t const nodeid, uint64_t const instid) const
	{
		Node const & node = nodes[nodeid];
		return SP[node.spo + instid].pos;
	}

	// get kmer with maximum frequency for position p
	uint64_t maxForPos(uint64_t const p) const
	{
		uint64_t maxv = 0;
		uint64_t maxc = 0;

		for ( uint64_t i = 0; i < numnodes; ++i )
		{
			Node const & node = nodes[i];

			uint64_t c = 0;
			for ( uint64_t j = 0; j < node.freq; ++j )
				if ( SP[node.spo+j].pos == p )
					++c;
			if ( c > maxc )
			{
				maxc = c;
				maxv = node.v;
			}
		}

		return maxv;
	}

	// get kmer with maximum frequency for position p
	uint64_t maxForPosList(uint64_t const p, libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > & PL) const
	{
		uint64_t PLo = 0;

		for ( uint64_t i = 0; i < numnodes; ++i )
		{
			Node const & node = nodes[i];

			uint64_t c = 0;
			for ( uint64_t j = 0; j < node.freq; ++j )
				if ( SP[node.spo+j].pos == p )
					++c;

			if ( c )
				PL.push(PLo,std::pair<uint64_t,uint64_t>(c,node.v));
		}

		std::sort(
			PL.begin(),
			PL.begin() + PLo,
			std::greater< std::pair<uint64_t,uint64_t> >()
		);

		return PLo;
	}

	void printLast(std::ostream & out) const
	{
		uint64_t l = 0;
		while ( l < lastn )
		{
			// frequency
			uint64_t c = 1;
			// interval upper bound
			uint64_t h = l+1;
			// count
			while ( h < lastn && (last[h]>>getKmerShift()) == (last[l]>>getKmerShift()) )
			{
				++c;
				++h;
			}

			out << decode( last[l] >> getKmerShift() ) << " " << h-l << std::endl;

			l = h;
		}

	}

	uint64_t maxLastWord() const
	{
		uint64_t maxv = 0;
		uint64_t maxc = 0;

		uint64_t l = 0;
		while ( l < lastn )
		{
			// frequency
			uint64_t c = 1;
			// interval upper bound
			uint64_t h = l+1;
			// count
			while ( h < lastn && (last[h]>>getKmerShift()) == (last[l]>>getKmerShift()) )
			{
				++c;
				++h;
			}

			if ( c > maxc )
			{
				maxv = (last[l]>>getKmerShift());
				maxc = c;
			}

			l = h;
		}

		return maxv;
	}

	uint64_t maxLastList(libmaus2::autoarray::AutoArray< std::pair<uint64_t,uint64_t> > & PL) const
	{
		uint64_t PLo = 0;

		uint64_t l = 0;
		while ( l < lastn )
		{
			// frequency
			uint64_t c = 1;
			// interval upper bound
			uint64_t h = l+1;
			// count
			while ( h < lastn && (last[h]>>getKmerShift()) == (last[l]>>getKmerShift()) )
			{
				++c;
				++h;
			}

			uint64_t const v = (last[l]>>getKmerShift());
			PL.push(PLo,std::pair<uint64_t,uint64_t>(c,v));

			l = h;
		}

		std::sort(
			PL.begin(),
			PL.begin() + PLo,
			std::greater< std::pair<uint64_t,uint64_t> >()
		);

		return PLo;
	}

	// add a single character at the end of the consensus
	void consPush(uint8_t const c)
	{
		Acons.push(conso,c);
	}

	// add a string add the end of the consensus
	void consPush(std::string const & s)
	{
		for ( uint64_t i = 0; i < s.size(); ++i )
			consPush(s[i]);
	}

	template<typename iterator>
	void consPush(iterator it, uint64_t const n)
	{
		Acons.ensureSize(conso + n);
		std::copy(it,it+n,Acons.begin()+conso);
		conso += n;
	}

	// get the currenct consensus as a string
	std::string getConsensus() const
	{
		return std::string(Acons.begin(),Acons.begin()+conso);
	}

	// decode a kmer (unpacked) to A,C,G,T
	std::string decode(uint64_t const v) const
	{
		return decode(v,kmersize);
	}

	// decode a kmer (unpacked) to A,C,G,T
	static std::string decode(uint64_t const v, unsigned int const k)
	{
		std::ostringstream ostr;

		unsigned int shift = 2*(k-1);
		for ( unsigned int i = 0; i < k; ++i, shift -= 2 )
			ostr.put(libmaus2::fastx::remapChar((v >> shift)&3));

		return ostr.str();
	}

	void consPushWord(uint64_t const w)
	{
		unsigned int shift = 2*(kmersize-1);
		for ( unsigned int i = 0; i < kmersize; ++i, shift -= 2 )
			consPush(libmaus2::fastx::remapChar((w >> shift)&3));
	}

	void consPushWord(uint64_t const w, libmaus2::autoarray::AutoArray<uint8_t> & cons, uint64_t & o) const
	{
		unsigned int shift = 2*(kmersize-1);
		for ( unsigned int i = 0; i < kmersize; ++i, shift -= 2 )
			cons.push(o,libmaus2::fastx::remapChar((w >> shift)&3));
	}

	// print contents
	std::ostream & print(std::ostream & out) const
	{
		for ( uint64_t j = 0; j < numnodes; ++j )
		{
			out << decode(nodes[j].v) << "\t" << nodes[j].freq << "\t";
			for ( uint64_t i = 0; i < nodes[j].freq; ++i )
				out << "("
					<< SP[nodes[j].spo+i].seq << ","
					<< SP[nodes[j].spo+i].pos << ")";
			out << "\n";

		}

		return out;
	}

	// print positions where needle appears in haystack
	static std::ostream & printPositions(
		std::ostream & out,
		std::string const & needle,
		std::string const & haystack
	)
	{
		std::string::size_type p = haystack.find(needle,0);

		while ( p != std::string::npos )
		{
			out << p << ";";
			p = haystack.find(needle,p+1);
		}

		return out;
	}

	std::ostream & print(
		std::ostream & out,
		// likelidhood of having a given number of instances of a kmer based on the number of instances seen
		std::vector< std::vector<double> > const & VINST,
		// likelihood vectors for positions
		OffsetLikely const & offsetLikely,
		// reference text for window
		std::string const & reftext
	) const
	{
		std::map < std::string, std::vector< std::pair<double,uint64_t> > > wordPos;

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node const & node = nodes[z];
			out << decode(node.v) << "\t" << node.freq;

			// print positions where k-mer actually appears in reference (may be empty)
			out << "\trpos ";
			printPositions(out,decode(node.v),reftext);

			// iterate over number of copies
			for ( uint64_t i = 0; i < VINST.size(); ++i )
				if (
					// if vector contains count of kmers
					node.freq < VINST[i].size()
					&&
					// and probability is at least 1e-5
					VINST[i][node.freq] >= 1e-5
				)
					out << "\t" << i << ":" << VINST[i][node.freq];

			// print sequence ids and positions of read kmers
			out << "\t";
			for ( uint64_t i = 0; i < node.freq; ++i )
				out << "(" << SP[node.spo+i].seq << "," << SP[node.spo+i].pos << ")";

			// end of line
			out << std::endl;

			// get position histogram
			std::map < uint64_t, uint64_t > posmap;
			for ( uint64_t i = 0; i < node.freq; ++i )
				posmap [ pos(z,i) ]++;

			// turn posmap into vector
			std::vector < double > posvec;
			for ( std::map < uint64_t, uint64_t >::const_iterator ita = posmap.begin();
				ita != posmap.end(); ++ita )
			{
				uint64_t const idx = ita->first;
				while ( ! (idx < posvec.size() ) )
					posvec.push_back(0);
				posvec[idx] = ita->second;
			}

			// iterate over positions in offsetLikely
			for ( uint64_t i = 0; i < offsetLikely.size(); ++i )
			{
				// compute dot product
				double const v = offsetLikely[i].dotproduct(posvec);

				// if dot product is not zero
				if ( v )
				{
					// position 0
					if ( i == 0 )
					{
						// is this the end of offsetLikely?
						if ( i+1 == offsetLikely.size() )
							// consider this a maximum
							wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						// other mark if this position is more likely than the next one
						else if ( v > offsetLikely[i+1].dotproduct(posvec) )
							wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
					}
					else
					{
						// is this the end of offsetLikely?
						if ( i+1 == offsetLikely.size() )
						{
							// mark this position if it is more likely than the last one
							if ( v > offsetLikely[i-1].dotproduct(posvec) )
								wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						}
						else
						{
							// mark this position if it is larger than the ones for the previous and next position
							if (
								v > offsetLikely[i-1].dotproduct(posvec)
								&&
								v > offsetLikely[i+1].dotproduct(posvec)
							)
								wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						}
					}

					// out << "dot product i=" << i << " = " <<  offsetLikely[i].dotproduct(posvec) << std::endl;
				}
			}
		}

		// sort by likelihood for each kmer
		for ( std::map < std::string, std::vector< std::pair<double,uint64_t> > >::iterator ita = wordPos.begin(); ita != wordPos.end(); ++ita )
		{
			// get vector
			std::vector< std::pair<double,uint64_t> > & V = ita->second;
			// sort by likelihood
			std::sort(V.begin(),V.end(),std::greater< std::pair<double,uint64_t> >());
		}

		// number of reference k-mers
		uint64_t const numrefk = (reftext.size() >= kmersize) ? (reftext.size()-kmersize+1) : 0;
		// iterate over reference k-mers
		for ( uint64_t i = 0; i < numrefk; ++i )
		{
			// get window/k-mer
			std::string const sub = reftext.substr(i,kmersize);

			out << "ref " << i << " " << sub << " ";

			// if there are no occurences
			if ( wordPos.find(sub) == wordPos.end() )
				out << "NO";
			else
			{
				// get vector for k-mer
				std::vector< std::pair<double,uint64_t> > const V = wordPos.find(sub)->second;
				// output by decreasing likelihood of position
				for ( uint64_t j = 0; j < V.size(); ++j )
					out << V[j].second << "," << V[j].first << ";";
			}

			out << std::endl;
		}

		return out;
	}

	// same as above, only without support for likelidhood of having a given number of instances of a kmer based on the number of instances seen
	std::ostream & print(
		std::ostream & out,
		// dot products for computing likelihood of offset for kmer by list of positions where it occurs
		OffsetLikely const & offsetLikely,
		// reference text window
		std::string const & reftext
	) const
	{
		std::map < std::string, std::vector< std::pair<double,uint64_t> > > wordPos;

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node const & node = nodes[z];

			out << decode(node.v) << "\t" << node.freq;

			out << "\trpos ";
			printPositions(out,decode(node.v),reftext);

			std::map < uint64_t, uint64_t > posmap;

			out << "\t";
			for ( uint64_t i = 0; i < node.freq; ++i )
				out << "(" << seq(z,i) << "," << pos(z,i) << ")";

			out << std::endl;

			for ( uint64_t i = 0; i < node.freq; ++i )
				posmap [ pos(z,i) ]++;

			std::vector < double > posvec;
			for ( std::map < uint64_t, uint64_t >::const_iterator ita = posmap.begin();
				ita != posmap.end(); ++ita )
			{
				uint64_t const idx = ita->first;
				while ( ! (idx < posvec.size() ) )
					posvec.push_back(0);
				posvec[idx] = ita->second;
			}

			for ( uint64_t i = 0; i < offsetLikely.size(); ++i )
			{
				double const v = offsetLikely[i].dotproduct(posvec);

				if ( v )
				{
					if ( i == 0 )
					{
						if ( i+1 == offsetLikely.size() )
							wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						else if ( v > offsetLikely[i+1].dotproduct(posvec) )
							wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
					}
					else
					{
						if ( i+1 == offsetLikely.size() )
						{
							if ( v > offsetLikely[i-1].dotproduct(posvec) )
								wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						}
						else
						{
							if (
								v > offsetLikely[i-1].dotproduct(posvec)
								&&
								v > offsetLikely[i+1].dotproduct(posvec)
							)
								wordPos[decode(node.v)].push_back(std::pair<double,uint64_t>(v,i));
						}
					}

					// out << "dot product i=" << i << " = " <<  offsetLikely[i].dotproduct(posvec) << std::endl;
				}
			}
		}

		for ( std::map < std::string, std::vector< std::pair<double,uint64_t> > >::iterator ita = wordPos.begin(); ita != wordPos.end(); ++ita )
		{
			std::vector< std::pair<double,uint64_t> > & V = ita->second;
			std::sort(V.begin(),V.end(),std::greater< std::pair<double,uint64_t> >());
		}

		uint64_t const numrefk = (reftext.size() >= kmersize) ? (reftext.size()-kmersize+1) : 0;
		for ( uint64_t i = 0; i < numrefk; ++i )
		{
			std::string const sub = reftext.substr(i,kmersize);

			out << "ref " << i << " " << sub << " ";

			if ( wordPos.find(sub) == wordPos.end() )
				out << "NO";
			else
			{
				std::vector< std::pair<double,uint64_t> > const V = wordPos.find(sub)->second;
				for ( uint64_t j = 0; j < V.size(); ++j )
					out << V[j].second << "," << V[j].first << ";";
			}

			out << std::endl;
		}

		return out;
	}

	void check()
	{
		for ( uint64_t z = 1; z < numnodes; ++z )
			assert ( nodes[z-1].v < nodes[z].v );

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Links L;
			Node const & node = nodes[z];
			getSuccessors(node.v,L);
			assert ( L.size() == node.numsucc );
			assert ( node.numsuccactive <= node.numsucc );

			getActiveSuccessors(node.v,L);
			// assert ( L.size() == node.numsuccactive );
		}
	}

	void setupAddHeap()
	{
		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node & node = nodes[z];
			node.numsucc = 0;
			node.numsuccactive = 0;
		}

		Links L;

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node & node = nodes[z];
			// get successors for v
			getSuccessors(node.v,L);

			// check number of successors
			if ( L.size() )
			{
				node.numsucc = L.size();
				node.numsuccactive = 1;

				while (
					node.numsuccactive < L.size()
					&&
					L.getFreq(node.numsuccactive) >= L.getFreq(0)/2
				)
					++node.numsuccactive;

				#if 0
				assert ( node.numsuccactive <= L.size() );
				assert ( node.numsuccactive <= node.numsucc );
				#endif
			}
			else
			{
				node.numsucc = 0;
				node.numsuccactive = 0;
			}
		}

		EAH.clear();

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node const & node = nodes[z];
			// get successors for v
			getSuccessors(node.v,L);

			for ( uint64_t i = node.numsuccactive; i < L.size(); ++i )
			{
				EAH.pushBump(
					EdgeActivationElement(
						L.getFreq(i),
						z /* nodeid */,
						i /* edgeid */
					)
				);
			}
		}
	}

	bool addNextFromHeap(std::ostream * /* logstr */)
	{
		if ( EAH.empty() )
		{
			return false;
		}
		else
		{
			#if 1
			uint64_t const topfreq = EAH.top().freq;

			#if 0
			if ( logstr )
				(*logstr) << "addNextFrom heap adding freq " << topfreq << std::endl;
			#endif

			while ( (!EAH.empty()) && (EAH.top().freq == topfreq) )
			{
				EdgeActivationElement EAE = EAH.pop();

				nodes[EAE.nodeid].numsuccactive += 1;

				#if 0
				assert ( nodes[EAE.nodeid].numsuccactive <= nodes[EAE.nodeid].numsucc );
				#endif
			}

			#else

			EdgeActivationElement EAE = EAH.pop();
			nodes[EAE.nodeid].numsuccactive += 1;

			#endif

			return true;
		}
	}

	void setupFeasBuckets(uint64_t const mpos)
	{
		if ( Afeasbuck.size() < mpos+1 )
		{
			uint64_t const oldsize = Afeasbuck.size();
			Afeasbuck.ensureSize(mpos+1);

			for ( uint64_t i = oldsize; i < Afeasbuck.size(); ++i )
			{
				libmaus2::autoarray::AutoArray<double>::shared_ptr_type A(new libmaus2::autoarray::AutoArray<double>(64,false));
				Afeasbuck[i] = std::pair < uint64_t, libmaus2::autoarray::AutoArray<double>::shared_ptr_type >(0,A);
			}
		}

		assert ( mpos < Afeasbuck.size() );

		markBV.ensureSize(mpos+1);
	}

	void setupNodes()
	{
		clearNodeCache();
		numnodes = 0;

		uint64_t l = 0;
		uint64_t spo = 0;
		uint64_t rspo = 0;
		uint64_t pfo = 0;
		uint64_t rpfo = 0;
		maxkmerpos = 0;

		while ( l < numprenodes )
		{
			uint64_t h = l;
			uint64_t li = l;
			uint64_t lp = posMask(prenodes[l]);
			uint64_t subfreq = 0;
			uint64_t const pfostart = pfo;
			uint64_t const rpfostart = rpfo;
			while ( h < numprenodes && kmerMask(prenodes[h]) == kmerMask(prenodes[l]) )
			{
				uint64_t const seq = seqMask(prenodes[h]);
				uint64_t const pos = posMask(prenodes[h]);

				if ( pos != lp )
				{
					uint64_t const pfreq = h-li;
					subfreq += pfreq;
					PF.push(pfo,PosFreq(lp,pfreq));

					lp = pos;
					li = h;
				}

				SP.push(spo,SeqPos(seq,pos));
				assert ( pos+kmersize <= seqlen[seq] );
				RSP.push(rspo,SeqPos(seq,seqlen[seq]-pos-kmersize));
				++h;
			}
			assert ( h > li );
			uint64_t const pfreq = h-li;
			subfreq += pfreq;
			PF.push(pfo,PosFreq(lp,pfreq));

			uint64_t const freq = h-l;
			assert ( freq == subfreq );

			std::sort(RSP.begin() + rspo-freq, RSP.begin() + rspo);
			for ( uint64_t i = 1; i < freq; ++i )
				assert ( RSP[rspo-freq+i-1].pos <= RSP[rspo-freq+i].pos );
			assert ( spo == rspo );

			uint64_t cl = rspo - freq;
			while ( cl < rspo )
			{
				uint64_t ch = cl+1;
				while ( ch < rspo && RSP[ch].pos == RSP[cl].pos )
					++ch;

				RPF.push(rpfo,PosFreq(RSP[cl].pos,ch-cl));

				cl = ch;
			}

			for ( uint64_t i = 1; i < rpfo-rpfostart; ++i )
				assert ( RPF[rpfostart + i - 1].pos < RPF[rpfostart + i].pos );

			uint64_t const mincpos = RPF[rpfostart].pos;
			uint64_t const maxcpos = RPF[rpfo-1].pos;

			uint64_t const kv = kmerMask(prenodes[l]);
			Node const node(kv,spo-freq,freq,pfostart,pfo-pfostart,rpfostart,rpfo-rpfostart,PF[pfostart].pos,PF[pfo-1].pos,mincpos,maxcpos);
			#if 0
			std::cerr << node << std::endl;
			for ( uint64_t i = 0; i < node.pfosize; ++i )
				std::cerr << PF[node.pfostart + i] << std::endl;
			#endif
			nodes.push(numnodes,node);

			maxkmerpos = std::max(maxkmerpos,std::max(node.cphigh,node.phigh));

			for ( uint64_t i = 1; i < node.pfosize; ++i )
				assert ( PF [ node.pfostart + i - 1].pos < PF [ node.pfostart + i ].pos );
			for ( uint64_t i = 1; i < node.cpfosize; ++i )
				assert ( RPF [ node.cpfostart + i - 1].pos <  RPF [ node.cpfostart + i].pos );

			assert ( PF [ node.pfostart + node.pfosize - 1 ].pos <= maxkmerpos );
			assert ( RPF [ node.cpfostart + node.cpfosize - 1 ].pos <= maxkmerpos );

			l = h;
		}

		setupFeasBuckets(maxkmerpos);

		setupNodeCache();
	}

	#define CONS_BUCK_SORT

	void setupPreNodes(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o)
	{
		numprenodes = 0;
		lastn = 0;
		seqlenn = 0;
		maxk = 0;

		// if kmer length is not zero
		if ( kmersize )
		{
			#if defined(CONS_BUCK_SORT)
			std::fill(kbuckethist.begin(),kbuckethist.end(),0ull);

			uint64_t const seqbits = o ? libmaus2::math::numbits(o-1) : 0;
			uint64_t const seqsortrounds = (seqbits + seqbucketsbits - 1)/seqbucketsbits;
			seqbuckethist.ensureSize(seqsortrounds * seqbuckets);
			std::fill(seqbuckethist.begin(),seqbuckethist.begin() + seqsortrounds * seqbuckets, 0ull);

			uint64_t seqlenhisto = 0;
			uint64_t seqlenhistlongo = 0;
			uint64_t seqactive = 0;
			#endif

			// iterator over strings
			for ( uint64_t j = 0; j < o; ++j )
			{
				// if string is sufficiently long
				if ( I[j].second >= kmersize )
				{
					// number of kmers
					uint64_t const numk = I[j].second - kmersize + 1;
					// base ptr
					uint8_t const * u = I[j].first;

					#if defined(CONS_BUCK_SORT)
					seqactive += 1;

					uint64_t vseq = j;
					uint64_t * pseqbuckethist = seqbuckethist.begin();
					for ( uint64_t i = 0; i < seqsortrounds; ++i )
					{
						uint64_t const b = vseq & seqbucketsmask;
						pseqbuckethist[b] += numk;
						vseq >>= seqbucketsbits;
						pseqbuckethist += seqbuckets;
					}

					if ( numk < seqlenhist.size() )
					{
						if ( ! seqlenhistset.get(numk) )
						{
							seqlenhistset.set(numk);
							seqlenhisto += 1;
						}
						seqlenhist[numk]++;
					}
					else
					{
						seqlenhistlong.push(seqlenhistlongo,numk);
					}
					#endif

					// values
					uint64_t v = 0;
					// process first k-1 symbols
					for ( unsigned int i = 0; i < kmersize-1; ++i )
					{
						v <<= 2;
						v |= libmaus2::fastx::mapChar(*(u++));
					}
					// process complete kmers
					for ( uint64_t i = 0; i < numk; ++i )
					{
						v <<= 2;
						v &= m;
						v |= libmaus2::fastx::mapChar(*(u++));
						prenodes.push(numprenodes,combine(v,j,i));
						// assert ( i + k <= I[j].second );

						#if defined(CONS_BUCK_SORT)
						uint64_t * pbuckhist = kbuckethist.begin();
						uint64_t vv = v;
						for ( uint64_t j = 0; j < kbucketrounds; ++j )
						{
							pbuckhist[vv & kbucketsmask]++;
							vv >>= kbucketsbits;
							pbuckhist += kbuckets;
						}
						#endif
					}
					last.push(lastn,combine(v,j,numk-1));

					maxk = std::max(maxk,numk);
				}

				seqlen.push(seqlenn,I[j].second);
			}

			#if defined(CONS_BUCK_SORT)
			uint64_t const posbits = maxk ? libmaus2::math::numbits(maxk-1) : 0;
			uint64_t const posrounds = (posbits + posbucketsbits - 1)/posbucketsbits;
			posbuckethist.ensureSize(posrounds * posbuckets);
			std::fill(posbuckethist.begin(),posbuckethist.begin() + posrounds * posbuckets,0ull);

			if ( seqlenhistlongo )
				std::sort(seqlenhistlong.begin(),seqlenhistlong.begin() + seqlenhistlongo);

			uint64_t nextp = 0;
			uint64_t ii = 0;
			for ( uint64_t i = 0; i < seqlenhisto; ++i )
			{
				ii = seqlenhistset.next1(ii);

				#if defined(CONS_BUCKET_DEBUG)
				assert ( ii < seqlenhist.size() );
				assert ( seqlenhistset.get(ii) );
				assert ( seqlenhist[ii] );
				assert ( seqactive >= seqlenhist[ii] );
				#endif

				for ( ; nextp < ii; ++nextp )
				{
					uint64_t * pbuckethist = posbuckethist.begin();
					uint64_t pv = nextp;
					for ( uint64_t z = 0; z < posrounds; ++z )
					{
						uint64_t const b = pv & posbucketsmask;
						pbuckethist[b] += seqactive;
						pv >>= posbucketsbits;
						pbuckethist += posbuckets;
					}
				}

				seqactive -= seqlenhist[ii];

				seqlenhistset.erase(ii);
				seqlenhist[ii] = 0;

				ii += 1;
			}
			#if defined(CONS_BUCKET_DEBUG)
			for ( uint64_t i = 0; i < seqlenhist.size(); ++i )
			{
				assert ( ! seqlenhist[i] );
				assert ( ! seqlenhistset.get(i) );
			}
			#endif
			uint64_t lo = 0;
			while ( lo < seqlenhistlongo )
			{
				uint64_t ho = lo+1;
				while ( ho < seqlenhistlongo && seqlenhistlong[ho] == seqlenhistlong[lo] )
					++ho;

				uint64_t const ii = seqlenhistlong[lo];
				uint64_t const freq = ho - lo;

				#if defined(CONS_BUCKET_DEBUG)
				assert ( seqactive >= freq );
				#endif

				for ( ; nextp < ii; ++nextp )
				{
					uint64_t * pbuckethist = posbuckethist.begin();
					uint64_t pv = nextp;
					for ( uint64_t z = 0; z < posrounds; ++z )
					{
						uint64_t const b = pv & posbucketsmask;
						pbuckethist[b] += seqactive;
						pv >>= posbucketsbits;
						pbuckethist += posbuckets;
					}
				}

				seqactive -= freq;
				lo = ho;
			}
			assert ( ! seqactive );

			for ( uint64_t i = 0; i < seqsortrounds; ++i )
			{
				uint64_t * const pbuckhist = seqbuckethist.begin() + i * seqbuckets;
				#if defined(CONS_BUCKET_DEBUG)
				assert (
					std::accumulate(
						pbuckhist,pbuckhist+seqbuckets,0ull
						) == numprenodes
				);
				#endif
				#if defined(CONS_BUCKET_DEBUG)
				uint64_t const s =
				#endif
					libmaus2::util::PrefixSums::prefixSums(pbuckhist,pbuckhist + seqbuckets);
				#if defined(CONS_BUCKET_DEBUG)
				assert ( s == numprenodes );
				#endif
				unsigned int const shift = i * seqbucketsbits + getSeqShift();

				prenodestmp.ensureSize(numprenodes);

				for ( uint64_t i = 0; i < numprenodes; ++i )
				{
					uint64_t const v = prenodes[i];
					uint64_t const bv = (v >> shift) & seqbucketsmask;
					prenodestmp [ pbuckhist[bv]++ ] = v;
				}

				prenodes.swap(prenodestmp);
			}

			for ( uint64_t i = 0;  i < posrounds; ++i )
			{
				uint64_t * pbuckethist = posbuckethist.begin() + i * posbuckets;
				#if defined(CONS_BUCKET_DEBUG)
				assert (
					std::accumulate(pbuckethist,pbuckethist+posbuckets,0ull
					) == numprenodes
				);
				#endif
				#if defined(CONS_BUCKET_DEBUG)
				uint64_t const s =
				#endif
					libmaus2::util::PrefixSums::prefixSums(pbuckethist,pbuckethist + posbuckets);
				#if defined(CONS_BUCKET_DEBUG)
				assert ( s == numprenodes );
				#endif
				unsigned int const shift = i * posbucketsbits + getPosShift();
				prenodestmp.ensureSize(numprenodes);

				for ( uint64_t i = 0; i < numprenodes; ++i )
				{
					uint64_t const v = prenodes[i];
					uint64_t const bv = (v >> shift) & posbucketsmask;
					prenodestmp [ pbuckethist[bv]++ ] = v;
				}

				prenodes.swap(prenodestmp);
			}

			for ( uint64_t i = 0; i < kbucketrounds; ++i )
			{
				#if defined(CONS_BUCKET_DEBUG)
				assert (
					std::accumulate(
						kbuckethist.begin() + i * kbuckets,
						kbuckethist.begin() + (i+1) * kbuckets,
						0ull
					)
					==
					numprenodes
				);
				#endif

				uint64_t * const pbuckhist = kbuckethist.begin() + i * kbuckets;
				#if defined(CONS_BUCKET_DEBUG)
				uint64_t const s =
				#endif
					libmaus2::util::PrefixSums::prefixSums(pbuckhist,pbuckhist + kbuckets);
				#if defined(CONS_BUCKET_DEBUG)
				assert ( s == numprenodes );
				#endif
				unsigned int const shift = i * kbucketsbits + getKmerShift();
				prenodestmp.ensureSize(numprenodes);

				for ( uint64_t i = 0; i < numprenodes; ++i )
				{
					uint64_t const v = prenodes[i];
					uint64_t const bv = (v >> shift) & kbucketsmask;
					prenodestmp [ pbuckhist[bv]++ ] = v;
				}

				prenodes.swap(prenodestmp);
			}

			// std::sort(prenodes.begin(),prenodes.begin()+numprenodes);

			for ( uint64_t i = 1; i < numprenodes; ++i )
				assert ( prenodes[i-1] <= prenodes[i] );
			#else
			// sort nodes by kmer
			std::sort(prenodes.begin(),prenodes.begin()+numprenodes);
			#endif

			// sort last nodes by kmer
			std::sort(last.begin(),last.begin()+lastn);
		}
	}

	// setup using list of (string,length) pairs
	void setup(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o)
	{
		stretcho = 0;
		clearNodeCache();
		numnodes = 0;
		EAH.clear();

		// number of kmers
		numprenodes = 0;
		// number of last kmers
		lastn = 0;
		seqlenn = 0;
		// maximum length of a string in number of kmers
		maxk = 0;

		setupPreNodes(I,o);

		setupNodes();

		setupAddHeap();

		#if 0
		check();
		#endif
	}

	#define CDH_SIZE 16
	#define CD_SIZE 16
	#define REVERSE_PATH_HEAP_SIZE 12
	#define PATH_HEAP_SIZE 12


	// constructor
	DebruijnGraph()
	:
	  kbuckethist(kbuckets * kbucketrounds), seqlenhist(256), seqlenhistset(seqlenhist.size()),
	  m(libmaus2::math::lowbits(2*kmersize)), numnodes(0), EAH(1024), stretchBV(0), markBV(0), conso(0),
	  nodecache(1ull << (2*kmersize),false),
	  PQ(1024), SIQ(1024), maxkmerpos(0), maxstretchlength(0), ARPHo(0), Apathheapo(0),
	  RPST(1024), CDH(CDH_SIZE), CH(CD_SIZE), maxsupto(0),
	  SNP(DebruijnGraphBase::getAligner()),
	  algn(libmaus2::lcs::AlignmentOneAgainstManyFactory::uconstruct())
	{
		std::fill(seqlenhist.begin(),seqlenhist.end(),0ull);
		std::fill(nodecache.begin(),nodecache.end(),-1);
	}

	virtual ~DebruijnGraph() {}

	// get count for kmer v
	uint64_t count(uint64_t const v) const
	{
		Node const * node = getNode(v);

		if ( node )
			return node->freq;
		else
			return 0;
	}


	void printRightLinks(std::ostream & out, uint64_t const v, Links const & L) const
	{
		out << "Node(" << decode(v) << ") links (";

		for ( uint64_t i = 0; i < L.size(); ++i )
		{
			uint64_t const ext = ((v << 2) & m) | L.getSym(i);
			out << "(" << decode(ext) << "," << L.getFreq(i) << ")";
		}

		out << ")";
	}

	// get all possible successors for v
	void getSuccessors(uint64_t const v, Links & L) const
	{
		L.reset();

		// prefix for successors
		uint64_t const masked = (v << 2)&m;

		for ( uint64_t i = 0; i < 4; ++i )
		{
			L.push(i,count(masked | i));
		}

		L.sort();
	}

	virtual void getSuccessorsVirtual(uint64_t const v, Links & L) const
	{
		return getSuccessors(v,L);
	}

	// get set of active successors of v
	void getActiveSuccessors(uint64_t const v, Links & L) const
	{
		L.reset();

		Node const * node = getNode(v);

		if ( node )
		{
			#if 0
			assert ( node );
			assert ( node->v == v );
			#endif

			getSuccessors(v,L);

			#if 0
			assert( node->numsucc == L.size() );
			assert( node->numsuccactive <= L.size() );
			#endif

			L.setSize(node->numsuccactive);
		}
		else
		{
			L.setSize(0);
		}
	}

	void getActiveSuccessorsVirtual(uint64_t const v, Links & L) const
	{
		getActiveSuccessors(v,L);
	}

	// get unique active successor
	uint64_t getUniqueActiveSuccessor(uint64_t const v) const
	{
		Links L;
		getActiveSuccessors(v,L);
		#if 0
		assert ( L.size() == 1 );
		#endif
		return ((v << 2) & m) | L.getSym(0);
	}

	// get number of active successors for word v
	uint64_t getNumActiveSuccessors(uint64_t const v) const
	{
		Links L;
		getActiveSuccessors(v,L);
		return L.size();
	}

	// return true if edge from -> to is active
	bool isEdgeActive(uint64_t const from, uint64_t const to) const
	{
		Links L;
		getActiveSuccessors(from,L);

		// prefix for successors
		uint64_t const masked = (from << 2)&m;

		for ( uint64_t i = 0; i < L.size(); ++i )
		{
			uint64_t const candto = masked | L.getSym(i);
			if ( to == candto )
				return true;
		}

		return false;
	}

	bool isEdgeActiveVirtual(uint64_t const from, uint64_t const to) const
	{
		return isEdgeActive(from,to);
	}

	// get predecessors for v
	void getPredecessors(uint64_t const v, Links & L) const
	{
		L.reset();

		// prefix for successors
		uint64_t const masked = (v >> 2)&m;
		#if 0
		assert ( k );
		#endif
		unsigned int const shift = 2*(kmersize-1);

		for ( uint64_t i = 0; i < 4; ++i )
			L.push(i,count(masked | (i << shift)));

		L.sort();
	}

	void getActivePredecessors(uint64_t const v, Links & L) const
	{
		L.reset();

		// get node for v
		Node const * node = getNode(v);

		if ( node )
		{
			#if 0
			assert ( node );
			assert ( node->v == v );
			#endif

			// prefix for predecessors
			uint64_t const masked = (v >> 2)&m;
			#if 0
			assert ( kmersize );
			#endif
			unsigned int const shift = 2*(kmersize-1);

			// get all existing predecessors
			getPredecessors(v,L);

			// iterate over all predecessors of v
			uint64_t o = 0;
			for ( uint64_t i = 0; i < L.size(); ++i )
			{
				// previous word
				uint64_t const prev = masked | (L.getSym(i) << shift);
				// if edge prev -> v is active
				if ( isEdgeActive(prev,v) )
					L.A[o++] = L.A[i];
			}

			// set number of active predecessors
			L.setSize(o);
		}
	}

	uint64_t getNumActivePredecessors(uint64_t const v) const
	{
		Links L;
		getActivePredecessors(v,L);
		return L.size();
	}

	void copyStretch(uint64_t const low, uint64_t const high)
	{
		uint64_t const first = stretchLinks[low];
		uint64_t const firstext = stretchLinks[low+1];
		uint64_t const last = stretchLinks[high-1];
		uint64_t const len = high-low;

		// uint64_t stretchLinksO = stretcho ? (stretches[stretcho-1].stretchO+stretches[stretcho-1].len) : 0;

		uint64_t stretchLinksStart = stretchLinksO;
		#if defined(SWEIGHT)
		uint64_t weight = std::numeric_limits<uint64_t>::max();
		#endif
		#if defined(AWEIGHT)
		uint64_t aweight = 0;
		#endif
		for ( uint64_t i = low; i < high; ++i )
		{
			uint64_t const link = stretchLinks[i];
			#if defined(SWEIGHT)
			weight = std::min(weight,getNode(link)->freq);
			#endif
			#if defined(AWEIGHT)
			aweight += getNode(link)->freq;
			#endif
			stretchLinks.push(stretchLinksO,link);
		}

		Stretch const stretch(first,firstext,last,len,
			#if defined(SWEIGHT)
			weight,
			#endif
			#if defined(AWEIGHT)
			aweight,
			#endif
			stretchLinksStart
		);
		stretches.push(stretcho,stretch);
	}

	void mergeStretches(Stretch A, Stretch B)
	{
		uint64_t const first = A.first;
		uint64_t const firstext = A.ext;
		uint64_t const last = B.last;
		uint64_t const len = A.len + B.len - 1;

		uint64_t stretchLinksStart = stretchLinksO;
		#if defined(SWEIGHT)
		uint64_t weight = std::min(A.weight,B.weight);
		#endif
		#if defined(AWEIGHT)
		uint64_t aweight = A.aweight + B.aweight - getNode(B.first)->freq;
		#endif

		#if defined(AWEIGHT)
		uint64_t raweight = 0;
		#endif

		for ( uint64_t i = 0; i < A.len; ++i )
		{
			uint64_t const link = stretchLinks[A.stretchO+i];
			stretchLinks.push(stretchLinksO,link);

			#if defined(AWEIGHT)
			raweight += getNode(link)->freq;
			#endif
		}
		for ( uint64_t i = 1; i < B.len; ++i )
		{
			uint64_t const link = stretchLinks[B.stretchO+i];
			stretchLinks.push(stretchLinksO,link);

			#if defined(AWEIGHT)
			raweight += getNode(link)->freq;
			#endif
		}

		#if defined(AWEIGHT)
		assert ( raweight==aweight );
		#endif

		Stretch const stretch(first,firstext,last,len,
			#if defined(SWEIGHT)
			weight,
			#endif
			#if defined(AWEIGHT)
			aweight,
			#endif
			stretchLinksStart
		);
		stretches.push(stretcho,stretch);
		maxstretchlength = std::max(maxstretchlength,len);
	}

	void mergeStretchesByIndex(uint64_t const i, uint64_t const j)
	{
		mergeStretches(stretches[i],stretches[j]);
	}


	static uint64_t encode(std::string const & ref)
	{
		uint64_t v = 0;
		for ( uint64_t i = 0; i < ref.size(); ++i )
		{
			v <<= 2;
			v |= libmaus2::fastx::mapChar(ref[i]);
		}

		return v;
	}

	void mergeStretches(uint64_t const va, uint64_t const vb)
	{
		bool merging = true;

		while ( merging )
		{
			merging = false;

			uint64_t const loopend = stretcho;

			for ( uint64_t i = 0; i < loopend; ++i )
				// if stretch is not a loop (first==last)
				if (
					! stretches[i].isLoop()
					&&
					(stretches[i].last != va)
					&&
					(stretches[i].last != vb)
				)
				{
					// search for possible continuations
					Stretch ref; ref.first = stretches[i].last;

					Stretch * a = stretches.begin();
					Stretch * e = stretches.begin()+loopend;
					std::pair<Stretch *,Stretch *> P = std::equal_range(a,e,ref,StretchesFirstComparator());

					if (
						// unique continuation
						P.second - P.first == 1
						&&
						// not merging to a loop
						(! P.first->isLoop())
						&&
						// not erased
						P.first->len
						&&
						// no loop produced
						P.first->last != stretches[i].first
					)
					{
						mergeStretchesByIndex(i, P.first - a);

						// erase first part
						stretches[i].len = 0;

						merging = true;
					}
				}

			uint64_t o = 0;
			for ( uint64_t i = 0; i < stretcho; ++i )
				if ( stretches[i].len )
					stretches[o++] = stretches[i];
			stretcho = o;
		}

		std::sort(stretches.begin(),stretches.begin()+stretcho);
	}

	void splitStretches(uint64_t const v)
	{
		uint64_t splito = 0;
		uint64_t const loopend = stretcho;

		for ( uint64_t z = 0; z < loopend; ++z )
		{
			Stretch stretch = stretches[z];

			uint64_t const start = stretch.stretchO;
			uint64_t const len = stretch.len;

			#if 0
			assert ( len >= 2 );
			#endif

			int64_t splitindex = -1;
			for ( uint64_t i = 1; i < len-1; ++i )
				if ( stretchLinks[start+i] == v )
				{
					splitindex=i;
					break;
				}
			if ( splitindex != -1 )
			{
				#if 0
				assert ( stretchLinks [ start + splitindex ] == v );
				#endif

				//uint64_t const spa = stretcho;
				copyStretch(start, start + splitindex + 1);

				//uint64_t const spb = stretcho;
				copyStretch(start + splitindex, start + len);

				splitA.push(splito,z);
			}
		}

		// input index on stretches
		uint64_t l = 0;
		// input index on splitA
		uint64_t idx = 0;
		// output index on steetches
		uint64_t o = 0;

		// while we have not reached the end of the filter list
		for ( ; idx < splito ; ++l )
		{
			#if 0
			assert ( l < loopend );
			#endif

			if ( l == splitA[idx] )
			{
				++idx;
			}
			else
			{
				stretches[o++] = stretches[l];
			}
		}

		while ( l < stretcho )
		{
			stretches[o++] = stretches[l++];
		}

		stretcho = o;
	}

	// use positions for stretch computation!!!
	void computeStretches(bool const checkpredecessors, bool const debug)
	{
		stretchBV.ensureSize(numnodes);
		stretchLinksO = 0;
		stretcho = 0;
		maxstretchlength = 0;

		for ( uint64_t z = 0; z < numnodes; ++z )
		{
			Node const & node = nodes[z];
			uint64_t const refk = node.v;

			uint64_t const numpred = getNumActivePredecessors(refk);
			uint64_t const numsucc = node.numsuccactive;

			if ( numsucc && (numpred != 1 || numsucc > 1) )
			{
				Links L;
				getActiveSuccessors(refk,L);

				// iterate over active successors
				for ( uint64_t i = 0; i < numsucc; ++i )
				{
					uint64_t stretchLinksStart = stretchLinksO;

					uint64_t const first = refk;
					uint64_t const firstext = ((refk << 2) & m) | L.getSym(i);

					// extension kmer
					uint64_t extk = firstext;

					#if defined(SWEIGHT)
					uint64_t weight = std::min(count(refk),count(firstext));
					#endif

					stretchLinks.push(stretchLinksO,refk);
					stretchBV.set(getNodeId(refk));

					stretchLinks.push(stretchLinksO,extk);
					stretchBV.set(getNodeId(extk));

					// std::cerr << "extk=" << decode(extk) << " count=" << count(extk) << std::endl;

					// length of path
					uint64_t len = 2;
					bool loop = (refk == extk);

					while (
						(!loop)
						&&
						(getNumActiveSuccessors(extk) == 1)
						&&
						(
							(!checkpredecessors)
							||
							(getNumActivePredecessors(extk) == 1)
						)
					)
					{
						// go to active successor
						extk = getUniqueActiveSuccessor(extk);
						#if defined(SWEIGHT)
						// update minimum path weight
						weight = std::min(weight,count(extk));
						#endif

						// push link
						stretchLinks.push(stretchLinksO,extk);
						// update path length
						len += 1;

						uint64_t const extid = getNodeId(extk);

						// loop?
						if ( stretchBV.get(extid) )
							loop = true;
						else
							stretchBV.set(extid);
					}

					uint64_t const last = extk;

					// erase stretchBV
					for ( uint64_t j = stretchLinksStart; j < stretchLinksStart + len; ++j )
						stretchBV.erase(getNodeId(stretchLinks[j]));

					if ( loop && (first != last) )
					{
						uint64_t j = 0;
						while ( stretchLinks[stretchLinksStart+j] != last )
							++j;

						j += 1;

						uint64_t const retract = len-j;

						len -= retract;
						stretchLinksO -= retract;

						#if defined(SWEIGHT)
						weight = std::min(count(refk),count(firstext));
						for ( uint64_t i = 0; i < len; ++i )
							weight = std::min(weight,getNode(stretchLinks[stretchLinksStart+i])->freq);
						#endif
					}

					#if defined(AWEIGHT)
					uint64_t aweight = 0;
					for ( uint64_t i = 0; i < len; ++i )
						aweight += getNode(stretchLinks[stretchLinksStart+i])->freq;
					#endif

					Stretch const stretch(
						first,firstext,last,len,
						#if defined(SWEIGHT)
						weight,
						#endif
						#if defined(AWEIGHT)
						aweight,
						#endif
						stretchLinksStart
					);
					maxstretchlength = std::max(maxstretchlength,len);

					if ( debug )
					{
						std::cerr << "refk=" << decode(refk) << " loop=" << (stretch.first==stretch.last) << " ";
						printStretch(stretch,std::cerr); // << " len " << len << " weight " << weight << " stretch " << stretch
					}

					stretches.push(stretcho,stretch);
				}
			}
		}

		#if ! defined(NDEBUG)
		for ( uint64_t i = 1; i < stretcho; ++i )
			assert ( stretches[i-1].first <= stretches[i].first );
		#endif

		if ( debug )
			std::cerr << "stretches computed" << std::endl;
	}

	#if 0
	// get frequency of predecessors in a 64 bit word (16 bits each symbol, top bits A, bottom bits T)
	uint64_t getPredecessors(uint64_t const v) const
	{
		// prefix for successors
		uint64_t const masked = (v >> 2)&m;
		// result vector
		uint64_t r = 0;

		#if 0
		assert ( k );
		#endif
		unsigned int const symshift = (k-1)<<2;

		for ( uint64_t i = 0; i < 4; ++i )
		{
			// shift
			r <<= 16;
			// get range of successor
			uint64_t const freq = count(masked | (i << symshift));

			// if range exceeds 0xFFFFULL then store 0xFFFFULL
			if ( freq & (~(0xFFFFULL)) )
				r |= 0xFFFFULL;
			// otherwise store frequency
			else
				r |= freq;
		}

		return r;
	}
	#endif

	// get number of predecessors for word v
	uint64_t getNumPredecessors(uint64_t const v) const
	{
		// prefix for successors
		uint64_t const masked = (v >> 2)&m;
		// result vector
		uint64_t r = 0;

		#if 0
		assert ( kmersize );
		#endif
		unsigned int const symshift = (kmersize-1)<<2;

		for ( uint64_t i = 0; i < 4; ++i )
			r += (count(masked | (i << symshift) ) != 0);

		return r;
	}

	// print successors for kmer v and their frequency
	std::ostream & printSuccessors(std::ostream & out, uint64_t const v) const
	{
		Links L;
		getActiveSuccessors(v,L);

		for ( uint64_t z = 0; z < L.size(); ++z )
		{
			uint64_t const freq = L.getFreq(z);
			uint64_t const sym = L.getSym(z);
			if ( freq )
				out << decode ( ((v << 2) & m) | sym ) << "\t" << freq << std::endl;
		}

		return out;
	}

	/**
	 * r: successor list as computed by getSuccessors
	 **/
	uint64_t getMaxFreqSucc(uint64_t r) const
	{
		uint64_t s = 3;
		uint64_t c = (r & 0xFFFFULL);

		for ( int i = 2; i >= 0; --i )
		{
			r >>= 16;
			uint64_t const lc = (r & 0xFFFFULL);
			if ( lc > c )
			{
				s = i;
				c = lc;
			}
		}

		return s;
	}

	/**
	 * r: predeccor list as computed by getPredecessors
	 **/
	uint64_t getMaxFreqPred(uint64_t r) const
	{
		return getMaxFreqSucc(r);
	}

	void stretchesUnique()
	{
		Stretch * a = stretches.begin();
		Stretch * e = a + stretcho;
		std::sort(a,e);
		Stretch * i = std::unique(a,e);
		stretcho = i-a;

		uint64_t l = 0;
		uint64_t o = 0;
		while ( l < stretcho )
		{
			uint64_t h = l+1;

			while (
				h < stretcho &&
				stretches[h].first == stretches[l].first &&
				stretches[h].ext == stretches[l].ext
			)
				++h;

			stretches[o++] = stretches[l];

			l = h;
		}

		stretcho = o;
	}


	void computeFeasibleKmerPositions(OffsetLikely const & offsetLikely, double const thres)
	{
		// feasible positions vector for k-mers
		// output pointer for Afeaspos, Acfeasposo
		Afeasposo = 0;
		Acfeasposo = 0;
		maxsupto = 0;

		for ( uint64_t i = 0; i < numnodes; ++i )
		{
			Node & node = nodes[i];
			node.feaspos = Afeasposo;
			node.cfeaspos = Acfeasposo;

			uint64_t const kmer = node.v;

			uint64_t const pfrom = offsetLikely.getSupportLow(node.plow);
			uint64_t const pto   = offsetLikely.getSupportHigh(node.phigh);

			assert ( (pfrom == 0) || (node.plow >= offsetLikely.DPnorm[pfrom-1].size()) );
			assert ( (pto == offsetLikely.DPnorm.size()) ||  (node.phigh < offsetLikely.DPnorm[pto].firstsign) );

			maxsupto = std::max(maxsupto,pto);

			// for ( uint64_t p = 0; p < offsetLikely.size(); ++p )
			for ( uint64_t p = pfrom; p < pto; ++p )
			{
				double const weight = getKmerPositionWeight(kmer,p,offsetLikely);

				if ( weight >= thres )
				{
					// std::cerr << "forward feasible " << p << " weight " << weight << " " << decode(kmer) << std::endl;
					Afeaspos.push(Afeasposo,std::pair<uint64_t,double>(p,weight));
				}
			}

			uint64_t const cpfrom = offsetLikely.getSupportLow(node.cplow);
			uint64_t const cpto   = offsetLikely.getSupportHigh(node.cphigh);
			maxsupto = std::max(maxsupto,cpto);

			assert ( (cpfrom == 0) || (node.cplow >= offsetLikely.DPnorm[cpfrom-1].size()) );
			assert ( (cpto == offsetLikely.DPnorm.size()) ||  (node.cphigh < offsetLikely.DPnorm[cpto].firstsign) );

			for ( uint64_t p = cpfrom; p < cpto; ++p )
			{
				double const weight = getKmerReversePositionWeight(kmer,p,offsetLikely);
				if ( weight >= thres )
				{
					// std::cerr << "reverse feasible " << p << " weight " << weight << " " << decode(kmer) << std::endl;
					Acfeaspos.push(Acfeasposo,std::pair<uint64_t,double>(p,weight));
				}
			}


			node.numfeaspos = Afeasposo - node.feaspos;
			node.numcfeaspos = Acfeasposo - node.cfeaspos;
		}
	}

	void computeFeasibleStretchPositions()
	{
		uint64_t Astretchfeaso = 0;
		uint64_t Acstretchfeaso = 0;

		for ( uint64_t i = 0; i < stretcho; ++i )
		{
			// get stretch
			Stretch & stretch = stretches[i];
			// length of stretch
			uint64_t const len = stretches[i].len;

			{
				// initialise pointers
				stretch.feasposO = Astretchfeaso;
				stretch.feasposL = 0;

				// iterate over stretch links (low to high)
				for ( uint64_t j = 0; j < len; ++j )
				{
					// get link word
					uint64_t const w = stretchLinks [ stretch.stretchO + j ];
					//get node
					Node const * node = getNode(w);
					assert ( node );

					uint64_t poff = len-j-1;

					// iterate over feasible positions for kmer
					for ( uint64_t k = 0; k < node->numfeaspos; ++k )
					{
						std::pair<uint64_t,double> const & FP = Afeaspos [ node->feaspos + k ];

						uint64_t const p = FP.first + poff;

						// push weight
						Afeasbuck[p].second->push(Afeasbuck[p].first,FP.second);
						markBV.set(p);
					}
				}

				uint64_t prank = markBV.getRank();

				uint64_t zz = 0;
				while ( prank-- )
				{
					zz = markBV.next1(zz);
					markBV.erase(zz);

					if ( (Afeasbuck[zz].first == len) && (zz >= len-1) )
					{
						// std::cerr << "zz=" << zz << std::endl;
						double weight = 0.0;

						libmaus2::autoarray::AutoArray<double> const & A = *(Afeasbuck[zz].second);

						for ( uint64_t i = 0; i < len; ++i )
							weight += A[i];

						Astretchfeas.push ( Astretchfeaso,
							StretchFeasObject(
								// stretch start position
								zz - (len-1),
								// stretch weight for position
								weight,
								// weight of first kmer for position
								A[0],
								A[len-1]
							)
						);

						// std::cerr << "forward feasible " << zz - (len-1) << " weight " << weight << std::endl;

						stretch.feasposL += 1;
					}

					Afeasbuck[zz].first = 0;
				}
			}

			{
				// initialise pointers
				stretch.cfeasposO = Acstretchfeaso;
				stretch.cfeasposL = 0;

				// iterate over stretch links (high to low)
				for ( uint64_t jj = 0; jj < len; ++jj )
				{
					uint64_t const j = len - jj - 1;

					// get link word
					uint64_t const w = stretchLinks [ stretch.stretchO + j ];
					//get node
					Node const * node = getNode(w);
					assert ( node );

					// offset to get to the last kmer on the stretch
					uint64_t poff = len-jj-1;

					// iterate over feasible positions for kmer
					for ( uint64_t k = 0; k < node->numcfeaspos; ++k )
					{
						std::pair<uint64_t,double> const & FP = Acfeaspos [ node->cfeaspos + k ];

						uint64_t const p = FP.first + poff;

						#if 0
						bool const ok = (p < Afeasbuck.size());
						if ( ! ok )
						{
							std::cerr << "fail p=" << p << " FP.first=" << FP.first << " poff=" << poff << " len=" << len << std::endl;
							std::cerr << "maxkmerpos=" << maxkmerpos << std::endl;
							std::cerr << "maxstretchlength=" << maxstretchlength << std::endl;
						}
						#endif

						// push weight
						Afeasbuck[p].second->push(Afeasbuck[p].first,FP.second);
						markBV.set(p);
					}
				}

				uint64_t prank = markBV.getRank();

				uint64_t zz = 0;
				while ( prank-- )
				{
					zz = markBV.next1(zz);
					markBV.erase(zz);

					if ( (Afeasbuck[zz].first == len) && (zz >= len-1) )
					{
						// std::cerr << "zz=" << zz << std::endl;
						double weight = 0.0;

						libmaus2::autoarray::AutoArray<double> const & A = *(Afeasbuck[zz].second);

						for ( uint64_t i = 0; i < len; ++i )
							weight += A[i];

						// std::cerr << "reverse feasible " << zz - (len-1) << " weight " << weight << std::endl;

						Acstretchfeas.push ( Acstretchfeaso, StretchFeasObject(zz - (len-1),weight,A[0],A[len-1]) );
						stretch.cfeasposL += 1;
					}

					Afeasbuck[zz].first = 0;
				}

				// check order of positions
				for ( uint64_t i = 1; i < stretch.cfeasposL; ++i )
					assert ( Acstretchfeas [ stretch.cfeasposO + i - 1 ].p < Acstretchfeas [ stretch.cfeasposO + i ].p );
			}
		}
	}

	#if 0
	double getStretchLinkWeight(Stretch const & A, Stretch const & B)
	{
		assert ( A.len );
		uint64_t const shift = A.len-1;

		for ( uint64_t i = 0; i < Afeasbuck.size(); ++i )
			assert ( Afeasbuck[i].first == 0 );

		assert ( B.first == A.last );

		// push shifted positions of first stretch
		StretchFeasObject const * Pa = Astretchfeas.begin() + A.feasposO;
		StretchFeasObject const * Pe = Pa + A.feasposL;
		for ( ; Pa != Pe; ++Pa )
		{
			uint64_t const po = Pa->p + shift;
			Afeasbuck[po].second->push(Afeasbuck[po].first,Pa->w);
			markBV.set(po);
		}
		Pa = Astretchfeas.begin() + B.feasposO;
		Pe = Pa + B.feasposL;
		for ( ; Pa != Pe; ++Pa )
		{
			uint64_t const po = Pa->p;
			Afeasbuck[po].second->push(Afeasbuck[po].first,Pa->w);
			markBV.set(po);
		}

		uint64_t prank = markBV.getRank();

		uint64_t zz = 0;
		double weight = 0.0;
		while ( prank-- )
		{
			zz = markBV.next1(zz);
			markBV.erase(zz);

			assert ( Afeasbuck[zz].first <= 2 );

			if ( Afeasbuck[zz].first == 2 )
			{
				// std::cerr << "zz=" << zz << std::endl;

				libmaus2::autoarray::AutoArray<double> const & A = *(Afeasbuck[zz].second);
				double const lweight = std::min(A[0],A[1]);
				weight = std::max(weight,lweight);
			}

			Afeasbuck[zz].first = 0;
		}

		return weight;
	}
	#endif

	double getReverseStretchLinkWeight(Stretch const & A, Stretch const & B)
	{
		assert ( A.len );
		uint64_t const shift = B.len-1;

		for ( uint64_t i = 0; i < Afeasbuck.size(); ++i )
			assert ( Afeasbuck[i].first == 0 );

		assert ( B.first == A.last );

		// push shifted positions of first stretch
		StretchFeasObject const * Pa = Acstretchfeas.begin() + B.cfeasposO;
		StretchFeasObject const * Pe = Pa + B.cfeasposL;
		for ( ; Pa != Pe; ++Pa )
		{
			uint64_t const po = Pa->p + shift;
			Afeasbuck[po].second->push(Afeasbuck[po].first,Pa->w);
			markBV.set(po);
		}
		Pa = Acstretchfeas.begin() + A.cfeasposO;
		Pe = Pa + A.cfeasposL;
		for ( ; Pa != Pe; ++Pa )
		{
			uint64_t const po = Pa->p;
			Afeasbuck[po].second->push(Afeasbuck[po].first,Pa->w - Pa->wf);
			markBV.set(po);
		}

		uint64_t prank = markBV.getRank();

		uint64_t zz = 0;
		double weight = 0.0;
		while ( prank-- )
		{
			zz = markBV.next1(zz);
			markBV.erase(zz);

			assert ( Afeasbuck[zz].first <= 2 );

			if ( Afeasbuck[zz].first == 2 )
			{
				// std::cerr << "zz=" << zz << std::endl;

				libmaus2::autoarray::AutoArray<double> const & A = *(Afeasbuck[zz].second);
				double const lweight = A[0] + A[1];
				weight = std::max(weight,lweight);
			}

			Afeasbuck[zz].first = 0;
		}

		return weight;
	}

	void computeStretchLinks()
	{
		reverseStretchLinksO = 0;

		// iterate over stretches
		for ( uint64_t i = 0; i < stretcho; ++i )
		{
			// get stretch
			Stretch & stretch = stretches[i];

			// get link candidates
			Stretch ref; ref.first = stretch.last;
			Stretch * a = stretches.begin();
			Stretch * e = stretches.begin()+stretcho;
			std::pair<Stretch *,Stretch *> ER = std::equal_range(a,e,ref,StretchesFirstComparator());

			for ( Stretch * p = ER.first; p != ER.second; ++p )
			{
				double const rweight = getReverseStretchLinkWeight(stretch,*p);

				#if 0
				double const weight = getStretchLinkWeight(stretch,*p);

				if ( (weight == 0 && rweight != 0) || (weight != 0 && rweight == 0) )
				{
					std::cerr << "weight=" << weight << " rweight=" << rweight << std::endl;
					printStretch(stretch,std::cerr);
					printStretch(*p,std::cerr);
				}
				#endif

				uint64_t const linkid = p - a;
				if ( rweight >= 1e-1 )
					reverseStretchLinks.push(reverseStretchLinksO,std::pair<uint64_t,uint64_t>(linkid,i));
			}
		}

		std::sort(reverseStretchLinks.begin(),reverseStretchLinks.begin()+reverseStretchLinksO);
	}

	double getPairScore(Path const & P, ReversePath const & RP) const
	{
		uint64_t const laststretchidP = AP [ P.off + P.len - 1 ];
		Stretch const & laststretchP = stretches[laststretchidP];
		assert ( laststretchP.len );
		assert ( P.pos >= (laststretchP.len-1) );
		uint64_t const laststretchPpos = P.pos - (laststretchP.len-1);
		StretchFeasObject const * SFO = getCachedStretchPositionWeight(laststretchidP, laststretchPpos);

		if ( SFO )
			return
				P.weight + RP.weight - SFO->wl;
		else
			return
				P.weight + RP.weight;
	}

	ScoreInterval getPrimaryScoreInterval(uint64_t const left, uint64_t const right, Path const & P) const
	{
		assert ( left != right );

		uint64_t const m = (*ARWR_RMQ)(left,right-1);

		for ( uint64_t i = left; i < right; ++i )
			assert ( (i == m) || (ARW[i] < ARW[m]) );

		ScoreInterval SI(left,right,m,getPairScore(P,ARP[m]),P);

		return SI;
	}

	bool nextScoreInterval(ScoreInterval & S)
	{
		uint64_t const v = ARW[S.current];

		if ( v )
		{
			uint64_t const u = ARWWT->rpv(S.left,S.right,v-1);

			if ( u == ::std::numeric_limits<uint64_t>::max() )
				return false;
			else
			{
				S.current = ARWWT->select(u,0);
				S.weight = getPairScore(S.P,ARP[S.current]);
				return true;
			}
		}
		else
		{
			return false;
		}
	}

	bool endsOn(std::string const & A, std::string const & B)
	{
		return A.size() >= B.size() && A.substr(A.size()-B.size()) == B;
	}

	void prepareTraverse(
		bool const checkpredecessors,
		bool const mergestretches,
		uint64_t const first, uint64_t const last,
		int64_t const
			#if 0
			lmin
			#endif
			,
		int64_t const lmax,
		bool const debug
	)
	{
		computeStretches(checkpredecessors, debug);

		if ( debug )
			std::cerr << "splitStretches(first)" << std::endl;
		splitStretches(first);
		if ( debug )
			std::cerr << "splitStretches(second)" << std::endl;
		splitStretches(last);
		if ( debug )
			std::cerr << "splitUnique()" << std::endl;
		stretchesUnique();
		if ( debug )
			std::cerr << "mergeStretches()" << std::endl;
		if ( mergestretches )
			mergeStretches(first,last);

		setupFeasBuckets(std::max(maxkmerpos,maxsupto) + maxstretchlength);

		computeFeasibleStretchPositions();

		computeStretchLinks();

		APRo = 0;
		ARPo = 0;

		for ( uint64_t i = 0; i < ARPHo; ++i )
			ARPH[i]->clear();

		if ( getNode(last) )
		{
			RPST.push(
				ReversePath(
					0 /* len */,
					0 /* off */,
					0 /* pos */,
					ReversePath::getDefaultWeight() /* weight */,
					last,
					kmersize /* baselen */
					#if defined(AWEIGHT)
					,getNode(last)->freq /* aweight */
					#endif
				)
			);
		}

		while ( !RPST.empty() )
		{
			ReversePath RP = RPST.top();
			RPST.pop();

			#if 0
			bool const pr = endsOn(
				"GCCGCTTTACAGGCAGGCGACAACCCTTCAGGCCCGCCAATCAGTAGACTGACGTCGCGACCATCCAGCTTCCAGCGTTCCAGCTCAGCGGC",
				printReversePathStringOnly(RP)
			);

			if ( pr )
				std::cerr << "handling " << printReversePathStringOnly(RP) << std::endl;
			#endif

			uint64_t const srcbaselen = RP.baselen;


			while ( expect_false(!(srcbaselen < ARPHo)) )
			{
				rph_ptr_type theap(new rph_type(REVERSE_PATH_HEAP_SIZE));
				ARPH.push(ARPHo,theap);
			}
			assert ( srcbaselen < ARPHo );

			if ( ARPH[srcbaselen]->full() )
			{
				#if 0
				if ( pr )
				{
					std::cerr << "heap is full" << std::endl;
				}
				#endif

				assert ( ! ARPH[srcbaselen]->empty() );

				if ( RP.weight <= ARPH[srcbaselen]->top().weight )
				{
					#if 0
					if ( pr )
					{
						std::cerr << "dropping " << RP.weight << " <= " << ARPH[srcbaselen]->top().weight << std::endl;

						std::vector < ReversePath > VV;
						while ( !ARPH[srcbaselen]->empty() )
							VV.push_back(ARPH[srcbaselen]->pop());

						for ( uint64_t i = 0; i < VV.size(); ++i )
							std::cerr << printReversePathStringOnly(VV[i]) << " " << VV[i].weight << std::endl;

						for ( uint64_t i = 0; i < VV.size(); ++i )
							ARPH[srcbaselen]->push(VV[i]);
					}
					#endif
					// std::cerr << "dropping " << RP << " dominated by " << ARPH[srcbaselen]->top() << std::endl;
					continue;
				}
				else
				{
					ARPH[srcbaselen]->pop();
				}
			}

			assert ( ! ARPH[srcbaselen]->full() );
			ARPH[srcbaselen]->push(RP);

			ARP.push(ARPo,RP);

			#if 0
			printReversePath(RP,std::cerr);
			#endif

			if ( RP.len == 0 )
			{
				for ( uint64_t i = 0; i < stretcho; ++i )
					if ( stretches[i].last == last )
					{
						ReversePath RPE = extendReversePath(RP, i);

						if ( checkReversePathFeasiblePosition(RPE) )
							RPST.pushBump(RPE);
					}
			}
			else if ( RP.baselen < static_cast<uint64_t>((lmax+1)/2) )
			{
				uint64_t const laststretchid = APR [ RP.linkoff + RP.len - 1 ];

				struct PairFirstComparator
				{
					bool operator()(std::pair<uint64_t,uint64_t> const & A, std::pair<uint64_t,uint64_t> const & B) const
					{
						return A.first < B.first;
					}
				};

				std::pair <
					std::pair<uint64_t,uint64_t> const *,
					std::pair<uint64_t,uint64_t> const *
				> EP = std::equal_range(
					reverseStretchLinks.begin(),
					reverseStretchLinks.begin()+reverseStretchLinksO,
					std::pair<uint64_t,uint64_t>(laststretchid,0),
					PairFirstComparator()
				);
				for ( std::pair<uint64_t,uint64_t> const * p = EP.first; p != EP.second; ++p )
				{
					assert ( p->first == laststretchid );
					uint64_t const extend = p->second;

					#if 0
					std::cerr << "extend by" << std::endl;
					printStretch(stretches[extend],std::cerr);
					#endif

					ReversePath RPE = extendReversePath(RP, extend);

					#if 0
					std::cerr << "extended" << std::endl;
					printReversePath(RPE,std::cerr);
					#endif

					if ( checkReversePathFeasiblePosition(RPE) )
						RPST.pushBump(RPE);
				}
			}
		}

		// reverse ReversePath objects by (first,baselen)
		std::sort(ARP.begin(),ARP.begin()+ARPo);

		ARWT.ensureSize(ARPo);
		for ( uint64_t i = 0; i < ARPo; ++i )
			// ARWT[i] = std::pair<double,uint64_t>(ARP[i].getAverageWeight(),i);
			ARWT[i] = std::pair<double,uint64_t>(ARP[i].weight,i);

		std::sort(ARWT.begin(),ARWT.begin() + ARPo);

		ARW.ensureSize(ARPo);
		ARWR.ensureSize(ARPo);
		for ( uint64_t i = 0; i < ARPo; ++i )
		{
			ARW [ ARWT[i].second ] = i;
			ARWR [ ARWT[i].second ] = ARPo - i - 1;
		}

		libmaus2::rmq::QuickDynamicRMQ < uint64_t const * >::unique_ptr_type tARWR_RMQ(
			new libmaus2::rmq::QuickDynamicRMQ < uint64_t const * >(ARWR.begin(),ARPo)
		);
		ARWR_RMQ = UNIQUE_PTR_MOVE(tARWR_RMQ);

		wt_ptr_type twt(new wt_type(ARW.begin(),ARPo));
		ARWWT = UNIQUE_PTR_MOVE(twt);

		#if 0
		ARWR_RMQ->regressionTest();
		#endif

		#if 0
		for ( uint64_t i = 0; i < ARPo; ++i )
			std::cerr << ARW[i] << ";";
		std::cerr << std::endl;
		for ( uint64_t i = 0; i < ARPo; ++i )
			std::cerr << ARP[i].getAverageWeight() << ";";
		std::cerr << std::endl;
		#endif

		#if 0
		std::cerr << std::string(80,'B') << std::endl;
		for ( uint64_t i = 0; i < ARPo; ++i )
			printReversePath(ARP[i],std::cerr);
		std::cerr << std::string(80,'E') << std::endl;
		#endif
	}

	/**
	 * trivially traverse graph and compute consensus
	 *
	 * returns true if we reached the last frequent kmer, false otherwise
	 **/
	bool traverseTrivial()
	{
		conso = 0;

		// maximum kmer for position 0
		uint64_t const first = maxForPos(0);
		// maximum kmer for last position
		uint64_t const last = maxLastWord();

		prepareTraverse(false /* checkpredecessors  */,false /* merge stretches */,first,last,0 /* lmin */,0 /* lmax */,false /* debug */);

		for ( uint64_t i = 0; i < stretcho; ++i )
		{
			if ( stretches[i].first == first && stretches[i].last == last )
			{
				// printStretch(stretches[i],std::cerr);

				consPushWord(first);

				for ( uint64_t j = 1; j < stretches[i].len; ++j )
				{
					uint64_t const w = stretchLinks[stretches[i].stretchO + j];
					consPush(libmaus2::fastx::remapChar(w & 3));
				}

				return true;
			}
		}

		return false;
	}

	double getKmerPositionWeight(uint64_t const kmer, uint64_t const p, OffsetLikely const & offsetLikely) const
	{
		if ( p >= offsetLikely.size() )
		{
			return 0;
		}

		Node const * node = getNode(kmer);

		if ( node )
		{
			DotProduct const & DP = offsetLikely.DPnormSquare[p];
			uint64_t const start = node->pfostart;
			uint64_t const len = node->pfosize;
			PosFreq const * p = PF.begin() + start;
			PosFreq const * pe = p + len;

			while ( p != pe && p->pos < DP.firstsign )
				++p;

			uint64_t const e = (DP.firstsign + DP.V.size());
			double dprr = 0.0;
			for ( ; p != pe && p->pos < e; ++p )
				dprr += p->freq * DP.V[ p->pos - DP.firstsign ];

			return dprr;
		}
		else
		{
			return 0.0;
		}
	}

	double getKmerReversePositionWeight(uint64_t const kmer, uint64_t const p, OffsetLikely const & offsetLikely) const
	{
		if ( p >= offsetLikely.size() )
		{
			return 0;
		}

		Node const * node = getNode(kmer);

		if ( node )
		{
			DotProduct const & DP = offsetLikely.DPnormSquare[p];
			uint64_t const start = node->cpfostart;
			uint64_t const len = node->cpfosize;
			PosFreq const * p = RPF.begin() + start;
			PosFreq const * pe = p + len;

			while ( p != pe && p->pos < DP.firstsign )
				++p;

			uint64_t const e = (DP.firstsign + DP.V.size());
			double dprr = 0.0;
			for ( ; p != pe && p->pos < e; ++p )
				dprr += p->freq * DP.V[ p->pos - DP.firstsign ];

			return dprr;
		}
		else
		{
			return 0.0;
		}
	}

	StretchFeasObject const * getCachedStretchPositionWeight(uint64_t const stretchid, uint64_t const p) const
	{
		StretchFeasObject const * a = Astretchfeas.begin() + stretches[stretchid].feasposO;
		StretchFeasObject const * e = a + stretches[stretchid].feasposL;
		// std::pair<uint64_t,double> const * c = std::lower_bound(a,e,std::pair<uint64_t,double>(P.pos,0.0),PairFirstComparator());
		while ( a != e && a->p < p )
			++a;

		if ( a != e && a->p == p )
			return a;
		else
			return 0;
	}

	StretchFeasObject const * getCachedStretchReversePositionWeight(uint64_t const stretchid, uint64_t const p) const
	{
		StretchFeasObject const * a = Acstretchfeas.begin() + stretches[stretchid].cfeasposO;
		StretchFeasObject const * e = a + stretches[stretchid].cfeasposL;
		// std::pair<uint64_t,double> const * c = std::lower_bound(a,e,std::pair<uint64_t,double>(P.pos,0.0),PairFirstComparator());
		while ( a != e && a->p < p )
			++a;

		if ( a != e && a->p == p )
			return a;
		else
			return 0;
	}

	Path copyPath(Path const & P)
	{
		uint64_t const o = APo;

		for ( uint64_t i = 0; i < P.len; ++i )
		{
			uint64_t const j = P.off + i;
			uint64_t const v = AP[j];
			AP.push(APo,v);
		}

		return Path(P.len,o,P.pos,P.weight,P.baselen
			#if defined(AWEIGHT)
			,P.aweight
			#endif
		);
	}

	ReversePath copyReversePath(ReversePath const & P)
	{
		uint64_t const o = APRo;

		for ( uint64_t i = 0; i < P.len; ++i )
		{
			uint64_t const j = P.linkoff + i;
			uint64_t const v = APR[j];
			APR.push(APRo,v);
		}

		return ReversePath(P.len, o, P.pos, P.weight, P.front, P.baselen
			#if defined(AWEIGHT)
			,P.aweight
			#endif
		);
	}

	#if defined(AWEIGHT)
	void checkPathAvgWeight(Path const & NP)
	{
		uint64_t caweight = 0;
		for ( uint64_t i = 0; i < NP.len; ++i )
		{
			Stretch const & stretch = stretches [ AP [ NP.off + i ] ];

			if ( i == 0 )
				caweight += getNode ( stretchLinks [ stretch.stretchO + 0 ] )->freq;

			for ( uint64_t j = 1; j < stretch.len; ++j )
				caweight += getNode ( stretchLinks [ stretch.stretchO + j ] )->freq;
		}

		assert ( caweight == NP.aweight );
	}
	#endif

	Path extendPath(Path const P, uint64_t const stretchid)
	{
		Path NP = copyPath(P);
		// extension weight
		// double eweight = getStretchPositionWeight(stretches[stretchid],P.pos,offsetLikely);
		StretchFeasObject const * SFO = getCachedStretchPositionWeight(stretchid, P.pos);

		AP.push(APo,stretchid);

		// update stretch length
		NP.len += 1;

		if ( NP.len == 1 )
		{
			NP.baselen = stretches[stretchid].len+kmersize-1;
			#if defined(AWEIGHT)
			NP.aweight = stretches[stretchid].aweight;
			#endif
			NP.weight = SFO ? SFO->w : 0;
		}
		else
		{
			NP.baselen += stretches[stretchid].len-1;
			#if defined(AWEIGHT)
			NP.aweight += stretches[stretchid].aweight -
				getNode(stretches[stretchid].first)->freq;
			#endif

			#if 0
			assert ( NP.len > 1 );
			assert ( AP [ APo - 1 ] == stretchid );
			uint64_t const prevlastid = AP[APo-2];
			Stretch const & prevlast = stretches[prevlastid];
			assert ( P.pos >= prevlast.len - 1 );
			uint64_t const prevpos = P.pos - (prevlast.len - 1);
			StretchFeasObject const * PSFO = getCachedStretchPositionWeight(prevlastid, prevpos);


			if ( PSFO && SFO )
			{
				std::cerr << "PSFO->wl=" << PSFO->wl << " SFO->wf=" << SFO->wf << std::endl;
			}

			if ( PSFO )
			{
				assert ( NP.weight >= PSFO->wl );
				NP.weight -= PSFO->wl;
			}
			#endif

			if ( SFO )
			{
				//std::cerr << "SFO->w=" << SFO->w << " SFO->wf=" << SFO->wf << std::endl;
				// NP.weight += (SFO->w - SFO->wf);
				NP.weight += SFO->w - SFO->wf;
			}
		}

		// update position
		NP.pos += (stretches[stretchid].len-1);

		#if 0
		uint64_t const stretchhash = stretches[stretchid].hash();
		//NP.blockhash.update(stretchhash,stretchhash);
		#endif

		return NP;
	}

	ReversePath extendReversePath(ReversePath const P, uint64_t const stretchid)
	{
		ReversePath NP = copyReversePath(P);

		StretchFeasObject const * SFO = getCachedStretchReversePositionWeight(stretchid,P.pos);

		APR.push(APRo,stretchid);

		// update stretch length
		NP.len += 1;

		if ( NP.len == 1 )
		{
			NP.baselen = stretches[stretchid].len + kmersize - 1;
			#if defined(AWEIGHT)
			NP.aweight = stretches[stretchid].aweight;
			#endif
			NP.weight = SFO ? SFO->w : 0.0;
		}
		else
		{
			NP.baselen += stretches[stretchid].len - 1;
			#if defined(AWEIGHT)
			NP.aweight += stretches[stretchid].aweight - getNode(stretches[stretchid].first)->freq;
			#endif

			if ( SFO )
				NP.weight += SFO->w - SFO->wf;
		}

		// update position
		NP.pos += stretches[stretchid].len-1;
		NP.front = stretches[stretchid].first;

		#if ! defined(NDEBUG)
		for ( uint64_t i = 1; i < NP.len; ++i )
		{
			uint64_t const firststretchid = APR [ NP.linkoff + i - 1];
			Stretch const & firststretch = stretches[firststretchid];
			uint64_t const secondstretchid = APR [ NP.linkoff + i - 0];
			Stretch const & secondstretch = stretches[secondstretchid];

			assert ( firststretch.first == secondstretch.last );
		}
		#endif

		return NP;
	}

	#if 0
	void compactPathData(std::priority_queue<Path> & PQ)
	{
		// std::cerr << "Compacting path data from " << APo << std::endl;

		std::vector<Path> PQV;
		while ( !PQ.empty() )
		{
			PQV.push_back(PQ.top());
			PQ.pop();
		}

		std::sort(PQV.begin(),PQV.end(),PathOffComparator());

		APo = 0;
		for ( uint64_t i = 0; i < PQV.size(); ++i )
			PQ.push(copyPath(PQV[i]));

		// std::cerr << "To " << APo << std::endl;
	}
	#endif


	bool checkReversePathFeasiblePosition(ReversePath const & RP) const
	{
		if ( RP.len )
		{
			// std::cerr << "checking " << RP << std::endl;

			// id of last stretch on path
			uint64_t const laststretchid = APR [ RP.linkoff + RP.len - 1];
			// last stretch on path
			Stretch const & laststretch = stretches[laststretchid];

			// sanity check
			assert ( RP.pos >= laststretch.len-1 );
			uint64_t const checkpos = RP.pos - (laststretch.len-1);

			for ( uint64_t i = 0; i < laststretch.cfeasposL; ++i )
			{
				StretchFeasObject const & SFO = Acstretchfeas[laststretch.cfeasposO+i];

				if ( SFO.p == checkpos && SFO.w >= 0.5 )
					return true;
			}

			return false;
		}
		else
		{
			return true;
		}
	}

	void printStretch(Stretch const & S, std::ostream & ostr) const
	{
		ostr << S << " first=" << decode(S.first) << " ext=" << decode(S.ext) << " last=" << decode(S.last) << " seq=" << decode(S.first);

		for ( uint64_t i = 1; i < S.len; ++i )
			ostr.put ( libmaus2::fastx::remapChar(stretchLinks[S.stretchO + i] & 3) );

		ostr.put('(');
		ostr.put('F');

		for ( uint64_t i = 0; i < S.feasposL; ++i )
		{
			ostr << "(" << Astretchfeas[S.feasposO+i].p << "," << Astretchfeas[S.feasposO+i].w << ")";
		};

		ostr.put(')');

		ostr.put('(');
		ostr.put('R');

		for ( uint64_t i = 0; i < S.cfeasposL; ++i )
		{
			ostr << "(" << Acstretchfeas[S.cfeasposO+i].p << "," << Acstretchfeas[S.cfeasposO+i].w << ")";
		};

		ostr.put(')');

		ostr.put('\n');
	}

	void printNode(Node const & S, std::ostream & ostr) const
	{
		Links L;
		getActiveSuccessors(S.v,L);

		ostr << S << " word=" << decode(S.v) << " succ ";

		for ( uint64_t i = 0; i < L.size(); ++i )
			ostr << decode ( ((S.v << 2) & m) | L.getSym(i) ) << ";";

		ostr << " feaspos ";
		for ( uint64_t i = 0; i < S.numfeaspos; ++i )
			ostr << "(" << Afeaspos[S.feaspos + i].first << "," << Afeaspos[S.feaspos+i].second << ");";
	}

	std::string printNode(Node const & S) const
	{
		std::ostringstream ostr;
		printNode(S,ostr);
		return ostr.str();
	}

	void printPath(Path const & P, std::ostream & ostr) const
	{
		ostr << "Path " << P << " length " << P.pos+kmersize << " " << std::endl;

		for ( uint64_t i = 0; i < P.len; ++i )
		{
			uint64_t const s = AP[P.off + i - 0];
			printStretch(stretches[s],ostr);
		}

		if ( P.len )
		{
			uint64_t const firststretch = AP[P.off];
			ostr << decode( stretches[firststretch].first );

			for ( uint64_t i = 0; i < P.len; ++i )
			{
				uint64_t const s = AP[P.off + i - 0];
				Stretch const & stretch = stretches[s];
				for ( uint64_t j = 1; j < stretch.len; ++j )
				{
					uint64_t const w = stretchLinks[stretch.stretchO + j];
					ostr.put(libmaus2::fastx::remapChar(w & 3));
				}
			}
		}

		ostr << std::endl;
	}

	uint64_t decodePathPair(
		Path const & P, ReversePath const & RP,
		libmaus2::autoarray::AutoArray<uint8_t> & A,
		uint64_t o = 0
	) const
	{
		assert ( P.len );

		uint64_t const firststretch = AP[P.off];
		// ostr << decode( stretches[firststretch].first );

		consPushWord(stretches[firststretch].first, A, o);

		for ( uint64_t i = 0; i < P.len; ++i )
		{
			uint64_t const s = AP[P.off + i - 0];
			Stretch const & stretch = stretches[s];
			for ( uint64_t j = 1; j < stretch.len; ++j )
				A.push(o,libmaus2::fastx::remapChar(stretchLinks[stretch.stretchO + j] & 3));
		}

		for ( uint64_t ii = 0; ii < RP.len; ++ii )
		{
			uint64_t i = RP.len - ii - 1;

			uint64_t const stretchid = APR [ RP.linkoff + i ];
			Stretch const & stretch = stretches[stretchid];

			for ( uint64_t j = 1; j < stretch.len; ++j )
				A.push(o,libmaus2::fastx::remapChar(stretchLinks[stretch.stretchO + j] & 3));
		}

		return o;
	}

	std::string printReversePath(ReversePath const & P) const
	{
		std::ostringstream ostr;
		printReversePath(P,ostr);
		return ostr.str();
	}

	std::string printReversePathStringOnly(ReversePath const & P) const
	{
		std::ostringstream ostr;
		printReversePathStringOnly(P,ostr);
		return ostr.str();
	}

	void printReversePath(ReversePath const & P, std::ostream & ostr) const
	{
		ostr << P << std::endl;

		for ( uint64_t ii = 0; ii < P.len; ++ii )
		{
			uint64_t const i = P.len-ii-1;
			uint64_t const s = APR[P.linkoff + i];
			printStretch(stretches[s],ostr);
		}

		if ( P.len )
		{
			uint64_t cbaselen = 0;

			for ( uint64_t ii = 0; ii < P.len; ++ii )
			{
				uint64_t i = P.len - ii - 1;

				uint64_t const stretchid = APR [ P.linkoff + i ];
				Stretch const & stretch = stretches[stretchid];

				if ( ii == 0 )
				{
					ostr << decode( stretch.first );
					cbaselen += kmersize;
				}

				for ( uint64_t j = 1; j < stretch.len; ++j )
				{
					uint64_t const w = stretchLinks[stretch.stretchO + j];
					ostr.put(libmaus2::fastx::remapChar(w & 3));
					cbaselen += 1;
				}
			}

			assert ( cbaselen == P.baselen );
		}

		ostr << std::endl;
	}

	void printReversePathStringOnly(ReversePath const & P, std::ostream & ostr) const
	{
		if ( P.len )
		{
			uint64_t cbaselen = 0;

			for ( uint64_t ii = 0; ii < P.len; ++ii )
			{
				uint64_t i = P.len - ii - 1;

				uint64_t const stretchid = APR [ P.linkoff + i ];
				Stretch const & stretch = stretches[stretchid];

				if ( ii == 0 )
				{
					ostr << decode( stretch.first );
					cbaselen += kmersize;
				}

				for ( uint64_t j = 1; j < stretch.len; ++j )
				{
					uint64_t const w = stretchLinks[stretch.stretchO + j];
					ostr.put(libmaus2::fastx::remapChar(w & 3));
					cbaselen += 1;
				}
			}

			assert ( cbaselen == P.baselen );
		}
	}


	/**
	 * traverse graph and compute consensus
	 *
	 * returns true if we reached the last frequent kmer, false otherwise
	 **/
	bool traverse(
		int64_t const lmin, int64_t const lmax,
		std::pair< uint8_t const *, uint64_t> const *
			#if ! defined(CANDIDATE_AVOID_ALIGN)
			MA
			#endif
			,
		uint64_t const
			#if ! defined(CANDIDATE_AVOID_ALIGN)
			MAo
			#endif
			,
		uint64_t const maxfrontpath,
		uint64_t const maxfullpath
	)
	{
		#if defined(TRAVERSE_TIME)
		libmaus2::timing::RealTimeClock rtc;
		#endif
		// std::cerr << "lmin=" << lmin << " lmax=" << lmax << std::endl;

		conso = 0;
		APo = 0;
		ACCo = 0;
		CDH.clear();

		// bool const debug = true;
		bool const debug = false;

		// get kmer with maximum frequency for position p
		#if defined(TRAVERSE_TIME)
		rtc.start();
		#endif
		uint64_t const maxFirstO = maxForPosList(0,maxFirst);
		uint64_t const maxLastO =  maxLastList(maxLast);
		#if defined(TRAVERSE_TIME)
		std::cerr << "maxForPosList/maxLastList " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;
		#endif

		if ( debug )
		{
			uint64_t maxFirstI = 0;
			for ( ; maxFirstI < maxFirstO && maxFirst[maxFirstI].first >= maxFirst[0].first/2; ++maxFirstI )
				std::cerr << "considering first " << maxFirst[maxFirstI].second << " " << decode(maxFirst[maxFirstI].second) << " freq " << maxFirst[maxFirstI].first << std::endl;
			for ( ; maxFirstI < maxFirstO; ++maxFirstI )
				std::cerr << "ignoring first " << maxFirst[maxFirstI].second << " " << decode(maxFirst[maxFirstI].second) << " freq " << maxFirst[maxFirstI].first << std::endl;

			uint64_t maxLastI = 0;
			for ( ; maxLastI < maxLastO && maxLast[maxLastI].first >= maxLast[0].first/2; ++maxLastI )
				std::cerr << "considering last " << maxLast[maxLastI].second << " " << decode(maxLast[maxLastI].second) << " freq " << maxLast[maxLastI].first << std::endl;
			for ( ; maxLastI < maxLastO; ++maxLastI )
				std::cerr << "ignoring last " << maxLast[maxLastI].second << " " << decode(maxLast[maxLastI].second) << " freq " << maxLast[maxLastI].first << std::endl;
		}

		#if 0
		std::cerr << "maxFirstO=" << maxFirstO << " maxLastO=" << maxLastO << std::endl;
		#endif

		uint64_t const threscnt = 3;
		uint64_t const thresden = 4;
		uint64_t const firstthres = maxFirstO ? (maxFirst[0].first*threscnt)/thresden : 0;
		uint64_t const lastthres = maxLastO ? (maxLast[0].first*threscnt)/thresden : 0;

		for ( uint64_t maxFirstI = 0; maxFirstI < maxFirstO && maxFirst[maxFirstI].first >= firstthres; ++maxFirstI )
			for ( uint64_t maxLastI = 0; maxLastI < maxLastO && maxLast[maxLastI].first >= lastthres; ++maxLastI )
			{
				if ( debug )
					std::cerr << "first freq " << maxFirst[maxFirstI].first <<  " word " << decode(maxFirst[maxFirstI].second) << " last freq " << maxLast[maxLastI].first <<  " word " << decode(maxLast[maxLastI].second) << std::endl;

				// maximum kmer for position 0
				uint64_t first = maxFirst[maxFirstI].second; // maxForPos(0);
				// maximum kmer for last position
				uint64_t last = maxLast[maxLastI].second;

				// bool const debug = (first == 11965); // || (first == 11951);

				if ( debug )
					std::cerr << "traverse first=" << first << " " << decode(first) << " freq " << maxFirst[maxFirstI].first << " last=" << last << " " << decode(last) << " freq " << maxLast[maxLastI].first << std::endl;

				libmaus2::timing::RealTimeClock preprtc; preprtc.start();
				prepareTraverse(true /* check predecessors */,false /* merge stretches */, first,last,lmin,lmax,debug /* debug */);
				double const preptime = preprtc.getElapsedSeconds();

				if ( debug )
					std::cerr << "prepareTraverse done" << std::endl;

				// std::cerr << "traverse first=" << decode(v) << std::endl;

				//std::cerr << "maxForPos()=" << decode(v) << std::endl;
				//std::cerr << "maxlast()=" << decode(maxLast()) << std::endl;

				libmaus2::timing::RealTimeClock treertc; treertc.start();
				// std::priority_queue<Path> PQ;
				PQ.clear();
				APo = 0;

				#if 0
				uint64_t APoclean = 0;
				#endif

				for ( uint64_t i = 0; i < stretcho; ++i )
					if ( stretches[i].first == first )
						PQ.pushBump(extendPath(Path(),i));

				SIQ.clear();

				for ( uint64_t i = 0; i < Apathheapo; ++i )
					Apathheap[i]->clear();

				for ( uint64_t numfrontpath = 0 ; (!PQ.empty()) && (numfrontpath < maxfrontpath) ; )
				{
					Path const P = PQ.top();
					PQ.pop();

					while ( expect_false(!(P.baselen < Apathheapo)) )
					{
						libmaus2::util::FiniteSizeHeap<Path,PathWeightHeapComparator>::shared_ptr_type NH(
							new libmaus2::util::FiniteSizeHeap<Path,PathWeightHeapComparator>(PATH_HEAP_SIZE)
						);
						Apathheap.push(Apathheapo,NH);
					}
					assert ( P.baselen < Apathheapo );
					if ( Apathheap[P.baselen]->full() )
					{
						assert ( ! Apathheap[P.baselen]->empty() );

						if ( P.weight <= Apathheap[P.baselen]->top().weight )
						{
							continue;
						}
						else
						{
							assert ( P.weight > Apathheap[P.baselen]->top().weight );
							Apathheap[P.baselen]->pop();
						}
					}
					assert ( ! Apathheap[P.baselen]->full() );
					Apathheap[P.baselen]->push(P);

					// printPath(P,std::cerr);

					assert ( P.len );

					// length of candidate stretch in bases
					int64_t const candlen = P.pos + kmersize;
					assert ( candlen == static_cast<int64_t>(P.baselen) );

					// reference object
					ReversePath RP;
					RP.front = stretches[AP[P.off + P.len - 1]].last;

					#if 0
					for ( uint64_t i = 1; i < ARPo; ++i )
						assert ( ARP[i-1].front <= ARP[i].front );
					#endif

					std::pair < ReversePath const *, ReversePath const *> RPP = std::equal_range(ARP.begin(),ARP.begin()+ARPo,RP,ReversePathFrontComparator());

					RP.baselen = std::max(lmin + static_cast<int64_t>(kmersize) - candlen,static_cast<int64_t>(0));
					ReversePath const * subRPP =std::lower_bound(RPP.first,RPP.second,RP,ReversePathBaseLenComparator());

					assert ( subRPP == RPP.first  || static_cast<int64_t>(subRPP[-1].baselen) + candlen - kmersize < lmin );
					assert ( subRPP == RPP.second || static_cast<int64_t>(subRPP->baselen) + candlen - kmersize >= lmin );

					RP.baselen = std::max(lmax + static_cast<int64_t>(kmersize) - candlen,static_cast<int64_t>(0));
					ReversePath const * supRPP = std::upper_bound(subRPP,RPP.second,RP,ReversePathBaseLenComparator());
					assert ( supRPP == RPP.second || static_cast<int64_t>(supRPP->baselen) + candlen - kmersize > lmax );

					if ( subRPP != supRPP )
					{
						SIQ.pushBump(getPrimaryScoreInterval(subRPP - ARP.begin(),supRPP - ARP.begin(),P));
						++numfrontpath;
					}

					if ( debug )
					{
						std::cerr << "processing " << P << " length " << P.pos+kmersize << " " << std::endl;
						printPath(P, std::cerr);
					}


					uint64_t const Plaststretchid = AP[P.off + P.len - 1];
					Stretch ref; ref.first = stretches[Plaststretchid].last;
					Stretch * a = stretches.begin();
					Stretch * e = stretches.begin()+stretcho;
					std::pair<Stretch *,Stretch *> ER = std::equal_range(a,e,ref,StretchesFirstComparator());

					if (
						P.baselen < kmersize ||
						(
							static_cast<int64_t>(P.baselen - kmersize) < ((lmax+1)/2)
						)
					)
					{
						for ( Stretch * p = ER.first; p != ER.second; ++p )
						{
							uint64_t const addstretchid = p-a;

							if ( debug )
								std::cerr << "adding stretch " << stretches[addstretchid]
									<< " first=" << decode(stretches[addstretchid].first)
									<< " ext=" << decode(stretches[addstretchid].ext)
									<< " last=" << decode(stretches[addstretchid].last)
									<< std::endl;

							StretchFeasObject const * SFO = getCachedStretchPositionWeight(addstretchid,P.pos);
							double const eweight = SFO ? SFO->w : 0.0;
							double const eweightthres = 0.1;

							if ( eweight > eweightthres )
							{
								Path EP = extendPath(P,addstretchid);

								//std::cerr << "EP.weight=" << EP.weight << " EP.len=" << EP.len << " eweight=" << eweight << std::endl;

								//assert ( EP.weight >= eweight );

								if ( EP.weight > eweightthres && static_cast<int64_t>(EP.pos + kmersize) <= lmax )
									PQ.pushBump(EP);
							}
						}
					}

					#if 0
					APoclean += P.len;
					if ( APoclean >= APo/2 )
					{
						compactPathData(PQ);
						APoclean = 0;
					}
					#endif
				}
				double const treetime = treertc.getElapsedSeconds();

				//std::cerr << std::string(80,'-') << std::endl;

				uint64_t prevo = 0, prevlen = std::numeric_limits<uint64_t>::max();

				libmaus2::timing::RealTimeClock evalrtc; evalrtc.start();
				for ( uint64_t numfullpath = 0; (! SIQ.empty()) && (numfullpath<maxfullpath); ++numfullpath )
				{
					ScoreInterval SI = SIQ.top();
					SIQ.pop();

					// copy and enque next (if any)
					ScoreInterval SIC = SI;
					if ( nextScoreInterval(SIC) )
						SIQ.pushBump(SIC);

					Path const & P = SI.P;
					ReversePath const & RP = ARP[SI.current];
					double const weight = SI.weight;

					if ( CDH.full() )
					{
						assert ( ! CDH.empty() );
						if ( weight <= CDH.top().weight )
							continue;
						else
							CDH.pop();
					}
					assert ( ! CDH.full() );

					uint64_t const consstart = conso;
					conso = decodePathPair(P,RP,Acons,conso);
					uint64_t const conslen = conso - consstart;

					if (
						conslen == prevlen &&
						std::equal(
							Acons.begin()+prevo,
							Acons.begin()+prevo+prevlen,
							Acons.begin()+consstart
						)
					)
						continue;

					prevo = consstart;
					prevlen = conslen;

					// construct candidate object
					CDH.push(ConsensusCandidate(consstart,conslen,weight,0.0 /* error */));
				}
				double const evaltime = evalrtc.getElapsedSeconds();

				if ( debug )
					std::cerr << "preptime " << preptime << " treetime " << treetime << " evaltime=" << evaltime << std::endl;
			}

		// std::cerr << std::string(80,'-') << std::endl;
		CH.clear();
		while ( ! CDH.empty() )
			CH.pushBump(CDH.pop());
		while ( ! CH.empty() )
		{
			ConsensusCandidate CC = CH.pop();

			#if ! defined(CANDIDATE_AVOID_ALIGN)

			// YYY
			#if 1

			double const canderr = getSimpleCandidateError(MA,MAo,
				std::pair<uint8_t const *, uint8_t const *>(Acons.begin()+CC.o,Acons.begin()+CC.o+CC.l)
			);

			CC.error = canderr;

			#else


			int64_t const score = getSumOfPairsCandidateScore(MA,MAo,
				std::pair<uint8_t const *, uint8_t const *>(Acons.begin()+CC.o,Acons.begin()+CC.o+CC.l)
			);

			CC.error = -score;

			// std::cerr << "weight " << CC.weight << " score " << score << " string " << std::string(Acons.begin()+CC.o,Acons.begin()+CC.o+CC.l) << std::endl;

			#endif


			#endif

			ACC.push(ACCo,CC);
		}

		for ( uint64_t i = 1; i < ACCo; ++i )
			assert ( ACC[i-1].weight >= ACC[i].weight );

		#if 0
		bool orderok = true;
		for ( uint64_t i = 1; i < ACCo; ++i )
			if ( ACC[i].error < ACC[0].error )
				orderok = false;
		if ( ! orderok )
			for ( uint64_t i = 0; i < ACCo; ++i )
				std::cerr << i << "\t" << ACC[i] << std::endl;
		#endif

		#if ! defined(CANDIDATE_AVOID_ALIGN)
		std::sort(ACC.begin(),ACC.begin()+ACCo,ConsensusCandidateErrorComparator());
		#endif

		// std::cerr << "ACCo=" << ACCo << std::endl;

		#if 0
		for ( uint64_t i = 0; i < ACCo; ++i )
		{
			std::cerr << ACC[i] << std::endl;
			std::cerr << std::string(Acons.begin()+ACC[i].o,Acons.begin()+ACC[i].o+ACC[i].l) << std::endl;
		}
		#endif

		return (ACCo != 0);
	}

	uint64_t getNumCandidates() const
	{
		return ACCo;
	}

	std::pair<uint8_t const *, uint8_t const *> getCandidate(uint64_t const i) const
	{
		uint8_t const * a = reinterpret_cast<uint8_t const *>(Acons.begin() + ACC[i].o);
		uint8_t const * e = a + ACC[i].l;
		return std::pair<uint8_t const *, uint8_t const *>(a,e);
	}

	double getCandidateWeight(uint64_t const i) const
	{
		return ACC[i].weight;
	}

	libmaus2::autoarray::AutoArray<uint64_t> E;

	struct MatrixElement
	{
		int32_t p;
		int32_t d;
		int32_t c;

		MatrixElement() {}
		MatrixElement(
			int32_t const rp,
			int32_t const rd,
			int32_t const rc
		) : p(rp), d(rd), c(rc) {}

		bool operator<(MatrixElement const & O) const
		{
			if ( p != O.p )
				return p < O.p;
			else if ( d != O.d )
				return d < O.d;
			else
				return c < O.c;
		}
	};
	libmaus2::autoarray::AutoArray<MatrixElement> AME;
	libmaus2::autoarray::AutoArray< std::pair<char,uint64_t> > SOP;

	int64_t getSumOfPairsCandidateScore(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	)
	{
		//ZZZ
		libmaus2::lcs::Aligner & NP = *SNP;

		uint64_t oAME = 0;
		for ( uint64_t i = 0; i < o; ++i )
		{
			NP.align(Ucand.first,Ucand.second-Ucand.first,I[i].first,I[i].second);

			libmaus2::lcs::AlignmentTraceContainer::step_type const * te = NP.getTraceContainer().te;
			libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = NP.getTraceContainer().ta;

			int64_t pu = Ucand.second-Ucand.first;
			int64_t pa = I[i].second;
			int64_t d = 0;

			while ( te != ta )
			{
				libmaus2::lcs::AlignmentTraceContainer::step_type const step = *(--te);

				switch ( step )
				{
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
					case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
						--pu;
						--pa;
						d = 0;
						AME.push(oAME,MatrixElement(pu,d,I[i].first[pa]));
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
						--pa;
						d -= 1;
						AME.push(oAME,MatrixElement(pu,d,I[i].first[pa]));
						break;
					case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
						--pu;
						d = 0;
						AME.push(oAME,MatrixElement(pu,d,'-'));
						break;
					default:
						break;
				}
			}
		}

		std::sort(AME.begin(),AME.begin()+oAME);

		uint64_t eAME = oAME;
		uint64_t l = 0;
		while ( l < oAME )
		{
			uint64_t h = l+1;
			while ( h < oAME && AME[l].p == AME[h].p && AME[l].d == AME[h].d )
				++h;

			uint64_t const n = h-l;
			assert ( n <= o );
			if ( n < o )
			{
				uint64_t const r = o-n;
				for ( uint64_t i = 0; i < r; ++i )
					AME.push(eAME,MatrixElement(AME[l].p,AME[l].d,'-'));
			}

			l = h;
		}

		oAME = eAME;
		std::sort(AME.begin(),AME.begin()+oAME);

		l = 0;
		int64_t gs = 0;
		while ( l < oAME )
		{
			uint64_t h = l+1;
			while ( h < oAME && AME[l].p == AME[h].p && AME[l].d == AME[h].d )
				++h;


			uint64_t oSOP = 0;
			uint64_t ll = l;
			while ( ll < h )
			{
				uint64_t hh = ll+1;
				while ( hh < h && AME[hh].c == AME[ll].c )
					++hh;

				SOP.push(oSOP,std::pair<char,uint64_t>(AME[ll].c,hh-ll));

				ll = hh;
			}

			int64_t s = 0;
			int64_t const score_match = 3;
			int64_t const score_mismatch = 2;
			int64_t const score_gap = 1;
			for ( uint64_t i = 0; i < oSOP; ++i )
			{
				if ( SOP[i].first != '-' )
					s += score_match * (SOP[i].second*(SOP[i].second-1)/2);

				for ( uint64_t j = i+1; j < oSOP; ++j )
				{
					if ( SOP[i].first != '-' && SOP[j].first != '-' )
					{
						s -= score_mismatch * static_cast<int64_t>(SOP[i].second * SOP[j].second);
					}
					else
					{
						s -= score_gap * static_cast<int64_t>(SOP[i].second * SOP[j].second);
					}
				}
			}

			#if 0
			std::cerr << std::string(80,'-') << std::endl;

			std::cerr << "s=" << s << std::endl;

			for ( uint64_t i = l; i < h; ++i )
			{
				std::cerr << "AME(" << AME[i].p << "," << AME[i].d << "," << static_cast<char>(AME[i].c) << std::endl;
			}
			#endif

			gs += s;

			l = h;
		}

		return gs;
	}

	double getSimpleCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	)
	{
		algn->process(Ucand.first,Ucand.second,I,o,E);
		return std::accumulate(E.begin(),E.begin()+o,0ull);
	}

	double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	)
	{
		// algn->process(Ucand.first,Ucand.second,I,o,E);

		#if 0
		std::cerr << "query " << std::string(Ucand.first,Ucand.second) << std::endl;
		for ( uint64_t i = 0; i < o; ++i )
			std::cerr << "seq " << std::string(I[i].first,I[i].first + I[i].second) << std::endl;
		#endif

		libmaus2::lcs::AlignmentStatistics AS;

		for ( uint64_t j = 0; j < o; ++j )
		{
			#if 0
			std::cerr << "aligning\n"
				<< std::string(Ucand.first,Ucand.second) << "\n"
				<< std::string(I[j].first,I[j].first + I[j].second) << "\n";
			#endif

			SNP->align(Ucand.first,Ucand.second-Ucand.first,I[j].first,I[j].second);
			AS += SNP->getTraceContainer().getAlignmentStatistics();

			#if 0
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				Ucand.first,Ucand.second-Ucand.first,I[j].first,I[j].second,
				80,
				SNP.ta,
				SNP.te
			);
			#endif
		}

		double const erate = AS.getErrorRate();

		return erate;
	}

	uint64_t getCandidateErrorU(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	)
	{
		#if 0
		std::cerr << "query " << std::string(Ucand.first,Ucand.second) << std::endl;
		for ( uint64_t i = 0; i < o; ++i )
			std::cerr << "seq " << std::string(I[i].first,I[i].first + I[i].second) << std::endl;
		#endif

		libmaus2::lcs::AlignmentStatistics AS;

		for ( uint64_t j = 0; j < o; ++j )
		{
			#if 0
			std::cerr << "aligning\n"
				<< std::string(Ucand.first,Ucand.second) << "\n"
				<< std::string(I[j].first,I[j].first + I[j].second) << "\n";
			#endif

			SNP->align(Ucand.first,Ucand.second-Ucand.first,I[j].first,I[j].second);
			AS += SNP->getTraceContainer().getAlignmentStatistics();

			#if 0
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				Ucand.first,Ucand.second-Ucand.first,I[j].first,I[j].second,
				80,
				SNP.ta,
				SNP.te
			);
			#endif
		}

		uint64_t const ecount = AS.getEditDistance();

		return ecount;
	}

	double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o,
		uint64_t const id
	)
	{
		std::pair<uint8_t const *, uint8_t const *> Ucand = getCandidate(id);
		return getCandidateError(I,o,Ucand);
	}

	uint64_t getCandidateErrorU(
		std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o,
		uint64_t const id
	)
	{
		std::pair<uint8_t const *, uint8_t const *> Ucand = getCandidate(id);
		return getCandidateErrorU(I,o,Ucand);
	}

	// setup using list of (string,length) pairs
	std::pair<uint64_t,double> checkCandidates(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o)
	{
		if ( getNumCandidates() )
			return std::pair<uint64_t,double>(0,getCandidateError(I,o,0));
		else
			return std::pair<uint64_t,double>(0,std::numeric_limits<double>::max());
	}

	std::pair<uint64_t,uint64_t> checkCandidatesU(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o)
	{
		if ( getNumCandidates() )
			return std::pair<uint64_t,uint64_t>(0,getCandidateErrorU(I,o,0));
		else
			return std::pair<uint64_t,uint64_t>(0,std::numeric_limits<double>::max());
	}
};
#endif
