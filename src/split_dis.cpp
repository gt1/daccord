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

#include <config.h>

#include <ChainSet.hpp>
#include <libmaus2/math/ClusteredBernoulli.hpp>
#include <libmaus2/aio/DebugLineOutputStream.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/dazzler/align/TrueOverlap.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/db/InqualContainer.hpp>
#include <libmaus2/fastx/KmerRepeatDetector.hpp>
#include <libmaus2/lcs/AlignerFactory.hpp>
#include <libmaus2/lcs/AlignmentOneAgainstManyFactory.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/SuffixArrayLCS.hpp>
#include <libmaus2/math/binom.hpp>
#include <libmaus2/math/Convolution.hpp>
#include <libmaus2/math/ipow.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/rank/ERank222B.hpp>
#include <libmaus2/rmq/QuickDynamicRMQ.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/util/PrefixSums.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/wavelet/WaveletTree.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/bambam/BamAccessor.hpp>
#include <libmaus2/clustering/KMeans.hpp>
#include <libmaus2/util/Enumerator.hpp>

// #define PHASE_DEBUG
//#define PHASE_SERIAL
//#define PHASER_HANDLE_SINGLE 534

static std::string getFastAIndexFileName(std::string const & consfn)
{
	std::string const consfnindex = consfn + ".phaser.index";
	return consfnindex;
}

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

static char const * getDefaultSequencingDepth()
{
	return "no default";
}

int generateindex(
	libmaus2::util::ArgParser const & arg,
	std::string const & consfn,
	uint64_t const fastaindexmod
)
{

	std::string const consfnindex = getFastAIndexFileName(consfn);
	std::string const tmpfilebase = getTmpFileBase(arg);

	if (
		! libmaus2::util::GetFileSize::fileExists(consfnindex)
		||
		libmaus2::util::GetFileSize::isOlder(consfnindex,consfn)
	)
	{
		std::string const tmpfn = tmpfilebase + ".faindextmp";
		libmaus2::fastx::FastAReader::enumerateOffsets(consfn,tmpfn,fastaindexmod);
		libmaus2::aio::OutputStreamFactoryContainer::rename(tmpfn,consfnindex);
	}

	return EXIT_SUCCESS;
}

typedef libmaus2::lcs::AlignmentTraceContainer trace_type;

struct TraceTypeInfo
{
	typedef trace_type element_type;
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

struct TraceAllocator
{
	typedef trace_type element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct ReadData
{
	typedef ReadData this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::autoarray::AutoArray<char> A;
	uint64_t l;

	ReadData()
	{
	}
};

struct ReadDataTypeInfo
{
	typedef ReadData element_type;
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

struct ReadDataAllocator
{
	typedef ReadData element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct ReadDecoder
{
	typedef ReadDecoder this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::dazzler::db::DatabaseFile const & DB;
	libmaus2::aio::InputStream::unique_ptr_type bpsfile;
	libmaus2::aio::InputStream::unique_ptr_type idxfile;

	ReadDecoder(libmaus2::dazzler::db::DatabaseFile const & rDB)
	: DB(rDB), bpsfile(DB.openBaseStream()), idxfile(DB.openIndexStream())
	{

	}

	size_t decodeRead(uint64_t const id, libmaus2::autoarray::AutoArray<char> & A) const
	{
		return DB.decodeReadAndRC(*idxfile,*bpsfile,id,A);
	}
};

struct ReadDecoderTypeInfo
{
	typedef ReadDecoder element_type;
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

struct ReadDecoderAllocator
{
	typedef ReadDecoder element_type;
	typedef element_type::shared_ptr_type pointer_type;

	libmaus2::dazzler::db::DatabaseFile const * DB;

	ReadDecoderAllocator(libmaus2::dazzler::db::DatabaseFile const * rDB)
	: DB(rDB)
	{

	}

	pointer_type operator()() const
	{
		return pointer_type(new element_type(*DB));
	}
};

struct OverlapTypeInfo
{
	typedef libmaus2::dazzler::align::Overlap element_type;
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

struct OverlapAllocator
{
	typedef libmaus2::dazzler::align::Overlap element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

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


struct DecodedReadContainer
{
	std::map<uint64_t,uint64_t> M;
	std::vector<ReadData::shared_ptr_type> V;

	libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & readDataFreeList;
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList;

	ReadDecoder::shared_ptr_type decoder;

	DecodedReadContainer(
		libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & rreadDataFreeList,
		libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & rreadDecoderFreeList
	) : readDataFreeList(rreadDataFreeList), readDecoderFreeList(rreadDecoderFreeList), decoder(readDecoderFreeList.get())
	{
	}

	~DecodedReadContainer()
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			readDataFreeList.put(V[i]);
		readDecoderFreeList.put(decoder);
	}

	uint64_t ensurePresent(uint64_t const id)
	{
		std::map<uint64_t,uint64_t>::iterator it = M.find(id);

		if ( it == M.end() )
		{
			ReadData::shared_ptr_type R = readDataFreeList.get();
			R->l = decoder->decodeRead(id,R->A);

			uint64_t const vid = V.size();
			V.push_back(R);
			M[id] = vid;
		}

		it = M.find(id);
		assert ( it != M.end() );

		return it->second;
	}

	char const * getForwardRead(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->A.begin();
	}

	char const * getReverseComplementRead(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->A.begin() + R->l;
	}

	uint64_t getReadLength(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->l;
	}
};

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static uint64_t getDefaultMaxInput()
{
	return 5000;
}

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

template<>
std::string formatRHS(std::string const & description, char const * def)
{
	std::ostringstream ostr;
	if ( def == std::string("no default") )
		ostr << description << " (" << def << ")";
	else
		ostr << description << " (default " << def << ")";
	return ostr.str();
}


#if defined(PHASE_DEBUG)
// load all alignments from a BAM file
std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> readBamAlignments(std::string const & bamfn)
{
	std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> V;
	libmaus2::bambam::BamDecoder bamdec(bamfn);
	libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
	while ( bamdec.readAlignment() )
		V.push_back(algn.sclone());
	return V;
}
#endif

struct AlignerBase
{
	static libmaus2::lcs::Aligner::unique_ptr_type getAligner()
	{
		std::set<libmaus2::lcs::AlignerFactory::aligner_type> const S = libmaus2::lcs::AlignerFactory::getSupportedAligners();

		libmaus2::lcs::AlignerFactory::aligner_type T[] = {
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8,
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8,
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_NP
		};

		for ( uint64_t i = 0; i < sizeof(T)/sizeof(T[0]); ++i )
			if ( S.find(T[i]) != S.end() )
			{
				libmaus2::lcs::Aligner::unique_ptr_type tptr(libmaus2::lcs::AlignerFactory::construct(T[i]));
				return UNIQUE_PTR_MOVE(tptr);
			}

		libmaus2::exception::LibMausException lme;
		lme.getStream() << "DebruijnGraph::getAligner: no suitable class found" << std::endl;
		lme.finish();
		throw lme;
	}
};

/**
 * pile element
 **/
struct MarkedPileElement
{
	// position on consensus
	int64_t aid;
	int64_t apos;
	int64_t apre;
	// symbol in A/B read
	char sym;
	// sequence ID
	int64_t chain;
	int64_t seq;
	// position on source sequence
	int64_t spos;
	int64_t spre;

	MarkedPileElement() {}
	MarkedPileElement(
		int64_t const raid,
		int64_t const rapos,
		int64_t const rapre,
		char const rsym,
		int64_t const rchain,
		int64_t const rseq,
		int64_t const rspos,
		int64_t const rspre
	) : aid(raid), apos(rapos), apre(rapre), sym(rsym), chain(rchain), seq(rseq), spos(rspos),spre(rspre) {}

	bool operator<(MarkedPileElement const & P) const
	{
		if ( aid != P.aid )
			return aid < P.aid;
		if ( apos != P.apos )
			return apos < P.apos;
		else if ( apre != P.apre )
			return apre < P.apre;
		else if ( sym != P.sym )
			return sym < P.sym;
		else if ( chain != P.chain )
			return chain < P.chain;
		else
			return seq < P.seq;
	}

	std::ostream & toString(std::ostream & out) const
	{
		return out << "MarkedPileElement(aid=" << aid << ",apos=" << apos << ",apre=" << apre << ",sym=" << sym << ",chain=" << chain << ",seq=" << seq << ",spos=" << spos << ",spre=" << spre << ")";
	}

	std::string toString() const
	{
		std::ostringstream ostr;
		toString(ostr);
		return ostr.str();
	}
};

std::ostream & operator<<(std::ostream & out, MarkedPileElement const & A)
{
	return A.toString(out);
}

struct IdentityPositionMapper
{
	int64_t operator()(int64_t const p) const
	{
		return p;
	}
};

struct MapPosPositionMapper
{
	bool inv;
	int64_t rlb;
	int64_t off;

	MapPosPositionMapper() {}
	MapPosPositionMapper(
		bool const rinv,
		int64_t const rrlb,
		int64_t const roff
	) : inv(rinv), rlb(rrlb), off(roff)
	{

	}

	int64_t operator()(int64_t const pos) const
	{
		int64_t r;

		if ( inv )
		{
			r = rlb - (pos+off) - 1;
		}
		else
		{
			r = (pos+off);
		}

		return r;
	}

};

struct ReadFragment
{
	uint64_t from;
	uint64_t to;
	std::string fragment;

	ReadFragment() {}
	ReadFragment(uint64_t const rfrom, uint64_t const rto, std::string const & rfragment)
	: from(rfrom), to(rto), fragment(rfragment)
	{

	}
};

struct SeqToChain
{
	std::vector<int64_t> const & seqtochain;

	SeqToChain(std::vector<int64_t> const & rseqtochain) : seqtochain(rseqtochain) {}

	int64_t operator()(int64_t const seq) const
	{
		if ( seq < 0 )
			return -1;
		else
			return seqtochain[seq];
	}
};

void completeLPV(std::vector<MarkedPileElement> & LPV)
{
	uint64_t o = 0;
	for ( uint64_t l = 0; l < LPV.size(); )
	{
		// interval of same aid,apos,apre
		uint64_t h = l+1;
		while (
			h < LPV.size()
			&&
			LPV[h].aid == LPV[l].aid
			&&
			LPV[h].apos == LPV[l].apos
			&&
			LPV[h].apre == LPV[l].apre
		)
			++h;

		for ( uint64_t i = l+1; i < h; ++i )
		{
			assert ( LPV[i-1].chain <= LPV[i].chain );
		}

		uint64_t ll = l;
		while ( ll < h )
		{
			uint64_t hh = ll+1;
			while ( hh < h && LPV[ll].chain == LPV[hh].chain )
				++hh;

			LPV[o++] = LPV[ll];

			ll = hh;
		}

		l = h;
	}
	assert ( o <= LPV.size() );

	if ( o < LPV.size() )
		LPV.resize(o);

	for ( uint64_t h = LPV.size(); h; )
	{
		// interval of same aid,apos,apre
		uint64_t l = h-1;
		while (
			l
			&&
			LPV[l-1].aid == LPV[l].aid
			&&
			LPV[l-1].apos == LPV[l].apos
			&&
			LPV[l-1].apre == LPV[l].apre
		)
			--l;

		for ( uint64_t i = l+1; i < h; ++i )
			assert ( LPV[i-1].chain < LPV[i].chain );

		h = l;
	}

	uint64_t l0 = 0, h0 = 0;
	libmaus2::autoarray::AutoArray< std::pair<int64_t,int64_t> > ASEQ;
	libmaus2::autoarray::AutoArray< std::pair<int64_t,int64_t> > ASEQt;
	std::vector<MarkedPileElement> TLPV;
	uint64_t aseqf = 0;
	// step through LPV backwards and add missing entries
	// (add D elements in columns where a sequence is active but not present)
	for ( uint64_t h = LPV.size(); h; )
	{
		// interval of same aid,apos,apre
		uint64_t l = h-1;
		while (
			l
			&&
			LPV[l-1].aid == LPV[l].aid
			&&
			LPV[l-1].apos == LPV[l].apos
			&&
			LPV[l-1].apre == LPV[l].apre
		)
			--l;

		// check in interval
		for ( uint64_t i = l; i < h; ++i )
			assert ( LPV[i].aid == LPV[l].aid && LPV[i].apos == LPV[l].apos && LPV[i].apre == LPV[l].apre );
		// check boundaries
		if ( l )
			assert ( LPV[l-1].aid != LPV[l].aid || LPV[l-1].apos != LPV[l].apos || LPV[l-1].apre != LPV[l].apre );
		if ( h < LPV.size() )
			assert ( LPV[h].aid != LPV[h-1].aid || LPV[h].apos != LPV[h-1].apos || LPV[h].apre != LPV[h-1].apre );

		// check sequence order
		for ( uint64_t i = l+1; i < h; ++i )
		{
			bool const ok = LPV[i-1].chain < LPV[i].chain;
			assert ( ok );
		}
		for ( uint64_t i = 1; i < aseqf; ++i )
			assert (
				ASEQ[i-1].first < ASEQ[i].first
				||
				(
					ASEQ[i-1].first == ASEQ[i].first
					&&
					ASEQ[i-1].second < ASEQ[i].second
				)
			);

		if ( LPV[l].apre == 0 )
		{
			l0 = l;
			h0 = h;

			// store sequence ids for this position
			aseqf = 0;
			for ( uint64_t i = l0; i < h0; ++i )
				ASEQ.push(aseqf,std::pair<int64_t,int64_t>(LPV[i].chain,LPV[i].seq));
		}
		else
		{
			uint64_t i0 = 0;
			uint64_t ii = l;
			uint64_t taseqf = 0;

			while ( i0 < aseqf && ii < h )
			{
				if ( ASEQ[i0].first < LPV[ii].chain )
				{
					// std::cerr << "add to LPV missing " << ASEQ[i0] << std::endl;
					TLPV.push_back(MarkedPileElement(LPV[l].aid,LPV[l].apos,LPV[l].apre,'D',ASEQ[i0].first,ASEQ[i0].second,-1,-1));
					ASEQt.push(taseqf,ASEQ[i0]);
					++i0;
				}
				else if ( LPV[ii].chain < ASEQ[i0].first )
				{
					// std::cerr << "add to ASEQ missing " << LPV[ii].seq << std::endl;
					ASEQt.push(taseqf,std::pair<int64_t,int64_t>(LPV[ii].chain,LPV[ii].seq));
					++ii;
				}
				else
				{
					assert ( LPV[ii].chain == ASEQ[i0].first );

					// std::cerr << "already present " << ASEQ[i0] << std::endl;

					ASEQt.push(taseqf,ASEQ[i0]);

					++ii;
					++i0;
				}
			}
			while ( i0 < aseqf )
			{
				// std::cerr << "add to LPV missing " << ASEQ[i0] << std::endl;
				TLPV.push_back(MarkedPileElement(LPV[l].aid,LPV[l].apos,LPV[l].apre,'D',ASEQ[i0].first,ASEQ[i0].second,-1,-1));
				ASEQt.push(taseqf,ASEQ[i0]);
				++i0;
			}
			while ( ii < h )
			{
				// std::cerr << "add to ASEQ missing " << LPV[ii].seq << std::endl;
				ASEQt.push(taseqf,std::pair<int64_t,int64_t>(LPV[ii].chain,LPV[ii].seq));
				++ii;
			}

			ASEQ.swap(ASEQt);
			aseqf = taseqf;
		}

		h = l;
	}

	// add new elements
	for ( uint64_t i = 0; i < TLPV.size(); ++i )
		LPV.push_back(TLPV[i]);
	// restore sort order after adding elements
	std::sort(LPV.begin(),LPV.end());
}

template<typename iterator_a, typename iterator_b, typename position_mapper_type>
void fillLPV(
	uint64_t const cid,
	libmaus2::lcs::AlignmentTraceContainer const & ATC, // NPcons
	std::vector<MarkedPileElement> & LPV,
	iterator_a ub,       // RC.getForwardRead(ita[0].aread)
	iterator_b SCO,      // SCO.c_str()
	uint64_t const bstart, // afirst
	uint64_t const bend,   // alast+1
	uint64_t const cstart, // 0
	uint64_t const cend,   // SCO.size()
	int64_t const seqid,
	position_mapper_type const & posmap,
	SeqToChain const & STC
)
{
	// fill LPV for A read
	libmaus2::lcs::AlignmentTraceContainer::step_type const * uta = ATC.ta;
	libmaus2::lcs::AlignmentTraceContainer::step_type const * ute = ATC.te;
	uint64_t cipos = cend; // SCO.size();
	uint64_t aipos = bend; // alast+1;
	int64_t cidif = 0;
	int64_t aidif = 0;
	// char const * ub = RC.getForwardRead(ita[0].aread);

	while ( ute != uta )
	{
		libmaus2::lcs::AlignmentTraceContainer::step_type const step = *(--ute);

		switch ( step )
		{
			case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
			case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
			{
				--cipos;
				--aipos;
				cidif = 0;
				aidif = 0;

				if ( step == libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH )
					assert ( ub[aipos] == SCO[cipos] );

				MarkedPileElement PE(cid,cipos,cidif,ub[aipos],STC(seqid),seqid,posmap(aipos),aidif);
				LPV.push_back(PE);

				break;
			}
			case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
			{
				--cipos;
				cidif = 0;
				--aidif;

				MarkedPileElement PE(cid,cipos,cidif,'D',STC(seqid),seqid,posmap(aipos),aidif);
				LPV.push_back(PE);

				break;
			}
			case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
			{
				--aipos;
				cidif -= 1;
				aidif = 0;

				MarkedPileElement PE(cid,cipos,cidif,ub[aipos],STC(seqid),seqid,posmap(aipos),aidif);
				LPV.push_back(PE);

				break;
			}
			default:
			{
				break;
			}
		}
	}

	assert ( cipos == cstart );
	assert ( aipos == bstart );
}


template<int _n0, int _n1>
struct PhaseQueueEntryMulti
{
	static int const n0 = _n0;
	static int const n1 = _n1;

	int64_t aseq;
	std::pair<int64_t,int64_t> VIa;
	std::pair<int64_t,int64_t> VIb;

	int64_t I0[n0];
	int64_t I1[n1];

	std::vector<int64_t> getIndices0() const
	{
		std::vector<int64_t> V;
		V.push_back(aseq);
		for ( unsigned int i = 0; i < n0; ++i )
			V.push_back(I0[i]);
		return V;
	}

	std::vector<int64_t> getIndices1() const
	{
		std::vector<int64_t> V;
		for ( unsigned int i = 0; i < n1; ++i )
			V.push_back(I1[i]);
		return V;
	}

	std::vector<int64_t> getChains0(std::vector<MarkedPileElement> const & LPV) const
	{
		std::vector<int64_t> V;
		for ( unsigned int i = 0; i < n0; ++i )
			V.push_back(LPV[I0[i]].chain);
		return V;
	}

	std::vector<int64_t> getChains1(std::vector<MarkedPileElement> const & LPV) const
	{
		std::vector<int64_t> V;
		for ( unsigned int i = 0; i < n1; ++i )
			V.push_back(LPV[I1[i]].chain);
		return V;
	}

	std::vector<int64_t> getIndices() const
	{
		std::vector<int64_t> V;
		V.push_back(aseq);
		for ( unsigned int i = 0; i < n0; ++i )
			V.push_back(I0[i]);
		for ( unsigned int i = 0; i < n1; ++i )
			V.push_back(I1[i]);
		return V;
	}

	#if defined(PHASE_DEBUG)
	bool splitOk(
		std::vector < uint64_t > const & chainlengths,
		std::vector < uint64_t > const & chains,
		BamAlignmentContainer & BAC,
		std::vector<MarkedPileElement> const & LPV,
		libmaus2::dazzler::align::Overlap const * const ita,
		libmaus2::dazzler::align::TrueOverlap & TO
	) const
	{
		for ( uint64_t z = 0; z < n1; ++z )
		{
			uint64_t const chainid = LPV[I1[z]].chain;

			for ( uint64_t chainsubid = chainlengths[chainid]; chainsubid < chainlengths[chainid+1]; ++chainsubid )
			{
				uint64_t const y = chains[chainsubid];
				libmaus2::dazzler::align::Overlap const & OVL = ita[y];
				bool const to = TO.trueOverlap(OVL,BAC[OVL.aread],BAC[OVL.bread]);
				if ( to )
					return false;
			}
		}

		return true;
	}

	bool splitOkRef(
		BamAlignmentContainer & BAC,
		std::vector<MarkedPileElement> const & LPV,
		libmaus2::dazzler::align::Overlap const * const ita
	) const
	{
		std::vector<int64_t> U0;
		std::vector<int64_t> U1;
		U0.push_back(aseq);
		for ( int i = 0; i < n0; ++i )
			U0.push_back(I0[i]);
		for ( int i = 0; i < n1; ++i )
			U1.push_back(I1[i]);

		for ( uint64_t i = 0; i < U0.size(); ++i )
		{
			MarkedPileElement const & MPE0 = LPV[U0[i]];
			int64_t const reada = MPE0.seq < 0 ? ita[0].aread : ita[MPE0.seq].bread;

			for ( uint64_t j = 0; j < U1.size(); ++j )
			{
				MarkedPileElement const & MPE1 = LPV[U1[j]];
				assert ( MPE1.seq >= 0 );
				int64_t const readb = ita[MPE1.seq].bread;

				if ( BAC[reada].getRefID() == BAC[readb].getRefID() )
					return false;
			}
		}

		return true;
	}

	void check(
		BamAlignmentContainer & BAC,
		std::vector<MarkedPileElement> const & LPV,
		libmaus2::dazzler::align::Overlap const * const ita
	) const
	{
		std::vector<int64_t> const V = getIndices();

		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			MarkedPileElement const & MPE = LPV[V[i]];
			int64_t const seqid = MPE.seq;
			int64_t const readid = seqid < 0 ? ita[0].aread : ita[seqid].bread;
			libmaus2::bambam::BamAlignment const & algn = BAC[readid];

			if ( MPE.spre == 0 )
			{
				bool const correct = algn.baseMatchPrime(MPE.spos);
				std::cerr << MPE.sym << (correct ? "+" : "-");
			}
			else
			{
				std::cerr << MPE.sym << "?";
			}
		}

		std::cerr << std::endl;
	}

	bool sameZero(
		BamAlignmentContainer & BAC,
		std::vector<MarkedPileElement> const & LPV,
		libmaus2::dazzler::align::Overlap const * const ita
	) const
	{
		if ( n0 )
		{
			int64_t ref = BAC[ita[LPV[I0[0]].seq].aread].getRefID();
			for ( int i = 0; i < n0; ++i )
				if ( BAC[ita[LPV[I0[i]].seq].bread].getRefID() != ref )
					return false;
			return true;
		}
		else
		{
			return true;
		}
	}

	bool sameOne(
		BamAlignmentContainer & BAC,
		std::vector<MarkedPileElement> const & LPV,
		libmaus2::dazzler::align::Overlap const * const ita
	) const
	{
		if ( n1 )
		{
			int64_t ref = BAC[ita[LPV[I1[0]].seq].bread].getRefID();
			for ( int i = 0; i < n1; ++i )
				if ( BAC[ita[LPV[I1[i]].seq].bread].getRefID() != ref )
					return false;
			return true;
		}
		else
		{
			return true;
		}
	}
	#endif

	bool sameChain(PhaseQueueEntryMulti<n0,n1> const & O, std::vector<MarkedPileElement> const & LPV) const
	{
		for ( int i = 0; i < n0; ++i )
			if ( LPV[I0[i]].chain != LPV[O.I0[i]].chain  )
				return false;
		for ( int i = 0; i < n1; ++i )
			if ( LPV[I1[i]].chain != LPV[O.I1[i]].chain )
				return false;

		return true;
	}

	static std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > getConsInterval(
		std::vector < uint64_t > const & chainlengths,
		std::vector < uint64_t > const & chains,
		std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > const & Gconsintervals,
		std::vector<uint64_t> const & VC
	)
	{
		// iterate over fragments
		std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > GIV;
		for ( uint64_t i = 0; i < Gconsintervals.size(); ++i )
		{
			std::vector < libmaus2::math::IntegerInterval<int64_t> > IV;

			{
				uint64_t const chainid = VC[0];
				for ( uint64_t j = chainlengths[chainid]; j < chainlengths[chainid+1]; ++j )
				{
					uint64_t const ovlid = chains[j];

					if ( !Gconsintervals[i][ovlid].isEmpty() )
						IV.push_back(Gconsintervals[i][ovlid]);
				}
			}
			for ( uint64_t z = 1; z < VC.size(); ++z )
			{
				std::vector < libmaus2::math::IntegerInterval<int64_t> > NIV;
				uint64_t const chainid = VC[z];

				for ( uint64_t j = chainlengths[chainid]; j < chainlengths[chainid+1]; ++j )
				{
					uint64_t const ovlid = chains[j];

					if ( !Gconsintervals[i][ovlid].isEmpty() )
					{
						for ( uint64_t k = 0; k < IV.size(); ++k )
						{
							libmaus2::math::IntegerInterval<int64_t> const C = Gconsintervals[i][ovlid].intersection(IV[k]);
							if ( !C.isEmpty() )
								NIV.push_back(C);
						}
					}
				}

				IV = NIV;
			}

			GIV.push_back(IV);
		}

		return GIV;

	}

	std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > getConsInterval(
		std::vector < uint64_t > const & chainlengths,
		std::vector < uint64_t > const & chains,
		std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > const & Gconsintervals,
		std::vector<MarkedPileElement> const & LPV
	) const
	{
		assert ( n0 + n1 );

		// chain IDs
		std::vector<uint64_t> VC;
		for ( uint64_t i = 0; i < n0; ++i )
			VC.push_back(LPV[I0[i]].chain);
		for ( uint64_t i = 0; i < n1; ++i )
			VC.push_back(LPV[I1[i]].chain);

		return getConsInterval(chainlengths,chains,Gconsintervals,VC);
	}

	PhaseQueueEntryMulti() {}
	PhaseQueueEntryMulti(
		int64_t const raseq,
		std::pair<int64_t,int64_t> const & rVIa,
		std::pair<int64_t,int64_t> const & rVIb
	) : aseq(raseq), VIa(rVIa), VIb(rVIb)
	{
		if ( VIa.second-VIa.first >= n0 && VIb.second-VIb.first >= n1 )
		{
			for ( int i = 0; i < n0; ++i )
				I0[i] = VIa.first + i;
			for ( int i = 0; i < n1; ++i )
				I1[i] = VIb.first + i;
		}
		else
		{
			for ( int i = 0; i < n0; ++i )
				I0[i] = VIa.second;
			for ( int i = 0; i < n1; ++i )
				I1[i] = VIb.second;
		}
	}

	std::string toString() const
	{
		std::ostringstream ostr;
		ostr << "PhaseQueueEntryMulti(";

		ostr << "I0=[";
		for ( int i = 0; i < n0; ++i )
			ostr << I0[i] << ";";
		ostr << "]";
		ostr << ";";
		ostr << "I1=[";
		for ( int i = 0; i < n1; ++i )
			ostr << I1[i] << ";";
		ostr << "]";
		ostr << ")";

		return ostr.str();
	}

	std::string toString(std::vector<MarkedPileElement> const & LPV, libmaus2::dazzler::align::Overlap const * const ita) const
	{
		std::ostringstream ostr;
		ostr << "PhaseQueueEntryMulti(";

		ostr << "I0=[";
		for ( int i = 0; i < n0; ++i )
			ostr << LPV[I0[i]].chain << ":" << LPV[I0[i]].seq << ":" << ita[LPV[I0[i]].seq].bread << ";";
		ostr << "]";
		ostr << ";";
		ostr << "I1=[";
		for ( int i = 0; i < n1; ++i )
			ostr << LPV[I1[i]].chain << ":" << LPV[I1[i]].seq << ":" << ita[LPV[I1[i]].seq].bread << ";";
		ostr << "]";
		ostr << ")";

		return ostr.str();
	}

	std::string toString(
		std::vector<MarkedPileElement> const & LPV, libmaus2::dazzler::align::Overlap const * const ita,
		BamAlignmentContainer & BAC
		) const
	{
		std::ostringstream ostr;
		ostr << "PhaseQueueEntryMulti(";

		ostr << "I0=[";
		for ( int i = 0; i < n0; ++i )
			ostr << LPV[I0[i]].chain << ":" << LPV[I0[i]].seq << ":" << ita[LPV[I0[i]].seq].bread << ":" << BAC[ita[LPV[I0[i]].seq].bread].getRefID() << ";";
		ostr << "]";
		ostr << ";";
		ostr << "I1=[";
		for ( int i = 0; i < n1; ++i )
			ostr << LPV[I1[i]].chain << ":" << LPV[I1[i]].seq << ":" << ita[LPV[I1[i]].seq].bread << ":" << BAC[ita[LPV[I0[i]].seq].bread].getRefID() << ";";
		ostr << "]";
		ostr << ")";

		return ostr.str();
	}

	bool next()
	{
		for ( int ii = 0; ii < n1; ++ii )
		{
			int const i = n1-ii-1;

			if ( I1[i]+1+ii < VIb.second )
			{
				++I1[i];

				for ( int j = i+1; j < n1; ++j )
					I1[j] = I1[i] + (j-i);

				return true;
			}
		}
		for ( int ii = 0; ii < n0; ++ii )
		{
			int const i = n0-ii-1;

			if ( I0[i]+1+ii < VIa.second )
			{
				++I0[i];

				for ( int j = i+1; j < n0; ++j )
					I0[j] = I0[i] + (j-i);

				for ( int j = 0; j < n1; ++j )
					I1[j] = VIb.first+j;

				return true;
			}
		}

		return false;
	}
};

template<int n0, int n1>
struct PhaseQueueEntryMultiComparator
{
	std::vector<MarkedPileElement> const * LPV;
	PhaseQueueEntryMultiComparator() : LPV(0) {}
	PhaseQueueEntryMultiComparator(std::vector<MarkedPileElement> const * rLPV) : LPV(rLPV) {}

	bool operator()(PhaseQueueEntryMulti<n0,n1> const & A, PhaseQueueEntryMulti<n0,n1> const & B) const
	{
		for ( int i = 0; i < n0; ++i )
			if ( (*LPV)[A.I0[i]].chain != (*LPV)[B.I0[i]].chain )
				return (*LPV)[A.I0[i]].chain < (*LPV)[B.I0[i]].chain;
		for ( int i = 0; i < n1; ++i )
			if ( (*LPV)[A.I1[i]].chain != (*LPV)[B.I1[i]].chain )
				return (*LPV)[A.I1[i]].chain < (*LPV)[B.I1[i]].chain;
		return false;
	}
};


void phaseChains(
	std::ostream & err,
	std::vector < ReadFragment > const & fragments,
	std::vector < uint64_t > const & chainlengths,
	std::vector < uint64_t > const & chains,
	uint64_t const numchains,
	// overlaps
	libmaus2::dazzler::align::Overlap const * const ita,
	libmaus2::dazzler::align::Overlap const * const ite,
	// reads
	DecodedReadContainer & RC,
	DecodedReadContainer & RC2,
	// trace free list
	libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> & traceFreeList,
	// aligner
	libmaus2::lcs::Aligner & NP,
	// trace point space for overlaps
	int64_t const tspace,
	libmaus2::dazzler::align::AlignmentWriter & AW,
	libmaus2::dazzler::align::AlignmentWriter & AK,
	uint64_t const /* d */,
	uint64_t const dthres,
	double const crate,
	double const drate,
	uint64_t const phasethres,
	bool const checkdisagreement,
	uint64_t const dcthres
	#if defined(PHASE_DEBUG)
		,
		BamAlignmentContainer & BAC,
		libmaus2::dazzler::db::DatabaseFile const & DB,
		libmaus2::dazzler::db::DatabaseFile const & /* DB2 */
	#endif
)
{
	libmaus2::timing::RealTimeClock rtc;
	rtc.start();

	if ( ita == ite )
		return;

	std::vector<int64_t> seqtochain(ite-ita,-1);
	std::vector<int64_t> chainbaselength(numchains,0);
	for ( uint64_t chainid = 0; chainid < numchains; ++chainid )
	{
		uint64_t const len = chainlengths[chainid+1]-chainlengths[chainid];
		for ( uint64_t j = 0; j < len; ++j )
		{
			uint64_t const z = chains[chainlengths[chainid]+j];
			seqtochain[z] = chainid;

			bool const ok = ita [ z ] . bread == ita [ chains[chainlengths[chainid]+0] ].bread;

			if ( ! ok )
			{
				std::cerr << "assertion failure for aread=" << ita[0].aread << " chain " << chainid << std::endl;
			}

			assert ( ok );

			chainbaselength[chainid] += ita[z].path.aepos-ita[z].path.abpos;
		}
	}

	#if defined(PHASE_DEBUG)
	int64_t const aread = ita->aread;
	for ( uint64_t i = 0; i < fragments.size(); ++i )
		err << "[V] aread=" << aread << " fragments[" << i << "]=[" << fragments[i].from << "," << fragments[i].to << ")" << std::endl;
	#endif

	#if defined(PHASE_DEBUG)
	libmaus2::dazzler::align::TrueOverlapStats TOS;
	libmaus2::dazzler::align::TrueOverlap TO(TOS,DB,tspace);
	std::vector<bool> Vtrue(ite-ita,false);
	std::vector<bool> Vctrue(numchains,false);
	std::vector<libmaus2::lcs::AlignmentStatistics> Vstat(ite-ita);
	std::vector<libmaus2::lcs::AlignmentStatistics> Vcstat(numchains);

	for ( libmaus2::dazzler::align::Overlap const * itc = ita; itc != ite; ++itc )
	{
		uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(itc->aread));
		uint8_t const * ub = reinterpret_cast<uint8_t const *>(itc->isInverse() ? RC2.getReverseComplementRead(itc->bread) : RC2.getForwardRead(itc->bread));
		libmaus2::bambam::BamAlignment const & abam = BAC[itc->aread];
		libmaus2::bambam::BamAlignment const & bbam = BAC[itc->bread];

		trace_type::shared_ptr_type Ptrace = traceFreeList.get();
		itc->computeTrace(ua,ub,tspace,*Ptrace,NP);
		bool const trueoverlap = libmaus2::dazzler::align::TrueOverlap::trueOverlap(abam,bbam,*Ptrace,itc->path.abpos,itc->path.bbpos,itc->isInverse(),RC2.getReadLength(itc->bread));
		traceFreeList.put(Ptrace);

		Vtrue[itc-ita] = trueoverlap;

		std::string const aref = abam.getReferenceUsedPrime(itc->path.abpos,itc->path.aepos);
		int64_t bb = itc->path.bbpos;
		int64_t be = itc->path.bepos;
		if ( itc->isInverse() )
		{
			std::swap(bb,be);
			bb = RC2.getReadLength(itc->bread) - bb;
			be = RC2.getReadLength(itc->bread) - be;
		}
		std::string bref = bbam.getReferenceUsedPrime(bb,be);
		if ( itc->isInverse() )
		{
			bref = libmaus2::fastx::reverseComplementUnmapped(bref);
		}

		libmaus2::lcs::NP np;
		np.np(aref.begin(),aref.end(),bref.begin(),bref.end());

		uint64_t clipfront = 0;
		while ( np.ta+clipfront != np.te && np.ta[clipfront] != libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH )
			++clipfront;
		uint64_t clipback = 0;
		while ( np.ta+clipback != np.te && np.te[-static_cast<int64_t>(clipback)-1] != libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH )
			++clipback;

		Vstat[itc-ita] = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(np.ta+clipfront,np.te-clipback);
	}
	for ( uint64_t i = 0; i < numchains; ++i )
	{
		uint64_t const len = chainlengths[i+1]-chainlengths[i];
		bool trueflag = true;
		for ( uint64_t j = 0; j < len; ++j )
		{
			uint64_t const z = chains[chainlengths[i]+j];
			trueflag = trueflag && Vtrue[z];
			Vcstat[i] += Vstat[z];
		}
		Vctrue[i] = trueflag;
	}
	#endif

	SeqToChain STC(seqtochain);
	std::vector<MarkedPileElement> LPV;

	std::vector < std::vector<libmaus2::math::IntegerInterval<int64_t> > > Gconsintervals(fragments.size());

	// iterate over consensus fragments
	for ( uint64_t cid = 0; cid < fragments.size(); ++cid )
	{
		ReadFragment const & F = fragments[cid];
		uint64_t const afirst = F.from;
		uint64_t const alast  = F.to-1;
		std::string const & SCO = F.fragment;


		// consensus start and end
		uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(SCO.c_str());
		uint8_t const * uSCOe = uSCO + SCO.size();

		uint8_t const * areadpart_a = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[0].aread)+afirst);
		uint8_t const * areadpart_e = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[0].aread)+alast+1);

		libmaus2::lcs::NP NPcons;
		libmaus2::lcs::NP NPbcons;
		// align A to consensus
		NPcons.np(areadpart_a,areadpart_e,uSCO,uSCOe);

		#if defined(PHASE_DEEP_DEBUG)
		libmaus2::lcs::AlignmentPrint::printAlignmentLines(
			err,
			areadpart_a,
			areadpart_e-areadpart_a,
			uSCO,
			uSCOe-uSCO,
			80,
			NPcons.ta,
			NPcons.te
		);
		#endif

		// interval on A read
		libmaus2::math::IntegerInterval<int64_t> IA(afirst,alast);

		// regions on consensus B reads align to
		std::vector < libmaus2::math::IntegerInterval<int64_t> > consintervals(ite-ita,libmaus2::math::IntegerInterval<int64_t>::empty());

		// fill LPV for A read
		fillLPV(cid,NPcons,LPV,RC.getForwardRead(ita[0].aread),SCO.c_str(),afirst,alast+1,0/*cstart*/,SCO.size(),-1,IdentityPositionMapper(),STC);

		// fill consintervals, LPV for B reads
		for ( uint64_t chainid = 0; chainid < numchains; ++chainid )
		{
			uint64_t const len = chainlengths[chainid+1]-chainlengths[chainid];
			for ( uint64_t j = 0; j < len; ++j )
			{
				uint64_t const z = chains[chainlengths[chainid]+j];
				libmaus2::dazzler::align::Overlap const * itc = ita + z;

				libmaus2::math::IntegerInterval<int64_t> const IB(itc->path.abpos,itc->path.aepos-1);
				libmaus2::math::IntegerInterval<int64_t> const IC = IB.intersection(IA);

				if ( ! IC.isEmpty() )
				{
					uint64_t const z = itc-ita;

					libmaus2::dazzler::align::Overlap const & OVL = *itc;
					#if defined(PHASE_DEEP_DEBUG)
					err << "[D] ovl " << IA << " " << libmaus2::math::IntegerInterval<int64_t>(OVL.path.abpos,OVL.path.aepos) << std::endl;
					#endif

					// compute alignment A->B
					trace_type::shared_ptr_type Ptrace = traceFreeList.get();
					trace_type & ATC = *Ptrace;
					{
					uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(OVL.aread));
					uint8_t const * ub = reinterpret_cast<uint8_t const *>(OVL.isInverse() ? RC2.getReverseComplementRead(OVL.bread) : RC2.getForwardRead(OVL.bread));
					OVL.computeTrace(ua,ub,tspace,ATC,NP);
					}

					libmaus2::lcs::AlignmentTraceContainer::step_type const * bta = ATC.ta;
					libmaus2::lcs::AlignmentTraceContainer::step_type const * bte = ATC.te;

					libmaus2::lcs::AlignmentTraceContainer::step_type const * cta = NPcons.ta;
					libmaus2::lcs::AlignmentTraceContainer::step_type const * cte = NPcons.te;

					uint64_t bstart = OVL.path.bbpos;
					uint64_t cstart = 0;

					if ( OVL.path.abpos < static_cast<int64_t>(afirst) )
					{
						// advance afirst-OVL.path.abpos on A->B
						std::pair<uint64_t,uint64_t> const P = libmaus2::lcs::AlignmentTraceContainer::advanceA(bta,bte,afirst - OVL.path.abpos);
						assert ( P.first == afirst - OVL.path.abpos );

						std::pair<uint64_t,uint64_t> const SL = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(bta,bta + P.second);
						assert ( SL.first == afirst - OVL.path.abpos );

						// advance respective bases on B
						bstart += SL.second;
						bta += P.second;

						assert ( OVL.path.abpos + SL.first == afirst );
					}
					else if ( static_cast<int64_t>(afirst) < OVL.path.abpos )
					{
						// advance on consensus
						std::pair<uint64_t,uint64_t> const P = libmaus2::lcs::AlignmentTraceContainer::advanceA(cta,cte,OVL.path.abpos-afirst);
						assert ( P.first == OVL.path.abpos - afirst );

						std::pair<uint64_t,uint64_t> const SL = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(cta,cta + P.second);
						assert ( SL.first == OVL.path.abpos - afirst );

						cstart += SL.second;
						cta += P.second;
					}

					// string length available on B
					std::pair<uint64_t,uint64_t> const SLB = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(bta,bte);
					// string length available on consensus
					std::pair<uint64_t,uint64_t> const SLC = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(cta,cte);

					// minimum of the two -> string length on A
					uint64_t const usea = std::min(SLB.first,SLC.first);

					// cut on B
					std::pair<uint64_t,uint64_t> const PB = libmaus2::lcs::AlignmentTraceContainer::advanceA(bta,bte,usea);
					assert ( PB.first == usea );
					std::pair<uint64_t,uint64_t> const PC = libmaus2::lcs::AlignmentTraceContainer::advanceA(cta,cte,usea);
					assert ( PC.first == usea );

					// string length on B
					std::pair<uint64_t,uint64_t> const UB = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(bta,bta+PB.second);
					assert ( UB.first == usea );
					// string length on C
					std::pair<uint64_t,uint64_t> const UC = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(cta,cta+PC.second);
					assert ( UC.first == usea );

					uint8_t const * ub = reinterpret_cast<uint8_t const *>(OVL.isInverse() ? (RC2.getReverseComplementRead(OVL.bread)+bstart) : (RC2.getForwardRead(OVL.bread)+bstart));
					int64_t const rlb = RC2.getReadLength(OVL.bread);

					#if defined(PHASE_DEBUG)
					assert ( rlb == BAC[OVL.bread].getLseq() );
					#endif

					//std::string const breadpart(ub,ub+UB.second);
					//std::string const conspart(uSCO + cstart, uSCO + cstart + UC.second);

					// interval on C
					consintervals[z] = libmaus2::math::IntegerInterval<int64_t>(cstart,cstart + UC.second - 1);

					// compute alignment B->C
					NPbcons.np(ub,ub+UB.second,uSCO + cstart, uSCO + cstart + UC.second);

					// fill LPV for B read
					MapPosPositionMapper MPPM(OVL.isInverse(),rlb,bstart);
					fillLPV(cid,NPbcons,LPV,ub,uSCO,0,UB.second,cstart,cstart + UC.second,z,MPPM,STC);

					traceFreeList.put(Ptrace);
				}
			}
		}

		Gconsintervals[cid] = consintervals;
	}

	// sort LPV vector by (cpos,cipos,seq)
	struct MarkedPileElementSeqComp
	{
		bool operator()(MarkedPileElement const & A, MarkedPileElement const & B) const
		{
			if ( A.aid != B.aid )
				return A.aid < B.aid;
			else if ( A.apos != B.apos )
				return A.apos < B.apos;
			else if ( A.apre != B.apre )
				return A.apre < B.apre;
			else if ( A.chain != B.chain )
				return A.chain < B.chain;
			else
				return A.seq < B.seq;
		}
	};
	std::sort(LPV.begin(),LPV.end(),MarkedPileElementSeqComp());

	#if defined(PHASE_DEBUG)
	// sanity check
	for ( uint64_t i = 0; i < LPV.size(); ++i )
		if ( LPV[i].spre == 0 && LPV[i].seq != -1 )
		{
			bool const inv = ita[LPV[i].seq].isInverse();
			char const r = RC2.getForwardRead(ita[LPV[i].seq].bread)[LPV[i].spos];
			assert ( ( inv ? libmaus2::fastx::invertUnmapped(r) : r ) == LPV[i].sym );
		}
	#endif

	completeLPV(LPV);

	// count number of distinct positions in LPV
	uint64_t zlow = 0;
	uint64_t lppos = 0;
	while ( zlow < LPV.size() )
	{
		uint64_t zhigh = zlow+1;
		while ( zhigh < LPV.size() && LPV[zhigh].aid == LPV[zlow].aid && LPV[zhigh].apos == LPV[zlow].apos && LPV[zhigh].apre == LPV[zlow].apre )
			++zhigh;

		lppos++;

		zlow = zhigh;
	}

	int64_t ihigh = LPV.size();
	uint64_t depth = 0;
	uint64_t ilppos = lppos;

	int const n0 = 2;
	int const n1 = 3;
	libmaus2::util::FiniteSizeHeap<PhaseQueueEntryMulti<n0,n1>,PhaseQueueEntryMultiComparator<n0,n1> > phasemultiFSH(0,PhaseQueueEntryMultiComparator<n0,n1>(&LPV));

	// compute disagreement positions
	// for this we traverse LPV from back to front
	while ( ihigh )
	{
		// consensus position ID
		--ilppos;
		int64_t ilow = --ihigh;

		while ( ilow >= 0 && LPV[ilow].aid == LPV[ihigh].aid && LPV[ilow].apos == LPV[ihigh].apos && LPV[ilow].apre == LPV[ihigh].apre )
			--ilow;

		++ilow;

		// sanity checks
		for ( int64_t i = ilow; i <= ihigh; ++i )
			assert ( LPV[i].aid == LPV[ilow].aid && LPV[i].apos == LPV[ilow].apos && LPV[i].apre == LPV[ilow].apre );
		assert ( ilow == 0 || LPV[ilow-1].aid != LPV[ilow].aid || LPV[ilow-1].apos != LPV[ilow].apos || LPV[ilow-1].apre != LPV[ilow].apre );
		assert ( ihigh+1 == static_cast<int64_t>(LPV.size()) || LPV[ihigh+1].aid != LPV[ihigh].aid || LPV[ihigh+1].apos != LPV[ihigh].apos || LPV[ihigh+1].apre != LPV[ihigh].apre );

		// local depth
		uint64_t const ldepth = ihigh-ilow+1;

		// if true position on consensus then set current depth
		if ( LPV[ilow].apre == 0 )
			depth = ldepth;
		else
		{
			// sanity check whether filling above worked
			assert ( ldepth >= depth );
		}

		// count symbols
		std::pair<uint64_t,uint64_t> C[] =
		{
			std::pair<uint64_t,uint64_t>(0,'A'),
			std::pair<uint64_t,uint64_t>(0,'C'),
			std::pair<uint64_t,uint64_t>(0,'G'),
			std::pair<uint64_t,uint64_t>(0,'T'),
			std::pair<uint64_t,uint64_t>(0,'D'),
			std::pair<uint64_t,uint64_t>(0,0)
		};

		for ( int64_t i = ilow; i <= ihigh; ++i )
		{
			switch ( LPV[i].sym )
			{
				case 'A':
					C[0].first++;
					break;
				case 'C':
					C[1].first++;
					break;
				case 'G':
					C[2].first++;
					break;
				case 'T':
					C[3].first++;
					break;
				case 'D':
					C[4].first++;
					break;
				default:
					C[5].first++;
					break;
			}
		}
		// there should be no symbols outside A,C,G,T,D
		assert ( C[5].first == 0 );

		// this should be empty (pad column with D until we reach the depth)
		for ( uint64_t i = ldepth; i < depth; ++i )
			C[4].first++;

		// sort by frequency
		std::sort(&C[0],&C[sizeof(C)/sizeof(C[0])],std::greater<std::pair<uint64_t,uint64_t> >());

		// disagreement?
		bool const disagreement = C[1].first >= dthres;
		if ( (checkdisagreement && disagreement) || (!checkdisagreement) )
		{
			int64_t tlow = ilow;
			std::vector < std::pair<int64_t,int64_t> > VI;
			int64_t VIa = -1;

			while ( tlow < (ihigh+1) )
			{
				int64_t thigh = tlow+1;

				while ( thigh < (ihigh+1) && LPV[thigh].sym == LPV[tlow].sym )
					++thigh;

				for ( int64_t i = tlow+1; i < thigh; ++i )
				{
					bool const ok = LPV[i-1].chain < LPV[i].chain;
					if (!ok )
					{
						std::cerr << i-1 << " " << LPV[i-1].toString() << " " << ita[LPV[i-1].seq].bread << std::endl;
						std::cerr << ita[LPV[i-1].seq] << std::endl;
						std::cerr << i << " " << LPV[i].toString() << " " << ita[LPV[i].seq].bread << std::endl;
						std::cerr << ita[LPV[i].seq] << std::endl;
						std::cerr << "LPV.size()=" << LPV.size() << std::endl;
					}
					assert ( ok );
				}

				// does this contain the A read?
				for ( int64_t i = tlow; i < thigh; ++i )
					if ( LPV[i].chain == -1 )
						VIa = VI.size();
				// push interval
				VI.push_back(std::pair<int64_t,int64_t>(tlow,thigh));

				tlow = thigh;
			}

			// if A read is in any of the intervals
			if ( VIa >= 0 )
			{
				// it should be in the first position
				assert ( LPV[VI[VIa].first].chain == -1 );

				// interval without A read
				std::pair<int64_t,int64_t> A(VI[VIa].first+1,VI[VIa].second);

				// if interval is sufficiently large
				if ( A.second-A.first >= n0 )
					// check other intervals
					for ( uint64_t i = 0; i < VI.size(); ++i )
						if ( static_cast<int64_t>(i) != VIa && VI[i].second-VI[i].first >= n1 )
						{
							PhaseQueueEntryMulti<n0,n1> PQ(VI[VIa].first,A,VI[i]);
							phasemultiFSH.pushBump(PQ);
						}
			}
		}

		ihigh = ilow;
	}
	// we should have seen every consensus position exactly once
	assert ( ilppos == 0 );


	int64_t prevchain = -1;
	libmaus2::autoarray::AutoArray<uint64_t> DC(numchains,false);
	std::fill(DC.begin(),DC.end(),0ull);

	struct MinReq
	{
		double const crate;
		double const crate6;
		double const drate;
		std::map<uint64_t,uint64_t> M;

		static uint64_t getMod()
		{
			return 50;
		}

		MinReq(double const rcrate, double const rdrate) : crate(rcrate), crate6(::std::pow(crate,6)), drate(rdrate) {}

		uint64_t operator[](uint64_t const i)
		{
			uint64_t const mod = getMod();
			uint64_t j = ((i + mod-1)/mod)*mod;

			if ( M.find(j) == M.end() )
				M[j] = libmaus2::math::Binom::binomRowUpperGmpFloatLimit(
					crate6,
					std::floor(j*drate+0.5),
					512,
					0.995
				);

			return M.find(j)->second;
		}
	};

	// double const drate = 0.01;
	// std::cerr << "drate=" << drate << std::endl;
	MinReq MR(crate,drate);
	uint64_t const cmin = MR[phasethres];
	while ( ! phasemultiFSH.empty() )
	{
		PhaseQueueEntryMulti<n0,n1> const R = phasemultiFSH.top();

		if ( LPV[R.I0[0]].chain != prevchain )
		{
			int64_t const chainid = LPV[R.I0[0]].chain;
			prevchain = chainid;
			err << "[V] chainid=" << chainid << "/" << numchains << std::endl;
		}

		uint64_t c = 0;

		#if defined(COMPUTE_VQ)
		std::vector < PhaseQueueEntryMulti<n0,n1> > VQ;
		#endif

		while (
			(! phasemultiFSH.empty()) && phasemultiFSH.top().sameChain(R,LPV)
		)
		{
			PhaseQueueEntryMulti<n0,n1> P = phasemultiFSH.pop();

			#if defined(COMPUTE_VQ)
			VQ.push_back(P);
			#endif

			if ( P.next() )
				phasemultiFSH.push(P);

			c += 1;
		}

		if ( c >= cmin )
		{
			std::vector < std::vector < libmaus2::math::IntegerInterval<int64_t> > > GCV = R.getConsInterval(chainlengths,chains,Gconsintervals,LPV);
			uint64_t gcsize = 0;
			for ( uint64_t i = 0; i < GCV.size(); ++i )
				for ( uint64_t j = 0; j < GCV[i].size(); ++j )
					gcsize += GCV[i][j].diameter();

			if ( gcsize >= phasethres && c >= MR[gcsize] )
			{
				std::vector<int64_t> const C1 = R.getChains1(LPV);

				for ( uint64_t i = 0; i < C1.size(); ++i )
					DC[C1[i]]++;
			}
		}
	}

	std::vector<bool> phasekeep(ite-ita,true);
	for ( uint64_t chainid = 0; chainid < numchains; ++chainid )
	{
		#if defined(PHASE_DEBUG)
		bool istrue = true;
		#endif

		uint64_t const len = chainlengths[chainid+1]-chainlengths[chainid];
		uint64_t olen = 0;
		for ( uint64_t j = 0; j < len; ++j )
		{
			uint64_t const z = chains[chainlengths[chainid]+j];
			#if defined(PHASE_DEBUG)
			istrue = istrue && Vtrue[z];
			#endif
			olen += ita[z].path.aepos-ita[z].path.abpos;
		}

		err << "DC[" << chainid << "]=" << DC[chainid];
		#if defined(PHASE_DEBUG)
		err << " istrue=" << istrue;
		err << " dif=" << Vcstat[chainid].getErrorRate();
		#endif
		err << " " << olen << std::endl;

		if ( DC[chainid] >= dcthres )
			for ( uint64_t j = 0; j < len; ++j )
				phasekeep [ chains[chainlengths[chainid]+j] ] = false;
	}

	std::ostringstream keepostr;
	std::ostringstream dropostr;
	keepostr << "[V] keep " << ita[0].aread << " ";
	dropostr << "[V] drop " << ita[0].aread << " ";
	std::map<uint64_t,uint64_t> keepset;
	std::map<uint64_t,uint64_t> dropset;

	for ( uint64_t i = 0; i < phasekeep.size(); ++i )
		if ( phasekeep[i] )
		{
			AW.put(ita[i]);
			keepset[ita[i].bread]++;
		}
		else
		{
			AK.put(ita[i]);
			dropset[ita[i].bread]++;
		}
	for ( std::map<uint64_t,uint64_t>::const_iterator ita = keepset.begin(); ita != keepset.end(); ++ita )
	{
		keepostr << ita->first;
		if ( ita->second > 1 )
			keepostr << '(' << ita->second << ')';
		keepostr << ';';
	}
	err << keepostr.str() << '\n';
	for ( std::map<uint64_t,uint64_t>::const_iterator ita = dropset.begin(); ita != dropset.end(); ++ita )
	{
		dropostr << ita->first;
		if ( ita->second > 1 )
			dropostr << '(' << ita->second << ')';
		dropostr << ';';
	}
	err << dropostr.str() << '\n';

	#if defined(PHASE_DEBUG)
	BAC.clear();
	#endif

	err << "[V] read id " << ita[0].aread << " time " << rtc << std::endl;
}

struct OverlapPosComparator
{
	bool operator()(libmaus2::dazzler::align::Overlap const & lhs, libmaus2::dazzler::align::Overlap const & rhs) const
	{
		return lhs.path.abpos < rhs.path.abpos;
	}
};

void patternToUpper(libmaus2::fastx::Pattern & pattern)
{
	for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
		pattern.spattern[i] = ::toupper(pattern.spattern[i]);

	pattern.pattern = pattern.spattern.c_str();
}

struct ReadAccessor
{
	std::string const consfn;
	std::string const consfnindex;
	uint64_t const mod;

	ReadAccessor(std::string const & rconsfn, std::string const & rconsfnindex, uint64_t const rmod)
	: consfn(rconsfn), consfnindex(rconsfnindex), mod(rmod)
	{

	}

	libmaus2::fastx::FastAReader::pattern_type operator[](uint64_t const i)
	{
		uint64_t const b = i / mod;
		libmaus2::aio::InputStreamInstance IISI(consfnindex);
		IISI.seekg(b * sizeof(libmaus2::fastx::FastInterval));
		libmaus2::fastx::FastInterval FI(IISI);
		libmaus2::fastx::FastAReader FARE(consfn,FI);
		uint64_t o = i - b * mod;
		for ( uint64_t j = 0; j < o; ++j )
			FARE.skipPattern();
		libmaus2::fastx::FastAReader::pattern_type P;
		bool const ok = FARE.getNextPatternUnlocked(P);
		assert ( ok );
		patternToUpper(P);
		return P;
	}

	std::string getName(uint64_t const i)
	{
		libmaus2::fastx::FastAReader::pattern_type const P = (*this)[i];
		return P.sid;
	}

	static uint64_t getReadId(std::string name)
	{
		if ( name.find('/') != std::string::npos )
			name = name.substr(0,name.find('/'));

		std::istringstream istr(name);
		uint64_t id;
		istr >> id;

		assert ( istr && istr.peek() == std::istream::traits_type::eof() );

		assert ( id );

		return id-1;
	}

	uint64_t getReadId(uint64_t const i)
	{
		std::string name = getName(i);
		return getReadId(name);
	}
};

uint64_t getId(std::string sid, libmaus2::fastx::FastAReader::pattern_type const & pattern)
{
	if ( sid.find('/') == std::string::npos )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unparsable read name (no slash) " << sid << std::endl;
		lme.finish();
		throw lme;
	}

	sid = sid.substr(0,sid.find('/'));

	for ( uint64_t i = 0; i < sid.size(); ++i )
		if ( ! ::isdigit(sid[i]) )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unparsable read name (string " << sid << " before slash not numerical) " << pattern.sid << std::endl;
			lme.finish();
			throw lme;
		}

	std::istringstream istr(sid);
	uint64_t id;
	istr >> id;

	assert ( istr.peek() == std::istream::traits_type::eof() );

	assert ( id );

	return id - 1;
}

std::vector<uint64_t> countConsensus(std::string const & consin)
{
	libmaus2::fastx::FastAReader FA(consin);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	std::vector<uint64_t> V;

	while ( FA.getNextPatternUnlocked(pattern) )
	{
		// >2/24/500_15948 A=[500,17349]

		std::string sid = pattern.sid;

		uint64_t const id = getId(sid,pattern);

		while ( ! (id<V.size()) )
			V.push_back(0);
		assert ( id < V.size() );

		V[id]++;
	}

	V.push_back(0);

	::libmaus2::util::PrefixSums::prefixSums(V.begin(),V.end());

	return V;
}

static uint64_t getDefaultPhaseThres()
{
	return 4000;
}

static bool getDefaultCheckDisagreement()
{
	return 0;
}

static double getDefaultDRate()
{
	return 0.01;
}

double poisson(double const lambda, uint64_t const k)
{
	double p = 1.0;
	for ( uint64_t i = 1; i <= k; ++i )
	{
		p *= lambda;
		p /= i;
	}
	p *= ::std::exp(-lambda);
	return p;
}

std::pair<uint64_t,double> poissonlimit(double const lambda, double const limit)
{
	double s = 0.0;
	double p = ::std::exp(-lambda);
	uint64_t k = 0;

	while ( true )
	{
		if ( s + p > limit )
			return std::pair<uint64_t,double>(k,s);

		s += p;
		k += 1;

		p *= lambda;
		p /= k;
	}
}

int phaser(
	libmaus2::util::ArgParser const & arg,
	std::string const & outlasfn,
	std::string const & consin,
	std::string const & inlasfn,
	std::string const & dbname,
	std::string const & db2name
	#if defined(PHASE_DEBUG)
		,
	std::string const & bamfn
	#endif
)
{
	arg.printArgs(std::cerr,std::string("[V] "));

	uint64_t const fastaindexmod = 1;
	generateindex(arg,consin,fastaindexmod);

	std::vector<uint64_t> VCNT = countConsensus(consin);

	std::string const eproffn = arg.uniqueArgPresent("E") ? arg["E"] : (inlasfn + ".eprof");
	libmaus2::lcs::AlignmentStatistics GAS;
	GAS.deserialise(eproffn);
	double const erate = GAS.getErrorRate();
	double const crate = 1.0 - erate;
	double const crate2 = crate * crate;

	if ( !arg.uniqueArgPresent("d") )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] required argument -d (depth) is missing" << std::endl;
		lme.finish();
		throw lme;
	}

	std::cerr << "[V] " << GAS << std::endl;
	std::cerr << "[V] p=" << crate << std::endl;
	std::cerr << "[V] p^2=" << crate2 << std::endl;

	uint64_t const d = arg.getParsedArg<uint64_t>("d");
	uint64_t const phasethres = arg.uniqueArgPresent("p") ? arg.getParsedArg<uint64_t>("p") : getDefaultPhaseThres();
	uint64_t const dthres = libmaus2::math::Binom::binomRowUpperGmpFloatLimit(crate2,d,512,0.99);
	bool const checkdisagreement = arg.uniqueArgPresent("c") ? arg.getParsedArg<uint64_t>("c") : getDefaultCheckDisagreement();
	double const drate = arg.uniqueArgPresent("drate") ? arg.getParsedArg<double>("drate") : getDefaultDRate();

	std::cerr << "[V] drate=" << drate << std::endl;

	std::cerr << "[V] checkdisagreement=" << checkdisagreement << std::endl;


	std::cerr << "[V] d=" << d << std::endl;
	std::cerr << "[V] dthres=" << dthres << std::endl;

	std::string const consfnindex = getFastAIndexFileName(consin);

	std::string const inlasindexname = libmaus2::dazzler::align::OverlapIndexer::getIndexName(inlasfn);
	std::string const indalignerlasindexname = libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(inlasfn);

	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	if (
		! libmaus2::util::GetFileSize::fileExists(inlasindexname)
		||
		libmaus2::util::GetFileSize::isOlder(inlasindexname,inlasfn)
	)
	{
		libmaus2::dazzler::align::OverlapIndexer::constructIndex(inlasfn,&std::cerr);
	}

	if (
		! libmaus2::util::GetFileSize::fileExists(indalignerlasindexname)
		||
		libmaus2::util::GetFileSize::isOlder(indalignerlasindexname,inlasfn)
	)
	{
		libmaus2::dazzler::align::OverlapIndexer::constructIndex(inlasfn,&std::cerr);
	}

	std::string const memlasindexname = "mem:.input.las.bidx";
	std::cerr << "[V] copying " << inlasindexname << " to " << memlasindexname << "...";
	libmaus2::util::GetFileSize::copy(inlasindexname,memlasindexname);
	std::cerr << "done." << std::endl;

	std::string const memdalignerlasindexname = "mem:.input.las.idx";
	std::cerr << "[V] copying " << indalignerlasindexname << " to " << memdalignerlasindexname << "...";
	libmaus2::util::GetFileSize::copy(indalignerlasindexname,memdalignerlasindexname);
	std::cerr << "done." << std::endl;

	int64_t minaread;
	int64_t maxaread;

	#if defined(PHASER_HANDLE_SINGLE)
	minaread = PHASER_HANDLE_SINGLE;
	maxaread = PHASER_HANDLE_SINGLE;
	#else

	minaread = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(inlasfn);
	maxaread = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(inlasfn);

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

		int64_t const toparead = maxaread + 1;
		int64_t const readspan = (toparead > minaread) ? (toparead-minaread) : 0;

		if ( readspan && ! Idiv )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] denominator of J argument cannot be zero" << std::endl;
			lme.finish();
			throw lme;
		}

		if ( toparead > minaread )
		{
			int64_t const partsize = Idiv ? (readspan + Idiv - 1)/Idiv : 0;

			int64_t const ilow = std::min(minaread + Icnt * partsize,toparead);
			int64_t const ihigh = std::min(ilow+partsize,toparead);

			if ( ihigh > ilow )
			{
				minaread = ilow;
				maxaread = ihigh-1;
			}
			else
			{
				minaread = 0;
				maxaread = -1;
			}
		}
	}
	else if ( arg.uniqueArgPresent("I") )
	{
		std::string const Is = arg["I"];
		std::istringstream istr(Is);
		int64_t Iminaread;
		istr >> Iminaread;

		if ( ! istr )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		int const c = istr.get();

		if ( ! istr || c == std::istream::traits_type::eof() || c != ',' )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		int64_t Imaxaread;
		istr >> Imaxaread;

		if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		minaread = std::max(Iminaread,minaread);
		maxaread = std::min(Imaxaread,maxaread);
	}
	#endif

	int64_t const toparead = maxaread >= 0 ? maxaread + 1 : maxaread;

	std::cerr << "[V] minaread=" << minaread << " toparead=" << toparead << std::endl;

	uint64_t const maxinput = arg.uniqueArgPresent("D") ? arg.getUnsignedNumericArg<uint64_t>("D") : getDefaultMaxInput();

	// number of threads
	uint64_t const numthreads =
		#if defined(PHASE_SERIAL)
		1
		#else
		arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads()
		#endif
		;

	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(numthreads);

	for ( uint64_t i = 0; i < Adalindex.size(); ++i )
	{
		libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
			new libmaus2::dazzler::align::DalignerIndexDecoder(inlasfn,memdalignerlasindexname)
		);
		Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
	}

	libmaus2::autoarray::AutoArray < libmaus2::lcs::Aligner::unique_ptr_type > ANP(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::lcs::Aligner::unique_ptr_type tptr(AlignerBase::getAligner());
		ANP[i] = UNIQUE_PTR_MOVE(tptr);
	}

	std::cerr << "[V] copying " << dbname << " to memory...";
	std::vector<std::string> tracklist;
	// tracklist.push_back("exqual");
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbname,&tracklist)
	);
	std::cerr << "done." << std::endl;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
	PDB->computeTrimVector();

	double const avgreadlen = PDB->getAverageReadLength();
	std::cerr << "[V] avgreadlen=" << avgreadlen << std::endl;

	uint64_t const a = (avgreadlen+d-1) / d;
	uint64_t const dprime = (avgreadlen - phasethres) / a;

	std::pair<uint64_t,double> const P = poissonlimit(dprime,0.01);

	std::cerr << "[V] P.first=" << P.first << std::endl;

	std::cerr << "[V] d'=" << dprime << std::endl;
	uint64_t const kd = libmaus2::math::Binom::binomialCoefficientInteger(2,P.first) * libmaus2::math::Binom::binomialCoefficientInteger(2,P.first);
	std::cerr << "[V] (d' 2)(d' 2)=" << kd << std::endl;

	uint64_t const dcthres = arg.uniqueArgPresent("dcthres") ? arg.getUnsignedNumericArg<uint64_t>("dcthres") : kd;

	std::cerr << "[V] using dcthres=" << dcthres << std::endl;

	libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type db2ptr;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB2;

	if ( db2name != dbname )
	{
		std::cerr << "[V] copying " << db2name << " to memory...";
		libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type tdbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(db2name,"mem:db2prefix"));
		std::cerr << "done." << std::endl;
		db2ptr = UNIQUE_PTR_MOVE(tdbptr);

		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type TDB(new libmaus2::dazzler::db::DatabaseFile(db2ptr->fn));
		TDB->computeTrimVector();

		PDB2 = UNIQUE_PTR_MOVE(TDB);
	}

	libmaus2::dazzler::db::DatabaseFile const & DB = *PDB;
	libmaus2::dazzler::db::DatabaseFile const & DB2 = PDB2 ? (*PDB2) : (*PDB);

	// load read/read alignments
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(inlasfn);

	libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> traceFreeList;
	libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> readDataFreeList;
	ReadDecoderAllocator RDA(&DB);
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo>::unique_ptr_type PreadDecoderFreeList(
		new libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo>(RDA)
	);
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList = *PreadDecoderFreeList;

	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo>::unique_ptr_type PreadDecoderFreeList2;
	ReadDecoderAllocator RDA2(&DB2);

	if ( &DB != &DB2 )
	{
		libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo>::unique_ptr_type TreadDecoderFreeList(
			new libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo>(RDA2)
		);
		PreadDecoderFreeList2 = UNIQUE_PTR_MOVE(TreadDecoderFreeList);
	}

	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList2 =
		PreadDecoderFreeList2 ? (*PreadDecoderFreeList2) : (*PreadDecoderFreeList);

	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > > AAOVL(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
		AAOVL[i] = libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap >(maxinput);
	libmaus2::autoarray::AutoArray< libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >::unique_ptr_type > AAH(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >::unique_ptr_type tptr(
			new libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >(maxinput)
		);
		AAH[i] = UNIQUE_PTR_MOVE(tptr);
	}

	#if defined(PHASE_DEBUG)
	uint64_t const bamindexmod = 1024;
	libmaus2::bambam::BamNumericalIndexGenerator::indexFileCheck(bamfn,bamindexmod,1 /* numthreads */);
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::bambam::BamAlignment,BamAlignmentAllocator,BamAlignmentTypeInfo> bamAlignmentFreeList;
	libmaus2::autoarray::AutoArray < BamAlignmentContainer::unique_ptr_type > ABAC(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		BamAlignmentContainer::unique_ptr_type tptr(new BamAlignmentContainer(bamfn,bamAlignmentFreeList));
		ABAC[i] = UNIQUE_PTR_MOVE(tptr);
	}
	#endif

	struct LogEntry
	{
		uint64_t id;
		std::string log;

		LogEntry() {}
		LogEntry(uint64_t const rid, std::string const & rlog)
		: id(rid), log(rlog) {}

		bool operator<(LogEntry const & L) const
		{
			return id > L.id;
		}
	};

	std::priority_queue<LogEntry> logQ;
	std::priority_queue<LogEntry> outQ;
	libmaus2::parallel::PosixMutex logQLock;
	uint64_t volatile logQprintnext = minaread;

	libmaus2::parallel::SynchronousCounter<uint64_t> wellcounter(0);

	#if defined(PHASER_HANDLE_SINGLE)
	libmaus2::dazzler::align::AlignmentWriter AW("debug_cons.las",tspace,false);
	#endif

	libmaus2::dazzler::align::AlignmentWriterArray AWA(tmpfilebase + "_phase_array",numthreads,tspace);
	libmaus2::dazzler::align::AlignmentWriterArray AWK(tmpfilebase + "_phase_drop_array",numthreads,tspace);

	#if defined(_OPENMP)
		#if ! defined(PHASE_SERIAL)
		#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
		#endif
	#endif
	for ( int64_t z = minaread; z < toparead; ++z )
	{
		try
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			std::vector < libmaus2::fastx::FastAReader::pattern_type > Vpat;
			if ( z + 1 < static_cast<int64_t>(VCNT.size()) )
			{
				// std::cerr << "z=" << z << " " << VCNT[z] << " " << VCNT[z+1] << std::endl;

				uint64_t const low = VCNT[z];
				uint64_t const high = VCNT[z+1];
				ReadAccessor RA(consin, consfnindex, fastaindexmod);

				for ( uint64_t i = low; i < high; ++i )
				{
					Vpat.push_back(RA[i]);

					assert (
						static_cast<int64_t>(getId(
							Vpat.back().sid,
							Vpat.back()
						))
						==
						z
					);

					// std::cerr << "got " << Vpat.back().sid << " for " << z << std::endl;
				}
			}

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & RH = *AAH[tid];
			RH.clear();
			libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

			libmaus2::dazzler::align::Overlap OVL;
			libmaus2::dazzler::align::AlignmentFileDecoder::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileAt(inlasfn,z,z+1,*Adalindex[tid]));
			// libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
			uint64_t f = 0;
			while ( pdec->getNextOverlap(OVL) )
			{
				// higher error rate -> higher score
				uint64_t const score =
					static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

				if ( RH.full() )
				{
					if ( score > RH.top().first )
					{
						continue;
					}
					else
					{
						uint64_t const p = RH.top().second;
						RH.pop();
						RO[p] = OVL;
						RH.push(std::pair<uint64_t,uint64_t>(score,p));
					}
				}
				else
				{
					uint64_t const p = f++;
					std::pair<uint64_t,uint64_t> P(score,p);
					RH.push(P);
					RO[p] = OVL;
				}

				// VOVL.push_back(OVL);
			}

			std::sort(RO.begin(),RO.begin()+f,libmaus2::dazzler::align::OverlapFullComparator());

			for ( uint64_t i = 1; i < f; ++i )
			{
				assert ( RO[i].aread == RO[0].aread );
				assert ( RO[i-1].bread <= RO[i].bread );
			}

			// std::sort(RO.begin(),RO.begin()+f,OverlapPosComparator());

			#if defined(PHASER_HANDLE_SINGLE)
			if ( f && RO[0].aread == PHASER_HANDLE_SINGLE )
			{
				for ( uint64_t i = 0; i < f; ++i )
				{
					// std::cerr << RO[i] << std::endl;
					AW.put(RO[i]);
				}
			}
			#endif

			std::ostringstream logostr;
			std::ostringstream outstr;

			#if defined(PHASER_HANDLE_SINGLE)
			for ( uint64_t i = 0; i < numchains; ++i )
			{
				uint64_t const len = chainlengths[i+1]-chainlengths[i];
				std::cerr << "chain " << i << " length " << len << std::endl;
				for ( uint64_t j = 0; j < len; ++j )
					std::cerr << "\t" << RO[chains[chainlengths[i]+j]].getHeader() << std::endl;
			}
			#endif

			//libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
			// DLOS << "[S] handling aread " << Vovl[low].aread << std::endl;
			if ( f )
				logostr << "[S] handling aread " << RO[0].aread << std::endl;

			// std::pair<unsigned char const *, unsigned char const *> const Q = inqual.getQualityForRead(Vovl[low].aread);

			std::vector < ReadFragment > fragments;
			for ( uint64_t i = 0; i < Vpat.size(); ++i )
			{
				std::string sid = Vpat[i].sid;

				assert ( sid.find(" A=[") != std::string::npos );
				sid = sid.substr(sid.find(" A=[") + strlen(" A=["));

				std::istringstream istr(sid);
				uint64_t first = 0;
				istr >> first;
				assert ( istr && istr.peek() == ',' );
				istr.get();
				uint64_t last = 0;
				istr >> last;
				assert ( istr && istr.peek() == ']' );

				std::string const SCO = Vpat[i].spattern;

				ReadFragment RF(first,last+1,SCO);

				fragments.push_back(RF);
			}

			if ( f )
			{
				ChainSet CS(RO.begin(),f);

				try
				{
					#if defined(PHASER_HANDLE_SINGLE)
					if ( RO[0].aread == PHASER_HANDLE_SINGLE )
					#endif
					{
						DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
						DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);
						libmaus2::lcs::Aligner & al = *(ANP[tid]);

						std::ostream & errchannel = (numthreads > 1) ? logostr : std::cerr;

						phaseChains(
							errchannel,
							fragments,
							CS.chainlengths,
							CS.chains,
							CS.numchains,
							RO.begin(),RO.begin()+f,
							RDC,
							RDC2,
							traceFreeList,
							al,
							tspace,
							AWA[tid],
							AWK[tid],
							d,
							dthres,
							crate,
							drate,
							phasethres,
							checkdisagreement,
							dcthres
							#if defined(PHASE_DEBUG)
								,
								*(ABAC[tid]),
								DB,
								DB2
							#endif
						);
					}
				}
				catch(std::exception const & ex)
				{
					try
					{
						libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
						DLOS << ex.what() << std::endl;
					}
					catch(...)
					{
						libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << "[E] exception handling failed" << std::endl;
					}
				}
			}

			#if !defined(PHASE_SERIAL)
			{
				libmaus2::parallel::ScopePosixMutex slock(logQLock);
				logQ.push(LogEntry(z,logostr.str()));
				outQ.push(LogEntry(z,outstr.str()));
			}
			#endif

			bool logrunning = true;

			while ( logrunning )
			{
				logrunning = false;

				std::vector < LogEntry > VLE;
				std::vector < LogEntry > VOE;

				{
					libmaus2::parallel::ScopePosixMutex slock(logQLock);
					uint64_t locallogQprintnext = logQprintnext;
					while (
						(!logQ.empty())
						&&
						logQ.top().id == locallogQprintnext
					)
					{
						assert (! outQ.empty() );
						VLE.push_back(logQ.top());
						VOE.push_back(outQ.top());
						logQ.pop();
						outQ.pop();
						++locallogQprintnext;
					}
				}

				if ( VLE.size() )
				{
					logrunning = true;

					for ( uint64_t i = 0; i < VLE.size(); ++i )
					{
						libmaus2::parallel::ScopePosixSpinLock spsl(libmaus2::aio::StreamLock::cerrlock);
						std::cerr << VLE[i].log;
					}
					for ( uint64_t i = 0; i < VOE.size(); ++i )
					{
						libmaus2::parallel::ScopePosixSpinLock spsl(libmaus2::aio::StreamLock::coutlock);
						std::cout << VOE[i].log;
					}

					libmaus2::parallel::ScopePosixMutex slock(logQLock);
					logQprintnext += VLE.size();
				}
			}
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
		}
	}

	AWA.merge(
		outlasfn,
		tmpfilebase + "_phase_merge_tmp"
	);

	std::string const outdroplasfn =
		libmaus2::util::OutputFileNameTools::clipOff(outlasfn,".las") + "_drop.las";

	AWK.merge(
		outdroplasfn,
		tmpfilebase + "_phase_merge_drop_tmp"
	);

	return EXIT_SUCCESS;
}

/*
 parameters:

 -t : default number of logical cores, threads
 */

static std::string helpMessage(libmaus2::util::ArgParser const & arg)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("d", formatRHS("sequencing depth",getDefaultSequencingDepth())));
	optionMap . push_back ( std::pair < std::string, std::string >("D", formatRHS("max align",getDefaultMaxInput())));
	optionMap . push_back ( std::pair < std::string, std::string >("p", formatRHS("phase threshold",getDefaultPhaseThres())));
	optionMap . push_back ( std::pair < std::string, std::string >("I", formatRHS("read interval",std::string("0,") + libmaus2::util::NumberSerialisation::formatNumber(std::numeric_limits<uint64_t>::max(),0))));
	optionMap . push_back ( std::pair < std::string, std::string >("J", formatRHS("reads part",std::string("0,1"))));
	optionMap . push_back ( std::pair < std::string, std::string >("E", formatRHS("error profile file name",std::string("in.las.eprof"))));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("drate", formatRHS("difference rate",getDefaultDRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("dcthres", formatRHS("splitting tuple count threshold","computed using d parameter")));

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
		#if defined(PHASE_DEBUG)
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 5 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " -d<depth> [options] out.las cons.fasta in.las in1.dam in2.dam in.bam\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::timing::RealTimeClock rtc;
			rtc.start();

			int r = EXIT_FAILURE;

			r = phaser(arg,arg[0] /* outlas */,arg[1] /* cons */,arg[2] /* in las */, arg[3] /* db1 */, arg[4] /* db2 */, arg[5] /* bam */);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			return r;
		}
		#else
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 4 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " -d<depth> [options] out.las cons.fasta in.las in.dam\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::timing::RealTimeClock rtc;
			rtc.start();

			int r = EXIT_FAILURE;

			if ( arg.size() == 4 )
				r = phaser(arg,arg[0],arg[1],arg[2],arg[3],arg[3]);
			else
				r = phaser(arg,arg[0],arg[1],arg[2],arg[3],arg[4]);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			return r;
		}
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[E] " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
