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

#include <libmaus2/dazzler/align/TrueOverlap.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferFragment.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NPL.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/bambam/BamRangeDecoder.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/util/TempFileNameGenerator.hpp>

struct BamLoader
{
	std::string const bamfn;
	std::string const indexfn;
	libmaus2::bambam::BamNumericalIndexDecoder indexdec;

	BamLoader(std::string const & rbamfn)
	: bamfn(rbamfn), indexfn(libmaus2::bambam::BamNumericalIndexBase::getIndexName(bamfn)), indexdec(indexfn) {}

	void loadAlignmentAt(uint64_t const i, libmaus2::bambam::BamAlignment & algn)
	{
		libmaus2::lz::BgzfInflateFile::unique_ptr_type tptr(indexdec.getStreamAt(bamfn,i));
		libmaus2::bambam::BamAlignmentDecoder::readAlignmentGz(*tptr,algn);
	}

	uint64_t getSplit(uint64_t i)
	{
		while ( i  )
		{
			libmaus2::bambam::BamAlignment A0, A1;
			loadAlignmentAt(i,A0);
			loadAlignmentAt(i-1,A1);

			if ( A0.getRefID() != A1.getRefID() || A0.getPos() != A1.getPos() )
				return i;
			else
				--i;
		}

		return i;
	}

	std::pair<int32_t,int32_t> getRefPos(uint64_t const i)
	{
		libmaus2::bambam::BamAlignment A;
		loadAlignmentAt(i,A);
		return std::pair<int32_t,int32_t>(A.getRefID(),A.getPos());
	}

	std::vector<uint64_t> getSplitVector(uint64_t const n, uint64_t const numthreads = 1)
	{
		std::vector<uint64_t> V;
		uint64_t const s = indexdec.getAlignmentCount();

		if ( ! s )
			return V;

		uint64_t const b = (s + n - 1)/n;
		uint64_t const v = (s + b - 1)/b;

		V.resize(v);

		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t i = 0; i < v; ++i )
			V[i] = getSplit(i*b);

		return V;
	}
};

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
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

struct AutoArrayTypeInfo
{
	typedef libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement> element_type;
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

struct AutoArrayAllocator
{
	typedef libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement> element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct RefPosComp
{
	bool operator()(libmaus2::bambam::PileVectorElement const & A, libmaus2::bambam::PileVectorElement const & B) const
	{
		if ( A.refpos != B.refpos )
			return A.refpos < B.refpos;
		else
			return A.predif < B.predif;
	}
};

int64_t getId(std::string const & readname)
{
	std::deque<std::string> V = libmaus2::util::stringFunctions::tokenize(readname,std::string("/"));

	if ( 1 < V.size() )
	{
		std::istringstream istr(V[1]);
		uint64_t id;
		istr >> id;

		if ( istr && istr.peek() == std::istream::traits_type::eof() )
			return id;
	}

	return -1;
}

void complete(
	libmaus2::bambam::BamAlignment const & algnA,
	libmaus2::bambam::BamAlignment const & algnB,
	libmaus2::lcs::AlignmentTraceContainer & ATC,
	int64_t & abpos,
	int64_t & aepos,
	int64_t & bbpos,
	int64_t & bepos
)
{
	std::string a = algnA.getRead();
	std::string b = algnB.getRead();

	bool reva = algnA.isReverse();
	bool revb = algnB.isReverse();
	bool const swap = reva;

	int64_t afrontclip = algnA.getFrontSoftClipping();
	int64_t abackclip = algnA.getBackSoftClipping();
	int64_t bfrontclip = algnB.getFrontSoftClipping();
	int64_t bbackclip = algnB.getBackSoftClipping();

	if ( swap )
	{
		reva = ! reva;
		revb = ! revb;
		std::reverse(ATC.ta,ATC.te);

		int64_t const nabpos = algnA.getLseq() - aepos;
		int64_t const naepos = algnA.getLseq() - abpos;
		int64_t const nbbpos = algnB.getLseq() - bepos;
		int64_t const nbepos = algnB.getLseq() - bbpos;

		assert ( nabpos >= 0 );
		assert ( naepos >= nabpos );
		assert ( naepos <= algnA.getLseq() );

		assert ( nbbpos >= 0 );
		assert ( nbepos >= nbbpos );
		assert ( nbepos <= algnB.getLseq() );

		abpos = nabpos;
		aepos = naepos;
		bbpos = nbbpos;
		bepos = nbepos;

		a = libmaus2::fastx::reverseComplementUnmapped(a);
		b = libmaus2::fastx::reverseComplementUnmapped(b);

		std::swap(afrontclip,abackclip);
		std::swap(bfrontclip,bbackclip);
	}

	std::string const origb = algnB.isReverse() ? algnB.getReadRC() : algnB.getRead();

	assert (  origb == (revb ? libmaus2::fastx::reverseComplementUnmapped(b) : b) );

	#if 0
	libmaus2::lcs::AlignmentPrint::printAlignmentLines(
		std::cerr,
		a.begin() + abpos,
		aepos-abpos,
		b.begin() + bbpos,
		bepos-bbpos,
		80,
		ATC.ta,
		ATC.te
	);
	#endif

	assert (
		libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
			ATC.ta,
			ATC.te,
			a.begin() + abpos,
			b.begin() + bbpos
		)
	);

	if ( abpos > afrontclip && bbpos > bfrontclip )
	{
		std::string fa = a.substr(afrontclip,abpos-afrontclip);
		std::string fb = b.substr(bfrontclip,bbpos-bfrontclip);

		std::reverse(fa.begin(),fa.end());
		std::reverse(fb.begin(),fb.end());

		libmaus2::lcs::NPL npl;
		npl.np(fa.begin(),fa.end(),fb.begin(),fb.end());

		std::reverse(npl.ta,npl.te);

		std::pair<int64_t,int64_t> SL = npl.getTraceContainer().getStringLengthUsed();
		assert ( SL.first <= abpos );
		assert ( SL.second <= bbpos );

		// assert ( SL.first == abpos || SL.second == bbpos );

		#if 0
		if ( SL.first || SL.second )
			std::cerr << "[V] complete front" << " abpos=" << abpos << " bbpos=" << bbpos << " SL.first=" << SL.first << " SL.second=" << SL.second << std::endl;
		#endif

		libmaus2::lcs::AlignmentTraceContainer TATC;
		TATC.push(npl);
		TATC.push(ATC);

		ATC.reset();
		ATC.push(TATC);

		abpos -= SL.first;
		bbpos -= SL.second;

		assert ( abpos == afrontclip || bbpos == bfrontclip );

		assert (
			libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
				ATC.ta,
				ATC.te,
				a.begin() + abpos,
				b.begin() + bbpos
			)
		);
	}
	if ( algnA.getLseq() - aepos > abackclip && algnB.getLseq()-bepos > bbackclip )
	{
		uint64_t const availa = algnA.getLseq() - abackclip - aepos;
		uint64_t const availb = algnB.getLseq() - bbackclip - bepos;

		std::string fa = a.substr(aepos,availa);
		std::string fb = b.substr(bepos,availb);

		libmaus2::lcs::NPL npl;
		npl.np(fa.begin(),fa.end(),fb.begin(),fb.end());

		std::pair<int64_t,int64_t> SL = npl.getTraceContainer().getStringLengthUsed();
		assert ( SL.first <= algnA.getLseq() - aepos );
		assert ( SL.second <= algnB.getLseq() - bepos );

		#if 0
		if ( SL.first || SL.second )
			std::cerr << "[V] complete back" << " aepos=" << aepos << " bepos=" << bepos << " SL.first=" << SL.first << " SL.second=" << SL.second << std::endl;
		#endif

		ATC.push(npl);

		aepos += SL.first;
		bepos += SL.second;

		// assert ( aepos == algnA.getLseq() || bepos == algnB.getLseq() );
		assert ( aepos == algnA.getLseq()-abackclip || bepos == algnB.getLseq()-bbackclip );

		assert (
			libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
				ATC.ta,
				ATC.te,
				a.begin() + abpos,
				b.begin() + bbpos
			)
		);
	}

	if ( swap )
	{
		std::reverse(ATC.ta,ATC.te);

		int64_t const nabpos = algnA.getLseq() - aepos;
		int64_t const naepos = algnA.getLseq() - abpos;
		int64_t const nbbpos = algnB.getLseq() - bepos;
		int64_t const nbepos = algnB.getLseq() - bbpos;

		assert ( nabpos >= 0 );
		assert ( naepos >= nabpos );
		assert ( naepos <= algnA.getLseq() );

		assert ( nbbpos >= 0 );
		assert ( nbepos >= nbbpos );
		assert ( nbepos <= algnB.getLseq() );

		abpos = nabpos;
		aepos = naepos;
		bbpos = nbbpos;
		bepos = nbepos;
	}

	//std::cerr << "[V] completed" << std::endl;
}


void fillOverlap(
	int64_t const tspace,
	libmaus2::bambam::BamAlignment const & algnA,
	libmaus2::bambam::BamAlignment const & algnB,
	libmaus2::lcs::AlignmentTraceContainer const & ATC,
	int64_t abpos,
	int64_t aepos,
	int64_t bbpos,
	int64_t bepos,
	libmaus2::dazzler::align::Overlap & OVL,
	uint64_t const ida,
	uint64_t const idb
)
{
	#if 0
	std::cerr << "ida=" << ida << " idb=" << idb << std::endl;
	#endif

	std::string a = algnA.getRead();
	std::string b = algnB.getRead();

	bool reva = algnA.isReverse();
	bool revb = algnB.isReverse();
	bool const swap = reva;

	if ( swap )
	{
		reva = ! reva;
		revb = ! revb;
		std::reverse(ATC.ta,ATC.te);

		int64_t const nabpos = algnA.getLseq() - aepos;
		int64_t const naepos = algnA.getLseq() - abpos;
		int64_t const nbbpos = algnB.getLseq() - bepos;
		int64_t const nbepos = algnB.getLseq() - bbpos;

		assert ( nabpos >= 0 );
		assert ( naepos >= nabpos );
		assert ( naepos <= algnA.getLseq() );

		assert ( nbbpos >= 0 );
		assert ( nbepos >= nbbpos );
		assert ( nbepos <= algnB.getLseq() );

		abpos = nabpos;
		aepos = naepos;
		bbpos = nbbpos;
		bepos = nbepos;

		a = libmaus2::fastx::reverseComplementUnmapped(a);
		b = libmaus2::fastx::reverseComplementUnmapped(b);
	}

	std::string const origb = algnB.isReverse() ? algnB.getReadRC() : algnB.getRead();

	assert (  origb == (revb ? libmaus2::fastx::reverseComplementUnmapped(b) : b) );

	#if 0
	libmaus2::lcs::AlignmentPrint::printAlignmentLines(
		std::cerr,
		a.begin() + abpos,
		aepos-abpos,
		b.begin() + bbpos,
		bepos-bbpos,
		80,
		ATC.ta,
		ATC.te
	);
	#endif

	assert (
		libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
			ATC.ta,
			ATC.te,
			a.begin() + abpos,
			b.begin() + bbpos
		)
	);

	#if 0
	uint64_t z = 0;
	while ( z < pa )
	{
		uint64_t const l = std::min(pa-z,static_cast<uint64_t>(80));

		std::cerr << std::string(PA.begin()+z,PA.begin()+z+l) << std::endl;
		std::cerr << std::string(PB.begin()+z,PB.begin()+z+l) << std::endl;

		z += l;
	}
	#endif

	OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
		revb ? libmaus2::dazzler::align::Overlap::getInverseFlag() : 0,
		ida,
		idb,
		abpos,
		aepos,
		bbpos,
		bepos,
		tspace,
		ATC
	);

	libmaus2::lcs::AlignmentTraceContainer DATC;
	libmaus2::lcs::NP NP;

	uint8_t const * aptr = reinterpret_cast<uint8_t const *>(a.c_str());
	std::string const ovlb = OVL.isInverse() ? libmaus2::fastx::reverseComplementUnmapped(origb) : origb;
	uint8_t const * bptr = reinterpret_cast<uint8_t const *>(b.c_str());

	OVL.computeTrace(aptr,bptr,tspace,DATC,NP);

	#if 0
	std::cerr << OVL << std::endl;
	#endif

	assert (
		libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
			DATC.ta,
			DATC.te,
			aptr + OVL.path.abpos,
			bptr + OVL.path.bbpos
		)
	);

	#if 0
	libmaus2::lcs::AlignmentPrint::printAlignmentLines(
		std::cerr,
		aptr + OVL.path.abpos,
		OVL.path.aepos-OVL.path.abpos,
		bptr + OVL.path.bbpos,
		OVL.path.bepos-OVL.path.bbpos,
		80,
		DATC.ta,
		DATC.te
	);

	std::cerr << DATC.getAlignmentStatistics() << std::endl;
	#endif


	#if 0
	std::cerr << OVL << std::endl;

	std::cerr << std::string(80,'-') << std::endl;
	#endif

	if ( swap )
		std::reverse(ATC.ta,ATC.te);
}

bool computeCommonPileVector(
	int64_t const tspace,
	libmaus2::bambam::BamAlignment const & algnA,
	libmaus2::bambam::BamAlignment const & algnB,
	libmaus2::bambam::PileVectorElement const * const A,
	uint64_t const oa,
	libmaus2::bambam::PileVectorElement const * const B,
	uint64_t const ob,
	libmaus2::autoarray::AutoArray<char> & PA,
	libmaus2::autoarray::AutoArray<char> & PB,
	libmaus2::lcs::AlignmentTraceContainer & ATC,
	libmaus2::dazzler::align::Overlap & OVLf,
	libmaus2::dazzler::align::Overlap & OVLr
)
{
	ATC.reset();

	int64_t const ida = getId(algnA.getName());
	int64_t const idb = getId(algnB.getName());

	#if 0
	std::cerr << algnA.getName() << " " << ida << std::endl;
	std::cerr << algnB.getName() << " " << idb << std::endl;
	#endif

	assert ( ida >= 0 && idb >= 0 );

	uint64_t la = 0, ha = oa;
	uint64_t lb = 0, hb = ob;

	RefPosComp comp;
	for ( uint64_t i = la+1; i < ha; ++i )
		assert ( comp(A[i-1],A[i]) );
	for ( uint64_t i = lb+1; i < hb; ++i )
		assert ( comp(B[i-1],B[i]) );

	while ( la < ha && lb < hb && (comp(A[la],B[lb]) || comp(B[lb],A[la]) || A[la].sym == '-' || B[lb].sym == '-' ) )
	{
		if ( comp(A[la],B[lb]) )
			++la;
		else if ( comp(B[lb],A[la]) )
			++lb;
		else
		{
			++la;
			++lb;
		}
	}

	while ( la < ha && lb < hb &&
		(
			comp(A[ha-1],B[hb-1])
			||
			comp(B[hb-1],A[ha-1])
			||
			A[ha-1].sym == '-'
			||
			B[hb-1].sym == '-'
		)
	)
	{
		if ( comp(A[ha-1],B[hb-1]) )
			--hb;
		else if ( comp(B[hb-1],A[ha-1]) )
			--ha;
		else
		{
			--ha;
			--hb;
		}
	}

	assert ( la <= ha );
	assert ( lb <= hb );
	assert ( la == ha || lb == hb || (A[la].refpos == B[lb].refpos && A[la].predif == B[lb].predif) );
	assert ( la == ha || lb == hb || (A[ha-1].refpos == B[hb-1].refpos && A[ha-1].predif == B[hb-1].predif) );

	if ( la != ha && lb != hb )
	{
		#if 0
		std::cerr << A[la] << std::endl;
		std::cerr << B[lb] << std::endl;
		std::cerr << A[ha-1] << std::endl;
		std::cerr << B[hb-1] << std::endl;
		#endif

		uint64_t pa = 0, pb = 0;

		uint64_t ia = la, ib = lb;

		while ( ia < ha && ib < hb )
		{
			if ( comp(A[ia],B[ib]) )
			{
				if ( A[ia].sym != '-' )
				{
					PA.push(pa,A[ia].sym);
					PB.push(pb,'-');
				}
				++ia;
			}
			else if ( comp(B[ib],A[ia]) )
			{
				if ( B[ib].sym != '-' )
				{
					PA.push(pa,'-');
					PB.push(pb,B[ib].sym);
				}
				++ib;
			}
			else
			{
				if ( A[ia].sym != '-' || B[ib].sym != '-' )
				{
					PA.push(pa,A[ia].sym);
					PB.push(pb,B[ib].sym);
				}
				++ia;
				++ib;
			}
		}

		assert ( ia == ha && ib == hb );

		assert ( pa == pb );

		uint64_t na = 0, nb = 0;
		for ( uint64_t i = 0; i < pa; ++i )
		{
			if ( PA[i] != '-' )
				++na;
			if ( PB[i] != '-' )
				++nb;
		}

		#if 0
		std::cerr << na << " " << A[ha-1].readpos - A[la].readpos << std::endl;
		std::cerr << nb << " " << B[hb-1].readpos - B[lb].readpos << std::endl;
		#endif

		assert ( static_cast<int64_t>(na) == static_cast<int64_t>(A[ha-1].readpos - A[la].readpos + 1) );
		assert ( static_cast<int64_t>(nb) == static_cast<int64_t>(B[hb-1].readpos - B[lb].readpos + 1) );

		int64_t abpos = A[la].readpos;
		int64_t aepos = A[ha-1].readpos+1;
		int64_t bbpos = B[lb].readpos;
		int64_t bepos = B[hb-1].readpos+1;

		assert ( aepos >= abpos );
		assert ( bepos >= bbpos );
		assert ( abpos >= 0 );
		assert ( bbpos >= 0 );

		// std::cerr << "abpos=" << abpos << " aepos=" << aepos << std::endl;
		// std::cerr << "bbpos=" << bbpos << " bepos=" << bepos << std::endl;

		assert ( aepos <= algnA.getLseq() );
		assert ( bepos <= algnB.getLseq() );

		if ( ATC.capacity() < pa )
			ATC.resize(pa);
		ATC.reset();

		ATC.ta -= pa;

		for ( uint64_t i = 0; i < pa; ++i )
		{
			if ( PA[i] == PB[i] )
			{
				assert ( PA[i] != '-' );
				ATC.ta[i] = libmaus2::lcs::BaseConstants::STEP_MATCH;
			}
			else
			{
				assert ( PA[i] != PB[i] );

				if ( PA[i] == '-' )
					ATC.ta[i] = libmaus2::lcs::BaseConstants::STEP_INS;
				else if ( PB[i] == '-' )
					ATC.ta[i] = libmaus2::lcs::BaseConstants::STEP_DEL;
				else
					ATC.ta[i] = libmaus2::lcs::BaseConstants::STEP_MISMATCH;
			}
		}

		complete(algnA,algnB,ATC,abpos,aepos,bbpos,bepos);

		fillOverlap(tspace,algnA,algnB,ATC,abpos,aepos,bbpos,bepos,OVLf,ida,idb);
		ATC.swapRoles();
		fillOverlap(tspace,algnB,algnA,ATC,bbpos,bepos,abpos,aepos,OVLr,idb,ida);

		return true;
	}
	else
	{
		// std::cerr << "nothing found" << std::endl;
		// assert ( false );
		return false;
	}

	// return 0;
}

// #define SEQDEBUG

static uint64_t getDefaultTSpace()
{
	return 100;
}

int generateperfectpiles(libmaus2::util::ArgParser const & arg)
{
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::bambam::BamAlignment,BamAlignmentAllocator,BamAlignmentTypeInfo> bamFL;
	libmaus2::parallel::LockedGrowingFreeList<libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>,AutoArrayAllocator,AutoArrayTypeInfo> AAFL;

	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	int64_t const tspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : getDefaultTSpace();

	std::cerr << "[V] using tspace=" << tspace << std::endl;

	std::string const lasfn = arg[0];
	std::string const bamfn = arg[1];

	libmaus2::bambam::BamHeader::unique_ptr_type Pheader;

	{
		libmaus2::bambam::BamDecoder dec(bamfn);
		libmaus2::bambam::BamHeader::unique_ptr_type Theader(dec.getHeader().uclone());
		Pheader = UNIQUE_PTR_MOVE(Theader);
	}

	#if defined(SEQDEBUG)
	std::string const dbfn = arg[2];
	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	DB.computeTrimVector();
	#endif

	uint64_t const genindexmod = 1024;
	libmaus2::bambam::BamNumericalIndexGenerator::indexFileCheck(bamfn,genindexmod,1 /* numthreads */);

	BamLoader BL(bamfn);
	uint64_t const mult = 64;
	std::vector<uint64_t> const splitV = BL.getSplitVector(mult*numthreads,numthreads);

	#if 0
	for ( uint64_t i = 0; i < splitV.size(); ++i )
		std::cerr << "split " << splitV[i] << std::endl;
	#endif

	std::string const tmplasfnbase = tmpfilebase + "_las_tmp";
	libmaus2::util::TempFileNameGenerator tmpgen(tmplasfnbase,3);
	std::vector<std::string> Vtmplasfn(splitV.size());
	for ( uint64_t i = 0; i < Vtmplasfn.size(); ++i )
	{
		std::ostringstream ostr;
		ostr << tmpgen.getFileName() << "." << std::setw(6) << std::setfill('0') << i;
		Vtmplasfn[i] = ostr.str();
	}

	int volatile failed = 0;
	libmaus2::parallel::PosixSpinLock failedlock;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t is = 0; is < splitV.size(); ++is )
	{
		try
		{
			uint64_t const ilow = splitV[is];
			uint64_t const ihigh = (is+1<splitV.size()) ? splitV[is+1] : BL.indexdec.getAlignmentCount();
			libmaus2::lz::BgzfInflateFile::unique_ptr_type tptr(BL.indexdec.getStreamAt(bamfn,ilow));

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > E(1024);
			int64_t prevrefid = std::numeric_limits<int64_t>::max();
			std::map<uint64_t,libmaus2::bambam::BamAlignment::shared_ptr_type> active;
			typedef libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement> pile_vector_array;
			typedef pile_vector_array::shared_ptr_type pile_vector_array_pointer;
			typedef std::pair < uint64_t, pile_vector_array_pointer > pvm_pair;
			std::map<uint64_t,pvm_pair> PIVM;
			libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigopin;
			libmaus2::autoarray::AutoArray<char> readdata;

			std::pair<int32_t,int32_t> const RP = BL.getRefPos(ilow);
			assert ( RP.first >= 0 );
			assert ( RP.second >= 0 );

			if ( RP.first >= 0 && RP.second >= 0 )
			{
				std::string const refname = Pheader->getRefIDName(RP.first);
				std::ostringstream rangestr;
				rangestr << refname << ":" << (RP.second+1) << "-" << (RP.second+1);
				std::string const range = rangestr.str();
				// std::cerr << "range=" << range << std::endl;
				libmaus2::bambam::BamRangeDecoder BRD(bamfn,range);
				libmaus2::bambam::BamAlignment const & algn = BRD.getAlignment();

				for ( uint64_t id = ihigh; BRD.readAlignment() && algn.getPos() < RP.second; ++id )
				{
					libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>::shared_ptr_type PV = AAFL.get();
					uint64_t const curo = algn.getPileVector(*PV,cigopin,readdata,id);

					// std::cerr << algn.getName() << " "  << algn.getRefID() << "," << algn.getPos() << "," << algn.getReferenceInterval() << std::endl;
					libmaus2::bambam::BamAlignment::shared_ptr_type sptr = bamFL.get();
					sptr->copyFrom(algn);
					active[id] = sptr;
					PIVM[id] = std::pair < uint64_t, libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>::shared_ptr_type >(curo,PV);
					E.pushBump(std::pair<uint64_t,uint64_t>(algn.getReferenceInterval().to,id));
					prevrefid = algn.getRefID();
				}
			}

			libmaus2::bambam::BamAlignment algn;
			std::string const tmplasfn = Vtmplasfn[is];

			libmaus2::util::TempFileRemovalContainer::addTempFile(tmplasfn);
			libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(new libmaus2::dazzler::align::AlignmentWriter(tmplasfn,tspace,false /* idx */));

			libmaus2::autoarray::AutoArray < libmaus2::bambam::PileVectorElement > Vpile;

			libmaus2::autoarray::AutoArray<char> PA;
			libmaus2::autoarray::AutoArray<char> PB;
			libmaus2::lcs::AlignmentTraceContainer ATC;
			libmaus2::dazzler::align::Overlap OVLf;
			libmaus2::dazzler::align::Overlap OVLr;

			uint64_t keep = 0, drop = 0;
			uint32_t previd = 0, prevpos = 0;

			for ( uint64_t id = ilow; id < ihigh; ++id )
			{
				libmaus2::bambam::BamAlignmentDecoder::readAlignmentGz(*tptr,algn);

				if ( (id-ilow) % 1024 == 0 )
				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << "[V] " << is << " " << (id-ilow) << " keep=" << keep << " drop=" << drop << std::endl;
				}

				uint32_t const curid = static_cast<uint32_t>(algn.getRefID());
				uint32_t const curpos = static_cast<uint32_t>(algn.getPos());

				bool const orderok = curid > previd || (curid == previd && curpos >= prevpos);

				if ( ! orderok )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] input file is not sorted by coordinate, aborting" << std::endl;
					lme.finish();
					throw lme;
				}

				previd = curid;
				prevpos = curpos;

				if ( algn.isUnmap() )
					continue;

				while (
					(!E.empty())
					&&
					(
						algn.getRefID() != prevrefid
						||
						algn.getPos() > static_cast<int64_t>(E.top().first)
					)
				)
				{
					std::pair<uint64_t,uint64_t> const P = E.top();
					E.pop();

					uint64_t const remid = P.second;
					assert ( active.find(remid) != active.end() );
					bamFL.put(active.find(remid)->second);
					active.erase(active.find(remid));
					assert ( PIVM.find(remid) != PIVM.end() );
					AAFL.put(PIVM.find(remid)->second.second);
					PIVM.erase(PIVM.find(remid));
				}

				libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>::shared_ptr_type PV = AAFL.get();
				uint64_t const curo = algn.getPileVector(*PV,cigopin,readdata,id);

				#if defined(SEQDEBUG)
				int64_t const rid = getId(algn.getName());

				if ( algn.isReverse() )
				{
					bool const ok = DB[rid] == algn.getReadRC();
					if ( ! ok )
					{
						std::cerr << algn.getName() << " " << rid << std::endl;
						std::cerr << DB[rid] << std::endl;
						std::cerr << algn.getReadRC() << std::endl;
						assert ( ok );
					}
				}
				else
				{
					bool const ok = DB[rid] == algn.getRead();
					if ( ! ok )
					{
						std::cerr << algn.getName() << " " << rid << std::endl;
						std::cerr << DB[rid] << std::endl;
						std::cerr << algn.getRead() << std::endl;
						assert ( ok );
					}
				}
				#endif

				#if 0
				// int64_t const rid = getId(algn.getName());
				if ( rid == 97776 )
				#endif
				{
					for ( std::map<uint64_t,libmaus2::bambam::BamAlignment::shared_ptr_type>::const_iterator ita = active.begin(); ita != active.end(); ++ita )
					{
						libmaus2::bambam::BamAlignment const & other = *(ita->second);
						assert ( PIVM.find(ita->first) != PIVM.end() );
						std::pair < uint64_t, libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>::shared_ptr_type > const & POPV =
							PIVM.find(ita->first)->second;

						if ( computeCommonPileVector(tspace,algn,other,PV->begin(),curo,POPV.second->begin(),POPV.first,PA,PB,ATC,OVLf,OVLr) )
						{
							AW->put(OVLf);
							AW->put(OVLr);
							keep++;
						}
						else
						{
							drop++;
						}
					}
				}

				libmaus2::bambam::BamAlignment::shared_ptr_type sptr = bamFL.get();
				sptr->copyFrom(algn);
				active[id] = sptr;
				PIVM[id] = std::pair < uint64_t, libmaus2::autoarray::AutoArray<libmaus2::bambam::PileVectorElement>::shared_ptr_type >(curo,PV);
				E.pushBump(std::pair<uint64_t,uint64_t>(algn.getReferenceInterval().to,id));
				prevrefid = algn.getRefID();
			}

			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[V] " << is << " keep=" << keep << " drop=" << drop << std::endl;
			}

			while ( (!E.empty()) )
			{
				std::pair<uint64_t,uint64_t> const P = E.top();
				E.pop();

				uint64_t const remid = P.second;
				assert ( active.find(remid) != active.end() );
				bamFL.put(active.find(remid)->second);
				active.erase(active.find(remid));
				assert ( PIVM.find(remid) != PIVM.end() );
				AAFL.put(PIVM.find(remid)->second.second);
				PIVM.erase(PIVM.find(remid));
			}

			AW.reset();
		}
		catch(std::exception const & ex)
		{
			std::cerr << "[E] " << ex.what() << std::endl;
			failedlock.lock();
			failed = 1;
			failedlock.unlock();
		}
	}

	if ( failed )
	{
		std::cerr << "[E] generation failed" << std::endl;
		return EXIT_FAILURE;
	}

	libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::sortAndMerge(
		Vtmplasfn,
		lasfn,
		tmpfilebase,
		libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::getDefaultMergeFanIn(),
		numthreads
	);
	for ( uint64_t i = 0; i < Vtmplasfn.size(); ++i )
		libmaus2::aio::FileRemoval::removeFile(Vtmplasfn[i]);

	return EXIT_SUCCESS;
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
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("prefix for temporary files",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("tspace", formatRHS("trace point spacing",getDefaultTSpace())));
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las in.bam\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return generateperfectpiles(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
