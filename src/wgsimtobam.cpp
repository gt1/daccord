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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/trie/TrieState.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>

#include <libmaus2/fastx/acgtnMap.hpp>
#include <libmaus2/autoarray/AutoArray.hpp>
#include <libmaus2/math/lowbits.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <queue>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/dazzler/align/TracePoint.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

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

static std::string kmerToString(uint64_t const v, uint64_t const k)
{
	std::ostringstream ostr;

	for ( uint64_t i = 0; i < k; ++i )
		ostr.put(libmaus2::fastx::remapChar((v >> (2*(k-i-1)))&3));

	return ostr.str();
}

uint64_t getKmers(std::string const & s, uint64_t const k, libmaus2::autoarray::AutoArray<Kmer> & K)
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

		// std::cerr << kmerToString(v,k) << std::endl;
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

uint64_t simpleprocess(
	std::string const & sa, std::string const & sb,
	uint64_t const k, int64_t const tspace,
	int64_t const aid, int64_t const bid,
	uint64_t const minlen
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

	libmaus2::lcs::NNP nnp;
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
		HeapTodo H = Q.top();
		Q.pop();

		if ( H.hasNext() )
			Q.push(H.getNext());

		int64_t const nd = H.getDiag();

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
				aid==bid
			);
		}
		else
		{
			res = nnp.align(
				ra.begin(),ra.end(),H.ka->p,
				rb.begin(),rb.end(),H.kb_c->p,
				nnptracecontainer,
				aid==bid
			);
		}

		// compute dense trace
		nnptracecontainer.computeTrace(ATC);

		libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = ATC.ta;
		libmaus2::lcs::AlignmentTraceContainer::step_type const * const te = ATC.te;

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

	uint64_t maxlen = 0;
	for ( std::map<uint64_t,libmaus2::dazzler::align::Overlap>::const_iterator ita = OVLseen.begin();
		ita != OVLseen.end(); ++ita )
	{
		libmaus2::dazzler::align::Overlap const & OVL = ita->second;
		uint64_t const len = OVL.path.aepos - OVL.path.abpos;
		maxlen = std::max(len,maxlen);
	}

	return maxlen;
}


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

static int wgsimtobam(libmaus2::util::ArgParser const & arg, libmaus2::util::ArgInfo const & arginfo)
{
	std::string const reffn = arg[0];

	std::map < std::string, uint64_t > M;
	std::vector < std::pair < std::string, uint64_t> > refmeta;
	std::vector < std::string > Vref;

	{
		libmaus2::fastx::FastAReader FAin(reffn);
		libmaus2::fastx::FastAReader::pattern_type pattern;

		for ( uint64_t id = 0 ; FAin.getNextPatternUnlocked(pattern); ++id )
		{
			std::string const sname = pattern.getShortStringId();

			if ( M.find(sname) != M.end() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] sequence names in " << reffn << " are not unique: " << sname << std::endl;
				lme.finish();
				throw lme;
			}

			M [ sname ] = id;

			std::string seq = pattern.spattern;
			for ( uint64_t i = 0; i < seq.size(); ++i )
				seq[i] = toupper(seq[i]);
			Vref.push_back(seq);

			refmeta.push_back(std::pair < std::string, uint64_t>(sname,seq.size()));
		}
	}

	// create SAM header
	std::ostringstream samheaderstr;
	samheaderstr << "@HD\tVN:1.5\tSO:unknown\n";
	for ( uint64_t i = 0; i < refmeta.size(); ++i )
	{
		samheaderstr << "@SQ\tSN:" << refmeta[i].first << "\tLN:" << refmeta[i].second << "\n";
	}
	libmaus2::bambam::BamHeader header(samheaderstr.str());

	std::vector < std::string > vnames;
	std::string prevname;
	bool prevnamevalid = false;
	for ( std::map < std::string, uint64_t >::const_iterator ita = M.begin(); ita != M.end(); ++ita )
	{
		if ( prevnamevalid && ita->first.size() >= prevname.size() &&
			ita->first.substr(0,prevname.size()) == prevname )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] sequence names in " << reffn << " are not prefix free" << std::endl;
			lme.finish();
			throw lme;
		}

		vnames.push_back(ita->first);

		prevname = ita->first;
		prevnamevalid = true;

	}
	::libmaus2::trie::Trie<char> trienofailure;
	trienofailure.insertContainer(vnames);

	// std::cerr << trienofailure << std::endl;

	::libmaus2::trie::LinearHashTrie<char,uint32_t>::unique_ptr_type LHTnofailure(trienofailure.toLinearHashTrie<uint32_t>());

	libmaus2::bambam::BamSeqEncodeTable const seqenc;
	::libmaus2::fastx::UCharBuffer UB;

	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(
		libmaus2::bambam::BamBlockWriterBaseFactory::construct(header, arginfo)
	);

	std::string const readsfn = arg[1];
	libmaus2::fastx::FastAReader FAin(readsfn);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	for ( uint64_t id = 0 ; FAin.getNextPatternUnlocked(pattern); ++id )
	{
		std::string const sid = pattern.getShortStringId();
		int64_t const seqindex = LHTnofailure->searchNoFailure(sid.begin(),sid.end());

		if ( seqindex < 0 )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] sequence id for " << sid << " not found" << std::endl;
			lme.finish();
			throw lme;
		}

		std::string const & refseqname = vnames[seqindex];
		uint64_t refid = M.find(refseqname)->second;

		std::string const restname = sid.substr(refseqname.size());

		std::istringstream istr(restname);
		assert ( istr.peek() == '_' );
		istr.get();

		int64_t from;
		istr >> from;

		assert ( istr.peek() == '_' );
		istr.get();

		int64_t to;
		istr >> to;

		std::string const reffrag = Vref[refid].substr(from-1,to-from+1);
		std::string readfrag = pattern.spattern;
		for ( uint64_t i = 0; i < readfrag.size(); ++i )
			readfrag[i] = toupper(readfrag[i]);
		std::string const readfragrc = libmaus2::fastx::reverseComplementUnmapped(readfrag);

		uint64_t const len_a = simpleprocess(reffrag,readfrag,14,100,0,1,100);
		uint64_t const len_b = simpleprocess(reffrag,readfragrc,14,100,0,1,100);

		// std::cerr << sid << " " << seqindex << " " << refid << " " << restname << " " << from << " " << to << " len_a=" << len_a << " len_b=" << len_b << std::endl;

		libmaus2::lcs::NP np;
		bool rc;

		if ( len_a > len_b )
		{
			np.np(reffrag.begin(),reffrag.end(),readfrag.begin(),readfrag.end());
			rc = false;
		}
		else
		{
			np.np(reffrag.begin(),reffrag.end(),readfragrc.begin(),readfragrc.end());
			rc = true;
		}

		libmaus2::autoarray::AutoArray< std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Aop;

		uint64_t const ncig = libmaus2::bambam::CigarStringParser::traceToCigar(
			np,
			Aopblocks,
			Aop,
			0,0,0,0);

		uint64_t delshift = 0;
		uint64_t z = 0;
		while ( z < ncig && Aop[z].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL )
			delshift += Aop[z++].second;

		std::ostringstream cigstr;
		for ( uint64_t i = 0; i < ncig; ++i )
		{
			switch ( Aop[i].first )
			{
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
					cigstr << Aop[i].second << 'D';
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
					cigstr << Aop[i].second << 'I';
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
					cigstr << Aop[i].second << 'X';
					break;
				case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
					cigstr << Aop[i].second << '=';
					break;
				default:
					break;
			}
		}
		std::string const cig = cigstr.str();

		std::vector < char > Vqual(readfrag.size(),'H');

		libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment
		(
			UB,
			seqenc,
			sid,
			refid,
			from - 1 + delshift,
			255,
			rc ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0,
			cig,
			-1,
			-1,
			0,
			rc ? readfragrc : readfrag,
			std::string(Vqual.begin(),Vqual.end()),
			33,
			true
		);
		::libmaus2::bambam::MdStringComputationContext mdnmcontext;
		libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(UB.buffer,UB.length,mdnmcontext,reffrag.begin() /* itref */);
		libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(UB,"MD",mdnmcontext.md.get());
		libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(UB,"NM",'i',mdnmcontext.nm);

		writer->writeBamBlock(UB.buffer,UB.length);

		// std::cerr << np.getAlignmentStatistics() << std::endl;

		std::cerr << id << "\t" << refid << "," << rc << ":" << from << "-" << to << std::endl;
	}

	writer.reset();
}

/*
 parameters:

 -t : default number of logical cores, threads
 */

static std::string helpMessage(libmaus2::util::ArgParser const & /* arg */)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	#if 0
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	#endif

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
			std::cerr << "usage: " << arg.progname << " [options] ref.fasta reads.fasta\n";
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

			libmaus2::util::ArgInfo const arginfo(argc,argv);
			r = wgsimtobam(arg,arginfo);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			return r;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
