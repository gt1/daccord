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

// #define HANDLE_INDEL_ESTIMATE_DEBUG
// #define HANDLE_DEBUG
// #define HANDLE_DEEP_DEBUG
// #define CONSENSUS_SERIAL
// #define WINDOWALIGNAVOID
// #define CANDIDATE_AVOID_ALIGN
// #define CONS_HANDLE_SINGLE 0

#include <libmaus2/util/U.hpp>
#include <HandleContext.hpp>
#include <libmaus2/aio/DebugLineOutputStream.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapParser.hpp>
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
#include <libmaus2/util/MemUsage.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/wavelet/WaveletTree.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>

#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/lcs/NP.hpp>
#endif

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}


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

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static uint64_t getDefaultVerbose()
{
	return std::numeric_limits<uint64_t>::max();
}

static uint64_t getDefaultK()
{
	return 8;
}

static uint64_t getDefaultMaxInput()
{
	return 5000;
}

static uint64_t getDefaultVarD()
{
	return 0;
}

static uint64_t getDefaultMaxAlign()
{
	return std::numeric_limits<uint64_t>::max();
}

static unsigned int getDefaultWindowSize()
{
	return 40;
}

static unsigned int getDefaultAdvanceSize()
{
	return 10;
}

static unsigned int getDefaultProduceFull()
{
	return 0;
}

static unsigned int getDefaultMinWindowCoverage()
{
	return 3;
}

static uint64_t getDefaultMinWindowError()
{
	return std::numeric_limits<uint64_t>::max();
}

static uint64_t getDefaultMinLen()
{
	return 0;
}

static uint64_t getDefaultMinFilterFreq()
{
	return 0;
}

static uint64_t getDefaultMaxFilterFreq()
{
	return 2;
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
	optionMap . push_back ( std::pair < std::string, std::string >("w", formatRHS("window size",getDefaultWindowSize())));
	optionMap . push_back ( std::pair < std::string, std::string >("a", formatRHS("advance size",getDefaultAdvanceSize())));
	optionMap . push_back ( std::pair < std::string, std::string >("d", formatRHS("max depth",getDefaultMaxAlign())));
	optionMap . push_back ( std::pair < std::string, std::string >("f", formatRHS("produce full sequences",getDefaultProduceFull())));
	optionMap . push_back ( std::pair < std::string, std::string >("V", formatRHS("verbosity",getDefaultVerbose())));
	optionMap . push_back ( std::pair < std::string, std::string >("I", formatRHS("read interval",std::string("0,") + libmaus2::util::NumberSerialisation::formatNumber(std::numeric_limits<uint64_t>::max(),0))));
	optionMap . push_back ( std::pair < std::string, std::string >("J", formatRHS("reads part",std::string("0,1"))));
	optionMap . push_back ( std::pair < std::string, std::string >("E", formatRHS("error profile file name",std::string("input.las.eprof"))));
	optionMap . push_back ( std::pair < std::string, std::string >("m", formatRHS("minimum window coverage",getDefaultMinWindowCoverage())));
	optionMap . push_back ( std::pair < std::string, std::string >("e", formatRHS("maximum window error",getDefaultMinWindowError())));
	optionMap . push_back ( std::pair < std::string, std::string >("l", formatRHS("minimum length of output",getDefaultMinLen())));
	optionMap . push_back ( std::pair < std::string, std::string >("minfilterfreq", formatRHS("minimum k-mer filter frequency",getDefaultMinFilterFreq())));
	optionMap . push_back ( std::pair < std::string, std::string >("maxfilterfreq", formatRHS("maximum k-mer filter frequency",getDefaultMaxFilterFreq())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("D", formatRHS("maximum number of alignments considered per read",getDefaultMaxInput())));
	optionMap . push_back ( std::pair < std::string, std::string >("vard", formatRHS("maximum number of alignments considered per read",getDefaultVarD())));
	optionMap . push_back ( std::pair < std::string, std::string >("eprofonly", formatRHS("compute error profile only",std::string("disable"))));
	optionMap . push_back ( std::pair < std::string, std::string >("deepprofileonly", formatRHS("compute error distribution estimate",std::string("disable"))));
	optionMap . push_back ( std::pair < std::string, std::string >("k", formatRHS("kmer size",getDefaultK())));

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

#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
// load reference text from a FastA file
static std::vector<std::string> loadTextVector(std::string const & textfn)
{
	libmaus2::aio::InputStreamInstance ISI(textfn);
	libmaus2::fastx::StreamFastAReaderWrapper SFARW(ISI);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
	std::vector<std::string> Vtext;

	while ( SFARW.getNextPatternUnlocked(pattern) )
		Vtext.push_back(pattern.spattern);

	return Vtext;
}

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


template<unsigned int k>
double handleIndelEstimate(
	std::ostream &
		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		err
		#endif
		,
	uint64_t const maxalign,
	libmaus2::dazzler::align::Overlap const * ita,
	libmaus2::dazzler::align::Overlap const * ite,
	// window size for consensus computation
	uint64_t const windowsize,
	uint64_t const advancesize,
	// reads
	DecodedReadContainer & RC,
	DecodedReadContainer & RC2,
	// trace free list
	libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> & traceFreeList,
	// aligner
	libmaus2::lcs::Aligner & NP,
	// trace point space for overlaps
	int64_t const tspace,
	// kmer length
	// unsigned int const k,
	libmaus2::lcs::AlignmentStatistics & RGAS,
	uint64_t & usable,
	uint64_t & unusable
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
	,
	// BAM alignment of OVL-A read to reference
	libmaus2::bambam::BamAlignment const & algn,
	// reference text
	std::string const & text,
	//
	libmaus2::lcs::AlignmentStatistics & TAS,
	libmaus2::parallel::PosixSpinLock & TASlock
	#endif
)
{
	// number of overlaps with A read
	uint64_t const nintv = ite-ita;

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

	int64_t const aread = nintv ? ita[0].aread : -1;
	for ( uint64_t i = 1; i < nintv; ++i )
		assert ( ita[i].aread == aread );

	typedef std::pair<uint64_t,uint64_t> upair;

	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
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

	{
		std::string const R = algn.getRead();
		uint8_t const * ru = reinterpret_cast<uint8_t const *>(R.c_str()) + algn.getFrontSoftClipping();
		libmaus2::lcs::NP NP;
		NP.np(
			reinterpret_cast<uint8_t const *>(psreftext.c_str()),
			reinterpret_cast<uint8_t const *>(psreftext.c_str())+psreftext.size(),
			ru,ru + (R.size() - (algn.getFrontSoftClipping() + algn.getBackSoftClipping()))
		);
		libmaus2::lcs::AlignmentStatistics const AS = NP.getAlignmentStatistics();

		{
			libmaus2::parallel::ScopePosixSpinLock slock(TASlock);
			TAS += AS;
		}

		// std::cerr << "minimal " << AS << std::endl;

		#if 0
		libmaus2::lcs::AlignmentStatistics const ASe = Pbamtrace->getAlignmentStatistics();
		std::cerr << "actual " << ASe << std::endl;
		#endif
	}
	#endif

	// traces
	std::map<uint64_t,trace_type::shared_ptr_type> Mtraces;

	// repeat detector for length k-1
	libmaus2::fastx::KmerRepeatDetector KRD(k-1);

	// compute traces
	uint64_t maxaepos = 0;

	if ( nintv )
	{
		uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[0].aread));

		for ( uint64_t z = 0; z < nintv; ++z )
		{
			// get trace from free list
			trace_type::shared_ptr_type Ptrace = traceFreeList.get();

			// update maximum aepos
			if ( ita[z].path.aepos > static_cast<int64_t>(maxaepos) )
				maxaepos = ita[z].path.aepos;

			// get pointers to data
			uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread));

			// compute the trace operations
			ita[z].computeTrace(ua,ub,tspace,*Ptrace,NP);

			// store trace object
			Mtraces[z] = Ptrace;
		}
	}

	// end of interval heap
	libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > E(1024);
	// set of active traces/reads
	std::map < uint64_t, ActiveElement > activeset;
	libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > MA;

	// de bruijn graph for k
	KmerLimit KL(0.85,0);
	DebruijnGraph<k> DG(0,KL);

	// insert pointer
	uint64_t z = 0;

	uint64_t const ylimit = (maxaepos + advancesize >= windowsize) ? ((maxaepos + advancesize - windowsize) / advancesize) : 0;

	assert ( ylimit * advancesize + windowsize > maxaepos );

	double esum = 0;
	uint64_t ecnt = 0;

	// window id
	for ( uint64_t y = 0; y < ylimit; ++y )
	{
		// start of window on A
		uint64_t const astart = y * advancesize;
		uint64_t const aend = astart + windowsize;
		assert ( aend <= maxaepos );

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
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
				uint64_t const aoff = astart - ita[z].path.abpos;
				// err << "aoff=" << aoff << " ab=" << ita[z].path.abpos << std::endl;

				assert ( (ita[z].path.abpos + aoff) % advancesize == 0 );

				// get trace
				trace_type::shared_ptr_type Ptrace = Mtraces.find(z)->second;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = Ptrace->ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * te = Ptrace->te;

				// see how many operations there are up to start of region we are interested in
				std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(ta,te,aoff);
				assert ( adv.first == aoff );

				uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)) + ita[z].path.abpos + aoff;
				uint64_t const uboff = ita[z].path.bbpos + libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ta,ta+adv.second).second;
				uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread)) + uboff;

				// advance in trace
				ta += adv.second;

				uint64_t const escore = static_cast<uint64_t>(((ita[z].getErrorRate() - minerate) / ediv) * std::numeric_limits<uint32_t>::max());
				assert ( escore <= std::numeric_limits<uint32_t>::max() );
				uint64_t const eindex = (escore<<32) | z;

				activeset[eindex] = ActiveElement(ua,ub,ta,te,uboff,ita[z].getErrorRate());
				// add to erase list
				E.pushBump(upair(ita[z].path.aepos,eindex));
			}

			z += 1;
		}
		// cleanup
		while ( (! (E.empty())) && E.top().first <= aend )
		{
			upair UP = E.pop();
			activeset.erase(UP.second);
		}

		uint64_t MAo = 0;

		// iterate over active read list
		for ( std::map<uint64_t,ActiveElement>::iterator s_ita = activeset.begin(); s_ita != activeset.end(); ++s_ita )
		{
			// active element
			ActiveElement & AE = s_ita->second;

			#if ! defined(NDEBUG)
			// read id
			uint64_t const key = s_ita->first & std::numeric_limits<uint32_t>::max();
			// get overlap
			libmaus2::dazzler::align::Overlap const & OVL = ita[key];
			// sanity checks
			//std::cerr << "check " << OVL << std::endl;
			assert ( OVL.path.abpos <= static_cast<int64_t>(astart ) );
			assert ( OVL.path.aepos >= static_cast<int64_t>(aend   ) );
			#endif

			// get end of region in trace
			std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,windowsize);
			assert ( adv.first == windowsize );
			// string length
			std::pair<uint64_t,uint64_t> const sl = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+adv.second);

			// push A read if this is the first instance of this loop
			if ( (! MAo) && (&RC != &RC2) )
				MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ua,windowsize));
			// push B read
			if ( MAo < maxalign )
				MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ub,sl.second));

			std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,advancesize);
			std::pair<uint64_t,uint64_t> const sladv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+advadv.second);

			// update active element by advancing to next window
			AE.ua += advancesize;
			AE.ta += advadv.second;
			AE.ub += sladv.second;
			AE.uboff += sladv.second;
		}

		if ( MAo >= 3 )
		{
			bool ghasrep = false;

			for ( uint64_t i = 0; i < MAo; ++i )
			{
				// does read have a repeating k-1 mer?
				bool const hasrep = KRD.detect(MA[i].first,MA[i].second);
				ghasrep = ghasrep || hasrep;
			}

			if ( ghasrep )
			{
				unusable += 1;
			}
			else
			{
				usable += 1;

				#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
				// check whether ground truth has repeating k-1 mer
				bool const refhasrep = KRD.detect(refu,slref.first);
				// err << "ghasrep=" << ghasrep << " refhasrep=" << refhasrep << std::endl;
				#endif

				// set up debruijn graph
				DG.setup(MA.begin(), MAo);
				// filter by freuency (need at least 2)
				DG.filterFreq(2,MAo);
				// DG.computeFeasibleKmerPositions(OffsetLikely const & offsetLikely, double const thres = 1e-3)
				//DB.check();

				// traverse the graph
				bool const consok = DG.traverseTrivial();

				#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
				if ( ghasrep != refhasrep )
				{
					err << "GHASREP != REFHASREP for aread = " << ita[0].aread << " " << consok << std::endl;
					KRD.printRepeats();
				}
				#endif

				// if consensus looks ok
				if ( consok )
				{
					// get the consensus
					std::string const consensus = DG.getConsensus();

					libmaus2::lcs::AlignmentStatistics GAS;

					// align all reads to the consensus
					for ( uint64_t i = 0; i < MAo; ++i )
					{
						NP.align(
							reinterpret_cast<uint8_t const *>(consensus.c_str()),
							consensus.size(),
							MA[i].first,
							MA[i].second
						);

						// update error profile
						GAS += NP.getTraceContainer().getAlignmentStatistics();
					}

					RGAS += GAS;

					esum += GAS.getErrorRate();
					ecnt += 1;
				}
			}
		}

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		uint64_t const refadv = frontsoftclip >= advancesize ? 0 : advancesize - frontsoftclip;
		std::pair<uint64_t,uint64_t> const advrefadv = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refadv);
		std::pair<uint64_t,uint64_t> const slrefadv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advrefadv.second);

		frontsoftclip -= std::min(frontsoftclip,advancesize);

		rta += advrefadv.second;
		refu += slrefadv.first;
		#endif
	}

	for ( std::map<uint64_t,trace_type::shared_ptr_type>::const_iterator m_ita =
		Mtraces.begin(); m_ita != Mtraces.end(); ++m_ita )
		traceFreeList.put(m_ita->second);
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
	traceFreeList.put(Pbamtrace);
	#endif

	return ecnt ? (esum / ecnt) : 0.0;
}

template<unsigned int k>
uint64_t handleIndelEstimateDeep(
	libmaus2::sorting::SerialisingSortingBufferedOutputFileArray< libmaus2::util::U<uint32_t> >::sorter_type & usorter,
	std::ostream &
		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		err
		#endif
		,
	uint64_t const maxalign,
	libmaus2::dazzler::align::Overlap const * ita,
	libmaus2::dazzler::align::Overlap const * ite,
	// window size for consensus computation
	uint64_t const windowsize,
	uint64_t const advancesize,
	// reads
	DecodedReadContainer & RC,
	DecodedReadContainer & RC2,
	// trace free list
	libmaus2::parallel::LockedGrowingFreeList<trace_type,TraceAllocator,TraceTypeInfo> & traceFreeList,
	// aligner
	libmaus2::lcs::Aligner & NP,
	// trace point space for overlaps
	int64_t const tspace,
	// kmer length
	// unsigned int const k,
	libmaus2::lcs::AlignmentStatistics & RGAS,
	uint64_t & usable,
	uint64_t & unusable
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
	,
	// BAM alignment of OVL-A read to reference
	libmaus2::bambam::BamAlignment const & algn,
	// reference text
	std::string const & text,
	//
	libmaus2::lcs::AlignmentStatistics & TAS,
	libmaus2::parallel::PosixSpinLock & TASlock
	#endif
)
{
	uint64_t ucnt = 0;
	// number of overlaps with A read
	uint64_t const nintv = ite-ita;

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

	int64_t const aread = nintv ? ita[0].aread : -1;
	for ( uint64_t i = 1; i < nintv; ++i )
		assert ( ita[i].aread == aread );

	typedef std::pair<uint64_t,uint64_t> upair;

	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
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

	{
		std::string const R = algn.getRead();
		uint8_t const * ru = reinterpret_cast<uint8_t const *>(R.c_str()) + algn.getFrontSoftClipping();
		libmaus2::lcs::NP NP;
		NP.np(
			reinterpret_cast<uint8_t const *>(psreftext.c_str()),
			reinterpret_cast<uint8_t const *>(psreftext.c_str())+psreftext.size(),
			ru,ru + (R.size() - (algn.getFrontSoftClipping() + algn.getBackSoftClipping()))
		);
		libmaus2::lcs::AlignmentStatistics const AS = NP.getAlignmentStatistics();

		{
			libmaus2::parallel::ScopePosixSpinLock slock(TASlock);
			TAS += AS;
		}

		// std::cerr << "minimal " << AS << std::endl;

		#if 0
		libmaus2::lcs::AlignmentStatistics const ASe = Pbamtrace->getAlignmentStatistics();
		std::cerr << "actual " << ASe << std::endl;
		#endif
	}
	#endif

	// traces
	std::map<uint64_t,trace_type::shared_ptr_type> Mtraces;

	// repeat detector for length k-1
	libmaus2::fastx::KmerRepeatDetector KRD(k-1);

	// compute traces
	uint64_t maxaepos = 0;

	if ( nintv )
	{
		uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[0].aread));

		for ( uint64_t z = 0; z < nintv; ++z )
		{
			// get trace from free list
			trace_type::shared_ptr_type Ptrace = traceFreeList.get();

			// update maximum aepos
			if ( ita[z].path.aepos > static_cast<int64_t>(maxaepos) )
				maxaepos = ita[z].path.aepos;

			// get pointers to data
			uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread));

			// compute the trace operations
			ita[z].computeTrace(ua,ub,tspace,*Ptrace,NP);

			// store trace object
			Mtraces[z] = Ptrace;
		}
	}

	// end of interval heap
	libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > E(1024);
	// set of active traces/reads
	std::map < uint64_t, ActiveElement > activeset;
	libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > MA;

	// de bruijn graph for k
	KmerLimit KL(0.85,0);
	DebruijnGraph<k> DG(0,KL);

	// insert pointer
	uint64_t z = 0;

	uint64_t const ylimit = (maxaepos + advancesize >= windowsize) ? ((maxaepos + advancesize - windowsize) / advancesize) : 0;

	assert ( ylimit * advancesize + windowsize > maxaepos );

	// window id
	for ( uint64_t y = 0; y < ylimit; ++y )
	{
		// start of window on A
		uint64_t const astart = y * advancesize;
		uint64_t const aend = astart + windowsize;
		assert ( aend <= maxaepos );

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
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
				uint64_t const aoff = astart - ita[z].path.abpos;
				// err << "aoff=" << aoff << " ab=" << ita[z].path.abpos << std::endl;

				assert ( (ita[z].path.abpos + aoff) % advancesize == 0 );

				// get trace
				trace_type::shared_ptr_type Ptrace = Mtraces.find(z)->second;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * ta = Ptrace->ta;
				libmaus2::lcs::AlignmentTraceContainer::step_type const * te = Ptrace->te;

				// see how many operations there are up to start of region we are interested in
				std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(ta,te,aoff);
				assert ( adv.first == aoff );

				uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)) + ita[z].path.abpos + aoff;
				uint64_t const uboff = ita[z].path.bbpos + libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(ta,ta+adv.second).second;
				uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread)) + uboff;

				// advance in trace
				ta += adv.second;

				uint64_t const escore = static_cast<uint64_t>(((ita[z].getErrorRate() - minerate) / ediv) * std::numeric_limits<uint32_t>::max());
				assert ( escore <= std::numeric_limits<uint32_t>::max() );
				uint64_t const eindex = (escore<<32) | z;

				activeset[eindex] = ActiveElement(ua,ub,ta,te,uboff,ita[z].getErrorRate());
				// add to erase list
				E.pushBump(upair(ita[z].path.aepos,eindex));
			}

			z += 1;
		}
		// cleanup
		while ( (! (E.empty())) && E.top().first <= aend )
		{
			upair UP = E.pop();
			activeset.erase(UP.second);
		}

		uint64_t MAo = 0;

		// iterate over active read list
		for ( std::map<uint64_t,ActiveElement>::iterator s_ita = activeset.begin(); s_ita != activeset.end(); ++s_ita )
		{
			// active element
			ActiveElement & AE = s_ita->second;

			#if ! defined(NDEBUG)
			// read id
			uint64_t const key = s_ita->first & std::numeric_limits<uint32_t>::max();
			// get overlap
			libmaus2::dazzler::align::Overlap const & OVL = ita[key];
			// sanity checks
			//std::cerr << "check " << OVL << std::endl;
			assert ( OVL.path.abpos <= static_cast<int64_t>(astart ) );
			assert ( OVL.path.aepos >= static_cast<int64_t>(aend   ) );
			#endif

			// get end of region in trace
			std::pair<uint64_t,uint64_t> const adv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,windowsize);
			assert ( adv.first == windowsize );
			// string length
			std::pair<uint64_t,uint64_t> const sl = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+adv.second);

			// push A read if this is the first instance of this loop
			if ( (! MAo) && (&RC != &RC2) )
				MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ua,windowsize));
			// push B read
			if ( MAo < maxalign )
				MA.push(MAo,std::pair< uint8_t const *, uint64_t>(AE.ub,sl.second));

			std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,advancesize);
			std::pair<uint64_t,uint64_t> const sladv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(AE.ta,AE.ta+advadv.second);

			// update active element by advancing to next window
			AE.ua += advancesize;
			AE.ta += advadv.second;
			AE.ub += sladv.second;
			AE.uboff += sladv.second;
		}

		if ( MAo >= 3 )
		{
			bool ghasrep = false;

			for ( uint64_t i = 0; i < MAo; ++i )
			{
				// does read have a repeating k-1 mer?
				bool const hasrep = KRD.detect(MA[i].first,MA[i].second);
				ghasrep = ghasrep || hasrep;
			}

			if ( ghasrep )
			{
				unusable += 1;
			}
			else
			{
				usable += 1;

				#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
				// check whether ground truth has repeating k-1 mer
				bool const refhasrep = KRD.detect(refu,slref.first);
				// err << "ghasrep=" << ghasrep << " refhasrep=" << refhasrep << std::endl;
				#endif

				// set up debruijn graph
				DG.setup(MA.begin(), MAo);
				// filter by freuency (need at least 2)
				DG.filterFreq(2,MAo);
				// DG.computeFeasibleKmerPositions(OffsetLikely const & offsetLikely, double const thres = 1e-3)
				//DB.check();

				// traverse the graph
				bool const consok = DG.traverseTrivial();

				#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
				if ( ghasrep != refhasrep )
				{
					err << "GHASREP != REFHASREP for aread = " << ita[0].aread << " " << consok << std::endl;
					KRD.printRepeats();
				}
				#endif

				// if consensus looks ok
				if ( consok )
				{
					// get the consensus
					std::string const consensus = DG.getConsensus();

					libmaus2::lcs::AlignmentStatistics GAS;

					// align all reads to the consensus
					for ( uint64_t i = 0; i < MAo; ++i )
					{
						NP.align(
							reinterpret_cast<uint8_t const *>(consensus.c_str()),
							consensus.size(),
							MA[i].first,
							MA[i].second
						);

						// update error profile
						GAS += NP.getTraceContainer().getAlignmentStatistics();
					}

					RGAS += GAS;

					double const e = GAS.getErrorRate();
					uint64_t const v = static_cast<uint64_t>(static_cast<double>(std::numeric_limits<uint32_t>::max()) * e + 0.5);
					uint32_t const v32 = std::min(v,static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()));
					usorter.put(v32);
					ucnt += 1;
				}
			}
		}

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		uint64_t const refadv = frontsoftclip >= advancesize ? 0 : advancesize - frontsoftclip;
		std::pair<uint64_t,uint64_t> const advrefadv = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refadv);
		std::pair<uint64_t,uint64_t> const slrefadv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advrefadv.second);

		frontsoftclip -= std::min(frontsoftclip,advancesize);

		rta += advrefadv.second;
		refu += slrefadv.first;
		#endif
	}

	for ( std::map<uint64_t,trace_type::shared_ptr_type>::const_iterator m_ita =
		Mtraces.begin(); m_ita != Mtraces.end(); ++m_ita )
		traceFreeList.put(m_ita->second);
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
	traceFreeList.put(Pbamtrace);
	#endif

	return ucnt;
}

struct OverlapPosComparator
{
	bool operator()(libmaus2::dazzler::align::Overlap const & lhs, libmaus2::dazzler::align::Overlap const & rhs) const
	{
		return lhs.path.abpos < rhs.path.abpos;
	}
};

struct OverlapErrorDescComparator
{
	bool operator()(libmaus2::dazzler::align::Overlap const & lhs, libmaus2::dazzler::align::Overlap const & rhs) const
	{
		return lhs.getErrorRate() < rhs.getErrorRate();
	}
};


struct ErrorInfo
{
	typedef ErrorInfo this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::lcs::AlignmentStatistics AS;

	ErrorInfo()
	{

	}
	ErrorInfo(libmaus2::lcs::AlignmentStatistics const & rAS)
	: AS(rAS)
	{

	}
	ErrorInfo(std::istream & in)
	{
		deserialise(in);
	}
	ErrorInfo(std::string const & fn)
	{
		deserialise(fn);
	}

	std::ostream & serialise(std::ostream & out) const
	{
		return AS.serialise(out);
	}

	void serialise(std::string const & fn) const
	{
		return AS.serialise(fn);
	}

	std::istream & deserialise(std::istream & in)
	{
		return AS.deserialise(in);
	}

	void deserialise(std::string const & fn)
	{
		AS.deserialise(fn);
	}
};

int daccord(
	libmaus2::util::ArgParser const & arg,
	std::string const & lasfn,
	std::string const & dbname,
	std::string const & db2name
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
		,
	std::string const & textfn,
	std::string const & bamfn
	#endif
)
{
	arg.printArgs(std::cerr,std::string("[V] "));

	std::string const lasindexname = libmaus2::dazzler::align::OverlapIndexer::getIndexName(lasfn);
	std::string const dalignerlasindexname = libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(lasfn);

	if (
		! libmaus2::util::GetFileSize::fileExists(lasindexname)
		||
		libmaus2::util::GetFileSize::isOlder(lasindexname,lasfn)
	)
	{
		libmaus2::dazzler::align::OverlapIndexer::constructIndex(lasfn,&std::cerr);
	}

	if (
		! libmaus2::util::GetFileSize::fileExists(dalignerlasindexname)
		||
		libmaus2::util::GetFileSize::isOlder(dalignerlasindexname,lasfn)
	)
	{
		libmaus2::dazzler::align::OverlapIndexer::constructIndex(lasfn,&std::cerr);
	}

	std::string const memlasindexname = "mem:.input.las.bidx";
	std::cerr << "[V] copying " << lasindexname << " to " << memlasindexname << "...";
	libmaus2::util::GetFileSize::copy(lasindexname,memlasindexname);
	std::cerr << "done." << std::endl;

	std::string const memdalignerlasindexname = "mem:.input.las.idx";
	std::cerr << "[V] copying " << dalignerlasindexname << " to " << memdalignerlasindexname << "...";
	libmaus2::util::GetFileSize::copy(dalignerlasindexname,memdalignerlasindexname);
	std::cerr << "done." << std::endl;

	int64_t minaread;
	int64_t maxaread;

	bool const producefull = arg.uniqueArgPresent("f");

	#if defined(CONS_HANDLE_SINGLE)
	minaread = CONS_HANDLE_SINGLE;
	maxaread = CONS_HANDLE_SINGLE;
	#else

	minaread = libmaus2::dazzler::align::OverlapIndexer::getMinimumARead(lasfn);
	maxaread = libmaus2::dazzler::align::OverlapIndexer::getMaximumARead(lasfn);

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

	uint64_t isize = (toparead < 0) ? 0 : (toparead - minaread);
	libmaus2::autoarray::AutoArray<uint64_t> RL(isize,true);

	std::cerr << "[V] minaread=" << minaread << " toparead=" << toparead << std::endl;

	uint64_t kmersizelow = getDefaultK();
	uint64_t kmersizehigh = getDefaultK();

	if ( arg.uniqueArgPresent("k") )
	{
		std::string const karg = arg["k"];

		if ( karg.find(",") != std::string::npos )
		{
			std::istringstream istr(karg);
			istr >> kmersizelow;
			if ( (!istr) || istr.peek() != ',' )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse k argument " << karg << std::endl;
				lme.finish();
				throw lme;
			}
			istr.get(); // ','
			istr >> kmersizehigh;
			if ( (!istr) || istr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse k argument " << karg << std::endl;
				lme.finish();
				throw lme;
			}
		}
		else
		{
			std::istringstream istr(karg);
			istr >> kmersizelow;

			if ( (!istr) || istr.peek() != std::istream::traits_type::eof() )
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "[E] unable to parse k argument " << karg << std::endl;
				lme.finish();
				throw lme;
			}

			kmersizehigh = kmersizelow;
		}
	}

	uint64_t const indel_estimate_windowsize = 40;
	uint64_t const indel_estimate_advancesize = 5;

	uint64_t const windowsize = arg.uniqueArgPresent("w") ? arg.getUnsignedNumericArg<uint64_t>("w") : getDefaultWindowSize(); // 100, 40;
	uint64_t const advancesize = arg.uniqueArgPresent("a") ? arg.getUnsignedNumericArg<uint64_t>("a") : getDefaultAdvanceSize(); // 25, 5;
	uint64_t const maxinput = arg.uniqueArgPresent("D") ? arg.getUnsignedNumericArg<uint64_t>("D") : getDefaultMaxInput();

	// number of threads
	uint64_t const numthreads =
		#if defined(CONSENSUS_SERIAL)
		1
		#else
		arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads()
		#endif
		;
	uint64_t const verbose = arg.uniqueArgPresent("V") ? arg.getUnsignedNumericArg<uint64_t>("V") : getDefaultVerbose();
	uint64_t const maxalign = arg.uniqueArgPresent("d") ? arg.getUnsignedNumericArg<uint64_t>("d") : getDefaultMaxAlign();
	uint64_t const keepeprof = arg.uniqueArgPresent("keepeprof") ? arg.getUnsignedNumericArg<uint64_t>("keepeprof") : 1;
	uint64_t const minwindowcov = arg.uniqueArgPresent("m") ? arg.getUnsignedNumericArg<uint64_t>("m") : getDefaultMinWindowCoverage();
	uint64_t const eminrate = arg.uniqueArgPresent("e") ? arg.getUnsignedNumericArg<uint64_t>("e") : getDefaultMinWindowError();
	uint64_t const minlen = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultMinLen();
	uint64_t const vard = arg.uniqueArgPresent("vard") ? arg.getUnsignedNumericArg<uint64_t>("vard") : getDefaultVarD();
	bool const eprofonly = arg.uniqueArgPresent("eprofonly");
	bool const deepprofileonly = arg.uniqueArgPresent("deepprofileonly");

	uint64_t const minfilterfreq = arg.uniqueArgPresent("minfilterfreq") ? arg.getUnsignedNumericArg<uint64_t>("minfilterfreq") : getDefaultMinFilterFreq();
	uint64_t const maxfilterfreq = arg.uniqueArgPresent("maxfilterfreq") ? arg.getUnsignedNumericArg<uint64_t>("maxfilterfreq") : getDefaultMaxFilterFreq();

	std::cerr << "[V] minfilterfreq=" << minfilterfreq << " maxfilterfreq=" << maxfilterfreq << std::endl;

	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type> Adalindex(numthreads);

	for ( uint64_t i = 0; i < Adalindex.size(); ++i )
	{
		libmaus2::dazzler::align::DalignerIndexDecoder::unique_ptr_type Pdalindex(
			new libmaus2::dazzler::align::DalignerIndexDecoder(lasfn,memdalignerlasindexname)
		);
		Adalindex[i] = UNIQUE_PTR_MOVE(Pdalindex);
	}

	libmaus2::autoarray::AutoArray < libmaus2::lcs::Aligner::unique_ptr_type > ANP(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::lcs::Aligner::unique_ptr_type tptr(DebruijnGraphBase::getAligner());
		ANP[i] = UNIQUE_PTR_MOVE(tptr);
	}

	std::cerr << "[V] copying " << dbname << " to memory...";
	// libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbname,"mem:db1prefix"));
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbname)
	);
	std::cerr << "done." << std::endl;
	std::cerr << "[V] loading and trimming " << dbname << "...";
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
	PDB->computeTrimVector();
	std::cerr << "done." << std::endl;

	libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type db2ptr;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB2;

	if ( db2name != dbname )
	{
		if ( vard )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[V] vard option is not supported for assyemtric (DB1 != DB2) alignments" << std::endl;
			lme.finish();
			throw lme;
		}

		std::cerr << "[V] copying " << db2name << " to memory...";
		libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type tdbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(db2name,"mem:db2prefix"));
		std::cerr << "done." << std::endl;
		db2ptr = UNIQUE_PTR_MOVE(tdbptr);

		std::cerr << "[V] loading and trimming " << db2name << "...";
		libmaus2::dazzler::db::DatabaseFile::unique_ptr_type TDB(new libmaus2::dazzler::db::DatabaseFile(db2ptr->fn));
		TDB->computeTrimVector();
		std::cerr << "done." << std::endl;

		PDB2 = UNIQUE_PTR_MOVE(TDB);
	}

	libmaus2::dazzler::db::DatabaseFile const & DB = *PDB;
	libmaus2::dazzler::db::DatabaseFile const & DB2 = PDB2 ? (*PDB2) : (*PDB);

	int64_t const avgreadlength = DB2.getAverageReadLength();

	if ( maxaread >= 0 )
		DB2.getReadLengthArrayParallel(minaread,toparead,RL.begin(),numthreads);

	//libmaus2::fastx::ReadContainer const RC(readsfn);
	// std::cerr << "opening " << dbptr->fn << std::endl;

	// load read/read alignments
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(lasfn);

	std::cerr << "[V] setting up read decoder allocators...";
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
	std::cerr << "done." << std::endl;

	// libmaus2::parallel::LockedGrowingFreeList<libmaus2::dazzler::align::Overlap,OverlapAllocator,OverlapTypeInfo> overlapFreeList;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > > AAOVL(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
		AAOVL[i] = libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap >(maxinput);

	struct PairFirstGreaterComp
	{
		bool operator()(std::pair<uint64_t,uint64_t> const & A, std::pair<uint64_t,uint64_t> const & B) const
		{
			return A.first > B.first;
		}
	};

	libmaus2::autoarray::AutoArray< libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp >::unique_ptr_type > AAH(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp >::unique_ptr_type tptr(
			new libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp >(maxinput)
		);
		AAH[i] = UNIQUE_PTR_MOVE(tptr);
	}

	std::string const eproffn = arg.uniqueArgPresent("E") ? arg["E"] : (lasfn + ".eprof");
	std::string const tmpfilebase = getTmpFileBase(arg);
	std::string const eproffntmp = tmpfilebase + "_eprof.tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(eproffntmp);
	libmaus2::lcs::AlignmentStatistics GAS;

	// reference sequence text
	#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
	// std::string const text = loadText(textfn);
	std::vector<std::string> const Vtext = loadTextVector(textfn);

	// read alignments
	std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> Vbam = readBamAlignments(bamfn);
	#endif

	libmaus2::autoarray::AutoArray<double> Aloc;
	uint64_t volatile oAloc = 0;

	// compute error profile if it does not already exist
	if ( deepprofileonly )
	{
		libmaus2::parallel::PosixSpinLock GASlock;
		libmaus2::parallel::PosixSpinLock ucntlock;
		uint64_t volatile ucnt = 0;

		libmaus2::sorting::SerialisingSortingBufferedOutputFileArray< libmaus2::util::U<uint32_t> >::unique_ptr_type pSSBOF(
			new libmaus2::sorting::SerialisingSortingBufferedOutputFileArray< libmaus2::util::U<uint32_t> >(tmpfilebase + "_deep",numthreads)
		);

		typedef std::map<uint64_t,uint64_t> copy_freq_type;
		typedef libmaus2::util::shared_ptr<copy_freq_type>::type copy_freq_ptr_type;
		std::vector <copy_freq_ptr_type> Vcopyfreq(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
			Vcopyfreq[i] = copy_freq_ptr_type(new copy_freq_type);

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		libmaus2::lcs::AlignmentStatistics TAS;
		libmaus2::parallel::PosixSpinLock TASlock;
		#endif

		// uint64_t volatile usable = 0;
		// uint64_t volatile unusable = 0;

		int64_t const top = std::min(toparead,minaread+4096);

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
		#endif
		for ( int64_t z = minaread; z < top; ++z )
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			libmaus2::sorting::SerialisingSortingBufferedOutputFileArray< libmaus2::util::U<uint32_t> >::sorter_type & usorter =
				(*pSSBOF)[tid];

			uint64_t const rl = RL[z-minaread];
			uint64_t const lmaxinput = vard ?
				std::max(
					std::max(2*vard,static_cast<uint64_t>(1)),
					static_cast<uint64_t>((2.0*static_cast<double>(vard)*static_cast<double>(rl)/static_cast<double>(avgreadlength))+0.5)
				)
				: maxinput;

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp > & RH = *AAH[tid];
			RH.clear();
			libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

			if ( lmaxinput > RH.H.size() )
				RH.ensureSize(lmaxinput);
			if ( lmaxinput > RO.size() )
				RO.ensureSize(lmaxinput);

			libmaus2::lcs::Aligner & al = *(ANP[tid]);

			// std::vector<libmaus2::dazzler::align::Overlap> VOVL;
			libmaus2::dazzler::align::Overlap OVL;

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
			uint64_t f = 0;
			while ( pdec->getNextOverlap(OVL) )
			{
				uint64_t const score =
					static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

				if ( RH.f == lmaxinput )
				{
					assert ( ! RH.empty() );

					// dump alignment if error rate is higher than highest in heap
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
					RH.push(std::pair<uint64_t,uint64_t>(score,p));
					RO[p] = OVL;
				}

				// VOVL.push_back(OVL);
			}

			// std::cerr << "f=" << f << std::endl;

			// std::cerr << "VOVL.size()=" << VOVL.size() << std::endl;

			if ( f )
			{
				#if 0
				if ( VOVL.size() > maxalign )
				{
					std::sort(VOVL.begin(),VOVL.end(),OverlapErrorDescComparator());
					VOVL.resize(maxalign);
				}
				#endif

				//std::sort(VOVL.begin(),VOVL.end(),OverlapPosComparator());
				std::sort(RO.begin(),RO.begin()+f,OverlapPosComparator());

				libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);

				DLOS << "[S] handling aread " << RO[0].aread << " depth " << f << std::endl;

				// std::pair<unsigned char const *, unsigned char const *> const Q = inqual.getQualityForRead(Vovl[low].aread);

				#if 0
				for ( unsigned char const * p = Q.first; p != Q.second; ++p )
					std::cerr << "block " << p-Q.first << " qual " << static_cast<int>(*p) << std::endl;
				#endif

				libmaus2::lcs::AlignmentStatistics LGAS;

				uint64_t lusable = 0, lunusable = 0;
				DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
				DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);

				// double eloc = 0.0;
				// if ( VOVL[0].aread == 137 )
				{
					#if 0
					char const * f = RDC.getForwardRead(VOVL[0].aread);
					uint64_t const l = RDC.getReadLength(VOVL[0].aread);
					std::cerr << std::string(f,f+l) << std::endl;
					#endif

					uint64_t const lcnt = handleIndelEstimateDeep<8>(
						usorter,
						DLOS,
						maxalign,
						// actual data
						RO.begin(),RO.begin()+f,
						indel_estimate_windowsize,indel_estimate_advancesize,
						RDC,
						RDC2,
						#if 0
						Vreads,Rreads,
						#endif
						traceFreeList,
						al,
						tspace,
						LGAS,
						lusable,
						lunusable
						#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
						,
						// for debugging
						*(Vbam[RO[0].aread]),
						Vtext[Vbam[RO[0].aread]->getRefID()],
						TAS,
						TASlock
						#endif
					);

					ucntlock.lock();
					ucnt += lcnt;
					ucntlock.unlock();
				}
			}
		}

		{
			libmaus2::sorting::SerialisingSortingBufferedOutputFileArray< libmaus2::util::U<uint32_t> >::merger_ptr_type pmerger(
				pSSBOF->getMergerParallel(numthreads)
			);
			libmaus2::util::U<uint32_t> D;
			int64_t prevD = -1;
			uint64_t n = 0;
			uint64_t s = 0;
			double const dcnt = ucnt;
			while ( pmerger->getNext(D) )
			{
				if ( static_cast<int64_t>(D.u) != prevD )
				{
					if ( n )
					{
						s += n;
						std::cout << "[deep]\t" << (prevD / static_cast<double>(std::numeric_limits<uint32_t>::max())) << "\t" << s/dcnt << "\n";
					}

					prevD = static_cast<int64_t>(D.u);
					n = 0;
				}

				n += 1;
			}

			if ( n )
			{
				s += n;
				std::cout << "[deep]\t" << (prevD / static_cast<double>(std::numeric_limits<uint32_t>::max())) << "\t" << s/dcnt << "\n";
			}
		}

		return EXIT_SUCCESS;
	}

	// compute error profile if it does not already exist
	if (
		(!libmaus2::util::GetFileSize::fileExists(eproffn))
		||
		(libmaus2::util::GetFileSize::isOlder(eproffn,lasfn) && (!keepeprof))
	)
	{
		libmaus2::parallel::PosixSpinLock GASlock;

		typedef std::map<uint64_t,uint64_t> copy_freq_type;
		typedef libmaus2::util::shared_ptr<copy_freq_type>::type copy_freq_ptr_type;
		std::vector <copy_freq_ptr_type> Vcopyfreq(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
			Vcopyfreq[i] = copy_freq_ptr_type(new copy_freq_type);

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		libmaus2::lcs::AlignmentStatistics TAS;
		libmaus2::parallel::PosixSpinLock TASlock;
		#endif

		uint64_t volatile usable = 0;
		uint64_t volatile unusable = 0;

		int64_t const top = std::min(toparead,minaread+1024);

		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
		#endif
		for ( int64_t z = minaread; z < top; ++z )
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			uint64_t const rl = RL[z-minaread];
			uint64_t const lmaxinput = vard ?
				std::max(
					std::max(2*vard,static_cast<uint64_t>(1)),
					static_cast<uint64_t>((2.0*static_cast<double>(vard)*static_cast<double>(rl)/static_cast<double>(avgreadlength))+0.5)
				)
				: maxinput;

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp > & RH = *AAH[tid];
			RH.clear();
			libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

			if ( lmaxinput > RH.H.size() )
				RH.ensureSize(lmaxinput);
			if ( lmaxinput > RO.size() )
				RO.ensureSize(lmaxinput);

			libmaus2::lcs::Aligner & al = *(ANP[tid]);

			// std::vector<libmaus2::dazzler::align::Overlap> VOVL;
			libmaus2::dazzler::align::Overlap OVL;

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
			uint64_t f = 0;
			while ( pdec->getNextOverlap(OVL) )
			{
				uint64_t const score =
					static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

				if ( RH.f == lmaxinput )
				{
					assert ( ! RH.empty() );

					// dump alignment if error rate is higher than highest in heap
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
					RH.push(std::pair<uint64_t,uint64_t>(score,p));
					RO[p] = OVL;
				}

				// VOVL.push_back(OVL);
			}

			// std::cerr << "f=" << f << std::endl;

			// std::cerr << "VOVL.size()=" << VOVL.size() << std::endl;

			if ( f )
			{
				#if 0
				if ( VOVL.size() > maxalign )
				{
					std::sort(VOVL.begin(),VOVL.end(),OverlapErrorDescComparator());
					VOVL.resize(maxalign);
				}
				#endif

				//std::sort(VOVL.begin(),VOVL.end(),OverlapPosComparator());
				std::sort(RO.begin(),RO.begin()+f,OverlapPosComparator());

				libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);

				DLOS << "[S] handling aread " << RO[0].aread << " depth " << f << std::endl;

				// std::pair<unsigned char const *, unsigned char const *> const Q = inqual.getQualityForRead(Vovl[low].aread);

				#if 0
				for ( unsigned char const * p = Q.first; p != Q.second; ++p )
					std::cerr << "block " << p-Q.first << " qual " << static_cast<int>(*p) << std::endl;
				#endif

				libmaus2::lcs::AlignmentStatistics LGAS;

				uint64_t lusable = 0, lunusable = 0;
				DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
				DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);

				double eloc = 0.0;
				// if ( VOVL[0].aread == 137 )
				{
					#if 0
					char const * f = RDC.getForwardRead(VOVL[0].aread);
					uint64_t const l = RDC.getReadLength(VOVL[0].aread);
					std::cerr << std::string(f,f+l) << std::endl;
					#endif

					eloc = handleIndelEstimate<8>(
						DLOS,
						maxalign,
						// actual data
						RO.begin(),RO.begin()+f,
						indel_estimate_windowsize,indel_estimate_advancesize,
						RDC,
						RDC2,
						#if 0
						Vreads,Rreads,
						#endif
						traceFreeList,
						al,
						tspace,
						LGAS,
						lusable,
						lunusable
						#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
						,
						// for debugging
						*(Vbam[RO[0].aread]),
						Vtext[Vbam[RO[0].aread]->getRefID()],
						TAS,
						TASlock
						#endif
					);
				}

				GASlock.lock();
				GAS += LGAS;
				usable += lusable;
				unusable += lunusable;
				uint64_t toAloc = oAloc;
				if ( eloc != 0.0 )
					Aloc.push(toAloc,eloc);
				oAloc = toAloc;
				GASlock.unlock();
			}
		}

		double esum = 0;
		uint64_t ecnt = 0;
		for ( uint64_t i = 0; i < oAloc; ++i )
		{
			esum += Aloc[i];
			ecnt += 1;
		}

		double eavg = 0.0;
		double edif = 0.0;
		if ( ecnt )
		{
			eavg = esum / ecnt;

			for ( uint64_t i = 0; i < oAloc; ++i )
				edif += (eavg - Aloc[i]) * (eavg - Aloc[i]);

			edif /= ecnt;
			edif = std::sqrt(edif);
		}

		std::cerr << "usable=" << usable << " unusable=" << unusable << " eavg=" << eavg << " edif=" << edif << std::endl;

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		std::cerr << "TAS=" << TAS << std::endl;
		#endif

		{
		libmaus2::aio::OutputStreamInstance OSI(eproffntmp);
		GAS.serialise(OSI);
		libmaus2::util::NumberSerialisation::serialiseDouble(OSI,eavg);
		libmaus2::util::NumberSerialisation::serialiseDouble(OSI,edif);
		}
		libmaus2::aio::OutputStreamFactoryContainer::rename(eproffntmp,eproffn);
	}
	else
	{
		GAS.deserialise(eproffn);
	}

	uint64_t const len = GAS.matches + GAS.mismatches + GAS.deletions;
	uint64_t const numerr = GAS.mismatches + GAS.deletions + GAS.insertions;

	double const est_erate = static_cast<double>(numerr) / len;
	double const est_cor = 1.0 - est_erate;
	double const est_i_frac = GAS.insertions / static_cast<double>(numerr);
	double const est_d_frac = GAS.deletions / static_cast<double>(numerr);
	double const est_s_frac = GAS.mismatches / static_cast<double>(numerr);

	double const p_i = static_cast<double>(GAS.insertions) / len;
	double const p_d = static_cast<double>(GAS.deletions) / len;
	double const p_s = static_cast<double>(GAS.mismatches) / len;

	std::cerr << GAS << std::endl;
	std::cerr << "len=" << len << std::endl;
	std::cerr << "error estimates" << std::endl;
	std::cerr << "erate=" << est_erate << std::endl;
	std::cerr << "cor=" << est_cor << std::endl;
	std::cerr << "ins=" << est_i_frac << " " << p_i << std::endl;
	std::cerr << "del=" << est_d_frac << " " << p_d << std::endl;
	std::cerr << "subst=" << est_s_frac << " " << p_s << std::endl;

	for ( uint64_t k = kmersizelow; k <= kmersizehigh; ++k )
	{
		double const est_cor_k = libmaus2::math::gpow(est_cor,k);
		std::cerr << "cor^" << k << "=" << est_cor_k << std::endl;
	}

	if ( eprofonly )
		return EXIT_SUCCESS;

	#if 0
	for ( uint64_t i = 0; i < 500; ++i )
	{
		std::cerr << i << "\t" << libmaus2::math::Binom::binomRowUpperLimit(est_cor_k, i, 0.98) << std::endl;
	}
	#endif

	#if 0
	double const i_ERATE=0.15;
	double const i_INSRATE=0.8;
	double const i_DELRATE=0.1333;
	double const i_MISRATE=0.0666;
	#endif

	// compute offset likelihood distributions
	OffsetLikely const DP = computeOffsetLikely(windowsize,p_i,p_d);
	//std::vector < DotProduct > const DP = computeOffsetLikely(windowsize,i_ERATE*i_INSRATE,i_ERATE*i_DELRATE);

	#if 0
	DP.plotData("theplot");

	for ( uint64_t i = 0; i <= 20; ++i )
	{
		std::ostringstream namestr;
		namestr << "theplot_" << i;
		DP.plotDataSingle(namestr.str(),i);
	}
	#endif

	#if 0
	for ( uint64_t i = 0; i < DP.size(); ++i )
	{
		std::cerr << "i=" << i << std::endl;
		std::cerr << DP[i] << std::endl;
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

	#if defined(CONS_HANDLE_SINGLE)
	libmaus2::dazzler::align::AlignmentWriter AW("debug_cons.las",libmaus2::dazzler::align::AlignmentFile::getTSpace(lasfn),false);
	#endif

	#if defined(HANDLE_DEBUG)
	typedef std::map<double,uint64_t> ed_type;
	typedef libmaus2::util::unique_ptr<ed_type>::type ed_ptr_type;
	libmaus2::autoarray::AutoArray < ed_ptr_type > AED(numthreads);
	for ( uint64_t i = 0; i < AED.size(); ++i )
	{
		ed_ptr_type tptr(new ed_type);
		AED[i] = UNIQUE_PTR_MOVE(tptr);
	}
	#endif

	#if 0
	std::vector<uint64_t> DGKV;
	for ( uint64_t k = kmersizelow; k <= kmersizehigh; ++k )
		DGKV.push_back(k);
	std::sort(DGKV.begin(),DGKV.end(),std::greater<uint64_t>());
	#endif

	std::cerr << "[V] using kmer range [" << kmersizelow << "," << kmersizehigh << "]" << std::endl;

	std::map < uint64_t, KmerLimit::shared_ptr_type > MKL;
	for ( uint64_t k = kmersizelow; k <= kmersizehigh; ++k )
	{
		std::cerr << "[V] creating KmerLimit for k=" << k << "...";
		KmerLimit::shared_ptr_type tptr(new KmerLimit(::std::pow(est_cor,k),100));
		MKL [ k ] = tptr;
		std::cerr << "done." << std::endl;
	}
	libmaus2::autoarray::AutoArray<HandleContext::unique_ptr_type> AHC(numthreads);
	#if defined(_OPENMP)
	#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
	#endif
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		HandleContext::unique_ptr_type tptr(
			new HandleContext(
				maxalign,
				wellcounter,
				windowsize,
				advancesize,
				readDataFreeList,
				readDecoderFreeList,
				readDecoderFreeList2,
				traceFreeList,
				tspace,
				DP,
				producefull,
				est_cor,
				kmersizelow,
				kmersizehigh,
				verbose,
				minwindowcov,
				eminrate,
				minlen,
				minfilterfreq,
				maxfilterfreq,
				MKL,
				1/*numthreads*/
			)
		);

		AHC[i] = UNIQUE_PTR_MOVE(tptr);
	}
	std::cerr << "[V] graph structures set up " << libmaus2::util::MemUsage() << std::endl;

	struct OverlapEntry
	{
		uint64_t score;
		uint64_t blockid;
		uint64_t entryid;
		uint64_t s;

		OverlapEntry() : score(0), blockid(0), entryid(0) {}
		OverlapEntry(
			uint64_t const rscore,
			uint64_t const rblockid,
			uint64_t const rentryid,
			uint64_t const rs
		) : score(rscore), blockid(rblockid), entryid(rentryid), s(rs) {}

		bool operator<(OverlapEntry const & O) const
		{
			return score < O.score;
		}
	};
	struct OverlapEntryBlockComparator
	{
		bool operator()(OverlapEntry const & A, OverlapEntry const & B) const
		{
			if ( A.blockid != B.blockid )
				return A.blockid < B.blockid;
			else
				return A.entryid < B.entryid;
		}
	};
	libmaus2::autoarray::AutoArray< libmaus2::util::FiniteSizeHeap< OverlapEntry >::unique_ptr_type > AAHO(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::util::FiniteSizeHeap< OverlapEntry >::unique_ptr_type tptr(
			new libmaus2::util::FiniteSizeHeap< OverlapEntry >(maxinput)
		);
		AAHO[i] = UNIQUE_PTR_MOVE(tptr);
	}

	uint64_t const inputbuffersize = 64*1024;
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray< uint8_t >::unique_ptr_type > Ainputbuffers(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::autoarray::AutoArray< uint8_t >::unique_ptr_type tptr(
			new libmaus2::autoarray::AutoArray< uint8_t >(inputbuffersize)
		);
		Ainputbuffers[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray< libmaus2::aio::InputStreamInstance::unique_ptr_type > Ainput(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::aio::InputStreamInstance::unique_ptr_type tptr(
			new libmaus2::aio::InputStreamInstance(lasfn)
		);
		Ainput[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray< libmaus2::dazzler::align::OverlapParser::unique_ptr_type > Aparsers(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::dazzler::align::OverlapParser::unique_ptr_type tptr(
			new libmaus2::dazzler::align::OverlapParser(tspace)
		);
		Aparsers[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray< uint8_t >::unique_ptr_type > Acopybuffers(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::autoarray::AutoArray< uint8_t >::unique_ptr_type tptr(
			new libmaus2::autoarray::AutoArray< uint8_t >()
		);
		Acopybuffers[i] = UNIQUE_PTR_MOVE(tptr);
	}
	libmaus2::autoarray::AutoArray< libmaus2::autoarray::AutoArray< libmaus2::dazzler::align::OverlapDataInterface >::unique_ptr_type > Acopypointers(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::autoarray::AutoArray< libmaus2::dazzler::align::OverlapDataInterface >::unique_ptr_type tptr(
			new libmaus2::autoarray::AutoArray< libmaus2::dazzler::align::OverlapDataInterface >()
		);
		Acopypointers[i] = UNIQUE_PTR_MOVE(tptr);
	}

	#if defined(_OPENMP)
		#if ! defined(CONSENSUS_SERIAL)
		#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
		#endif
	#endif
	for ( int64_t z = minaread; z < toparead; ++z )
	{
		#if defined(_OPENMP)
		uint64_t const tid = omp_get_thread_num();
		#else
		uint64_t const tid = 0;
		#endif

		uint64_t const rl = RL[z-minaread];
		uint64_t const lmaxinput = vard ?
			std::max(
				std::max(2*vard,static_cast<uint64_t>(1)),
				static_cast<uint64_t>((2.0*static_cast<double>(vard)*static_cast<double>(rl)/static_cast<double>(avgreadlength))+0.5)
			)
			: maxinput;


		libmaus2::util::FiniteSizeHeap< OverlapEntry > & RHO = *AAHO[tid];
		RHO.clear();
		if ( lmaxinput > RHO.H.size() )
			RHO.ensureSize(lmaxinput);
		uint64_t const ofrom = (*(Adalindex[tid]))[z];
		uint64_t const oto = (*(Adalindex[tid]))[z+1];
		std::istream & istr = *(Ainput[tid]);
		istr.clear();
		istr.seekg(ofrom);
		libmaus2::dazzler::align::OverlapParser & parser = *(Aparsers[tid]);
		assert ( parser.isIdle() );
		libmaus2::autoarray::AutoArray< uint8_t > & inputbuffer = *(Ainputbuffers[tid]);
		libmaus2::autoarray::AutoArray< uint8_t > & copybuffer = *(Acopybuffers[tid]);
		libmaus2::autoarray::AutoArray< libmaus2::dazzler::align::OverlapDataInterface > & copypointers = *(Acopypointers[tid]);

		uint64_t ocur = ofrom;
		uint64_t blockid = 0;
		for ( ; ocur < oto; ++blockid )
		{
			uint64_t const space = inputbuffer.size();
			uint64_t const todo = oto - ocur;
			uint64_t const toread = std::min(space,todo);
			istr.read(reinterpret_cast<char * >(inputbuffer.begin()),toread);
			assert ( istr.gcount() == static_cast<int64_t>(toread) );

			parser.parseBlock(
				inputbuffer.begin(),
				inputbuffer.begin() + toread,
				libmaus2::dazzler::align::OverlapParser::overlapparser_do_split
			);

			libmaus2::dazzler::align::OverlapData & data = parser.getData();
			for ( uint64_t j = 0; j < data.size(); ++j )
			{
				std::pair<uint8_t const *, uint8_t const *> const P = data.getData(j);
				libmaus2::dazzler::align::OverlapDataInterface const I(P.first);

				uint64_t const score =
					static_cast<uint64_t>(ldexp((static_cast<double>(I.diffs()) / static_cast<double>(I.aepos() - I.abpos())),30));

				if ( RHO.f == lmaxinput )
				{
					if ( score > RHO.top().score )
						RHO.popvoid();
				}
				if ( RHO.f < lmaxinput )
				{
					RHO.push(OverlapEntry(score,blockid,j,P.second-P.first));
				}
			}

			ocur += toread;
		}
		assert ( parser.isIdle() );

		uint64_t const f = RHO.f;

		std::sort(RHO.H.begin(),RHO.H.begin() + RHO.f,OverlapEntryBlockComparator());

		uint64_t gs = 0;
		for ( uint64_t i = 0; i < RHO.f; ++i )
			gs += RHO.H[i].s;

		if ( gs > copybuffer.size() )
		{
			copybuffer = libmaus2::autoarray::AutoArray< uint8_t >();
			copybuffer = libmaus2::autoarray::AutoArray< uint8_t >(gs,false);
		}

		uint64_t o_copypointers = 0;
		uint8_t * copypointer = copybuffer.begin();
		while ( RHO.f && RHO.H[RHO.f-1].blockid == blockid-1 )
		{
			uint64_t const j = RHO.f-1;

			libmaus2::dazzler::align::OverlapData & data = parser.getData();
			std::pair<uint8_t const *, uint8_t const *> const P = data.getData(RHO.H[j].entryid);

			libmaus2::dazzler::align::OverlapDataInterface const OI(copypointer);
			std::copy(P.first,P.second,copypointer);
			copypointer += (P.second-P.first);
			copypointers.push(o_copypointers,OI);

			uint64_t const score =
				static_cast<uint64_t>(ldexp((static_cast<double>(OI.diffs()) / static_cast<double>(OI.aepos() - OI.abpos())),30));
			assert ( RHO.H[j].score == score );
			assert ( static_cast < ::std::ptrdiff_t >(RHO.H[j].s) == P.second-P.first );

			RHO.f--;
		}

		if ( RHO.f )
		{
			istr.clear();
			istr.seekg(ofrom);
			ocur = ofrom;
			blockid = 0;

			OverlapEntry const * opa = RHO.H.begin();
			OverlapEntry const * ope = opa + RHO.f;

			for ( ; ocur < oto; ++blockid )
			{
				OverlapEntry const * opc = opa;
				while ( opc != ope && opc->blockid == blockid )
					++opc;

				uint64_t const space = inputbuffer.size();
				uint64_t const todo = oto - ocur;
				uint64_t const toread = std::min(space,todo);
				istr.read(reinterpret_cast<char * >(inputbuffer.begin()),toread);
				assert ( istr.gcount() == static_cast<int64_t>(toread) );

				parser.parseBlock(
					inputbuffer.begin(),
					inputbuffer.begin() + toread,
					libmaus2::dazzler::align::OverlapParser::overlapparser_do_split
				);

				libmaus2::dazzler::align::OverlapData & data = parser.getData();

				for ( OverlapEntry const * opq = opa; opq != opc; ++opq )
				{
					OverlapEntry const & OE = *opq;
					uint64_t const entryid = OE.entryid;
					std::pair<uint8_t const *, uint8_t const *> const P = data.getData(entryid);

					libmaus2::dazzler::align::OverlapDataInterface const OI(copypointer);
					std::copy(P.first,P.second,copypointer);
					copypointer += (P.second-P.first);
					copypointers.push(o_copypointers,OI);
					assert ( static_cast<int64_t>(OE.s) == P.second-P.first );

					uint64_t const score =
						static_cast<uint64_t>(ldexp((static_cast<double>(OI.diffs()) / static_cast<double>(OI.aepos() - OI.abpos())),30));
					assert ( OE.score == score );
				}

				ocur += toread;
				opa = opc;
			}
			assert ( parser.isIdle() );
		}
		assert ( copypointer == copybuffer.begin() + gs );
		assert ( o_copypointers == f );
		RHO.clear();

		struct OverlapDataInterfacePosComparator
		{
			bool operator()(libmaus2::dazzler::align::OverlapDataInterface const & A, libmaus2::dazzler::align::OverlapDataInterface const & B) const
			{
				return A.abpos() < B.abpos();
			}
		};

		std::sort(
			copypointers.begin(),
			copypointers.begin() + o_copypointers,
			OverlapDataInterfacePosComparator()
		);

		#if 0
		for ( uint64_t i = 0; i < o_copypointers; ++i )
		{
			libmaus2::dazzler::align::OverlapDataInterface const & OI = copypointers[i];
			std::cerr << OI << std::endl;
		}
		#endif

		#if 0
		// score,id on RO
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp > & RH = *AAH[tid];
		RH.clear();
		libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

		if ( lmaxinput > RH.H.size() )
			RH.ensureSize(lmaxinput);
		if ( lmaxinput > RO.size() )
			RO.ensureSize(lmaxinput);


		libmaus2::dazzler::align::Overlap OVL;
		libmaus2::dazzler::align::AlignmentFileDecoder::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileAt(lasfn,z,z+1,*Adalindex[tid]));
		// libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
		uint64_t f = 0;
		while ( pdec->getNextOverlap(OVL) )
		{
			uint64_t const score =
				static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

			if ( RH.f == lmaxinput )
			{
				// dump alignment if error rate is higher than highest in heap
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
				RH.push(std::pair<uint64_t,uint64_t>(score,p));
				RO[p] = OVL;
			}

			// VOVL.push_back(OVL);
		}

		#if 0
		{
			libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "[V] " << z << " " << VOVL.size() << std::endl;
		}
		#endif

		std::sort(RO.begin(),RO.begin()+f,OverlapPosComparator());

		#if defined(CONS_HANDLE_SINGLE)
		if ( f && RO[0].aread == CONS_HANDLE_SINGLE )
		{
			for ( uint64_t i = 0; i < f; ++i )
			{
				std::cerr << RO[i] << std::endl;
				AW.put(RO[i]);
			}
		}
		#endif
		#endif

		std::ostringstream logostr;
		std::ostringstream outstr;

		//libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
		// DLOS << "[S] handling aread " << Vovl[low].aread << std::endl;
		if ( f )
			logostr << "[S] handling aread " << copypointers[0].aread() << " depth " << f << std::endl;

		// std::pair<unsigned char const *, unsigned char const *> const Q = inqual.getQualityForRead(Vovl[low].aread);

		#if 0
		for ( unsigned char const * p = Q.first; p != Q.second; ++p )
			std::cerr << "block " << p-Q.first << " qual " << static_cast<int>(*p) << std::endl;
		#endif


		if ( f )
		{
			try
			{
				#if defined(CONS_HANDLE_SINGLE)
				if ( RO[0].aread == CONS_HANDLE_SINGLE )
				#endif
				{
					#if 0
					DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
					DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);
					assert ( tid < ANP.size() );
					assert ( ANP[tid] );
					libmaus2::lcs::Aligner & al = *(ANP[tid]);
					#endif

					HandleContext & context = *(AHC[tid]);

					std::ostream & handleout = (numthreads > 1) ? outstr : std::cout;
					std::ostream & handleerr = (numthreads > 1) ? logostr : std::cerr;

					context(
						handleout,handleerr,
						copypointers.begin(),
						copypointers.begin() + f
						#if defined(HANDLE_DEBUG)
							,
						// for debugging
						*(Vbam[RO[0].aread]),
						Vbam,
						Vtext[Vbam[RO[0].aread]->getRefID()],
						*(AED[tid])
						#endif
					);
					#if 0
					context(
						handleout,handleerr,RO.begin(),RO.begin()+f
						#if defined(HANDLE_DEBUG)
							,
						// for debugging
						*(Vbam[RO[0].aread]),
						Vbam,
						Vtext[Vbam[RO[0].aread]->getRefID()],
						*(AED[tid])
						#endif
					);
					#endif

					#if 0
					handle(
						handleout,
						handleerr,
						maxalign,
						wellcounter,
						// actual data
						RO.begin(),RO.begin()+f,windowsize,advancesize,
						RDC,
						RDC2,
						#if 0
						Vreads,Rreads,
						#endif
						traceFreeList,
						al,
						tspace,
						DP,
						producefull,
						ADG,
						verbose,
						minwindowcov,
						eminrate,
						minlen,
						minfilterfreq,
						maxfilterfreq
						#if defined(HANDLE_DEBUG)
							,
						// for debugging
						*(Vbam[RO[0].aread]),
						Vbam,
						Vtext[Vbam[RO[0].aread]->getRefID()],
						*(AED[tid])
						#endif
					);
					#endif
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

		#if !defined(CONSENSUS_SERIAL)
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

		#if 0
		if ( z >= 20 )
			exit(0);
		#endif
	}

	#if defined(HANDLE_DEBUG)
	for ( uint64_t i = 1; i < AED.size(); ++i )
	{
		ed_type const & I = *(AED[i]);
		ed_type & O = *(AED[0]);

		for ( ed_type::const_iterator ita = I.begin(); ita != I.end(); ++ita )
			O [ ita->first ] += ita->second;
	}
	if ( AED.size() )
	{
		ed_type const & ED = *(AED[0]);

		uint64_t s = 0;
		for ( ed_type::const_iterator ita = ED.begin(); ita != ED.end(); ++ita )
			s += ita->second;

		uint64_t t = 0;
		for ( ed_type::const_iterator ita = ED.begin(); ita != ED.end(); ++ita )
		{
			t += ita->second;
			std::cerr << "[ED]\t" << ita->first << " " << static_cast<double>(t) / s << std::endl;
		}
	}
	#endif


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
		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 5 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] reads.las reads.dam reads.dam ref.fasta reads.bam\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return daccord(arg,arg[0],arg[1],arg[2],arg[3],arg[4]);
		}
		#else
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] reads.las reads.dam\n";
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

			if ( arg.size() == 2 )
				r = daccord(arg,arg[0],arg[1],arg[1]);
			else
				r = daccord(arg,arg[0],arg[1],arg[2]);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			return r;
		}
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
