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
#include <ComputeOffsetLikely.hpp>
#include <DebruijnGraph.hpp>

#include <config.h>

#include <libmaus2/aio/DebugLineOutputStream.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
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

// #define HANDLE_INDEL_ESTIMATE_DEBUG
// #define HANDLE_DEBUG
// #define HANDLE_DEEP_DEBUG
// #define CONSENSUS_SERIAL
// #define WINDOWALIGNAVOID
// #define CANDIDATE_AVOID_ALIGN
// #define CONS_HANDLE_SINGLE 97776

#if defined(HANDLE_INDEL_ESTIMATE_DEBUG) || defined(HANDLE_DEBUG)
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/lcs/NP.hpp>
#endif

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
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

/**
 * active sequence in multiple alignment
 **/
struct ActiveElement
{
	// base data for a sequence
	uint8_t const * ua;
	// base data for b sequence
	uint8_t const * ub;
	// current trace pointer between a and b sequence
	libmaus2::lcs::AlignmentTraceContainer::step_type const * ta;
	libmaus2::lcs::AlignmentTraceContainer::step_type const * te;
	// offset on b sequence
	uint64_t uboff;
	// error rate of alignment
	double erate;

	ActiveElement() {}
	ActiveElement(
		uint8_t const * rua,
		uint8_t const * rub,
		libmaus2::lcs::AlignmentTraceContainer::step_type const * rta,
		libmaus2::lcs::AlignmentTraceContainer::step_type const * rte,
		uint64_t const ruboff,
		double const rerate
	) : ua(rua), ub(rub), ta(rta), te(rte), uboff(ruboff), erate(rerate) {}
};

template<unsigned int k>
void handleIndelEstimate(
	std::ostream &
		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		err
		#endif
		,
	uint64_t const maxalign,
	#if 0
	// overlaps
	std::vector<libmaus2::dazzler::align::Overlap>::const_iterator ita,
	std::vector<libmaus2::dazzler::align::Overlap>::const_iterator ite,
	#endif
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
	std::map<uint64_t,uint64_t> & copyFreq,
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
	DebruijnGraph<k> DG;

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

		if ( MAo )
		{
			copyFreq[MAo-1]++;
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
				DG.filterFreq(2);
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
}

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
	#endif

	#if defined(HANDLE_DEEP_DEBUG)
	for ( uint64_t i = 0; i < MAo; ++i )
	{
		logstr << "[filterfreq=" << kmersize << "," << filterfreq << "] S " << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
	}
	#endif
}
#endif

void handle(
	std::ostream & out,
	std::ostream & err,
	uint64_t const maxalign,
	libmaus2::parallel::SynchronousCounter<uint64_t> & wellcounter,
	// overlaps
	#if 0
	std::vector<libmaus2::dazzler::align::Overlap>::const_iterator ita,
	std::vector<libmaus2::dazzler::align::Overlap>::const_iterator ite,
	#endif
	libmaus2::dazzler::align::Overlap const * const ita,
	libmaus2::dazzler::align::Overlap const * const ite,
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
	//unsigned int const k,
	OffsetLikely const & offsetLikely,
	// produce full sequence (even if uncorrected)
	bool const producefull,
	libmaus2::autoarray::AutoArray < DebruijnGraphInterface::unique_ptr_type > & ADG,
	int const verbose,
	uint64_t const minwindowcov,
	uint64_t const eminrate,
	uint64_t const minlen,
	int64_t const minfilterfreq,
	int64_t const maxfilterfreq
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
	// number of overlaps with A read
	uint64_t const nintv = ite-ita;

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
	std::map<uint64_t,trace_type::shared_ptr_type> Mtraces;

	// compute traces
	for ( uint64_t z = 0; z < nintv; ++z )
	{
		if ( ! windowstracepointaligned )
		{
			uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread));
			uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread));

			trace_type::shared_ptr_type Ptrace = traceFreeList.get();

			ita[z].computeTrace(ua,ub,tspace,*Ptrace,NP);

			Mtraces[z] = Ptrace;
		}
	}

	// end of interval heap
	libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > E(1024);
	// set of active traces/reads
	std::map < uint64_t, ActiveElement > activeset;
	libmaus2::autoarray::AutoArray < std::pair< uint8_t const *, uint64_t> > MA;

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

	#if 0
	uint64_t const ylimit = (maxaepos + advancesize >= windowsize) ? ((maxaepos + advancesize - windowsize) / advancesize) : 0;
	assert ( ylimit * advancesize + windowsize > maxaepos );
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
	};

	std::vector < PileElement > PV;
	int failcount = 0;
	int insufcount = 0;
	std::vector<uint64_t> Vsuf;

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

	uint64_t z = 0;
	for ( uint64_t y = 0; y < W.size(); ++y )
	{
		bool printlogstr = false;
		std::ostringstream logstr;

		// start of window on A
		uint64_t const astart = W[y].first; // y * advancesize;
		uint64_t const aend   = W[y].second; // astart + windowsize;
		assert ( aend <= maxaepos );

		logstr << std::string(20,'=') << " aread=" << aid << " window y=" << y << "/" << W.size() << " [" << astart << "," << aend << ")" << std::endl;

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
					trace_type::shared_ptr_type Ptrace = Mtraces.find(z)->second;
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

			std::map<uint64_t,uint64_t> Ldist;
			for ( uint64_t i = 0; i < MAo; ++i )
				Ldist [ MA [ i ] . second ] ++;


			std::vector < double > VVVV;
			for ( std::map<uint64_t,uint64_t>::const_iterator ita = Ldist.begin(); ita != Ldist.end(); ++ita )
			{
				while ( ! ( ita->first < VVVV.size() ) )
					VVVV.push_back(0);
				VVVV[ita->first] = ita->second-1;
			}

			for ( uint64_t i = 0; i < offsetLikely.DPnormSquare.size(); ++i )
			{
				double const v = offsetLikely.DPnormSquare[i].dotproduct(VVVV);
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

		if ( MAo >= minwindowcov ) // && y == 335 )
		{
			int64_t const elength = maxvprodindex+1;

			Vsuf.push_back(y);

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

			for ( uint64_t adgi = 0; adgi < ADG.size(); ++adgi )
			{
				DebruijnGraphInterface & DG = *(ADG[adgi]);

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
					DG.filterFreq(std::max(filterfreq,static_cast<int64_t>(1)));
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
						DG.setupAddHeap();
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
						bool const consok = DG.traverse(elength - 5,elength + 5,MA.begin(),MAo,16 /* maxfronpath */,16 /* maxfullpath */);
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

				#if 0
				std::cerr << std::string(20,'=') << " aread=" << aid << " window y=" << y << "/" << ylimit << " MAo=" << MAo << " maxvprodindex=" << maxvprodindex
					<< " " << std::string(consensus.first,consensus.second)
					<< std::endl;
				#endif

				#if defined(HANDLE_DEBUG)
				libmaus2::lcs::AlignmentStatistics AS;
				for ( uint64_t i = 0; i < MAo; ++i )
				{
					#if 0
					std::cerr << "slref.first=" << slref.first << " MA[i].second=" << MA[i].second << std::endl;
					std::cerr << std::string(refu,refu+slref.first) << std::endl;
					std::cerr << std::string(MA[i].first,MA[i].first+MA[i].second) << std::endl;
					#endif
					NP.align(refu,slref.first,MA[i].first,MA[i].second);
					AS += NP.getTraceContainer().getAlignmentStatistics();
				}

				// err << "minrate=" << minrate << " AS.getErrorRate()=" << AS.getErrorRate() << std::endl;

				if ( minrate > AS.getErrorRate() )
				{
					logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] consensus error " << minrate << " ref error " << AS.getErrorRate() << std::endl;

					for ( uint64_t i = 0; i < minDG->getNumCandidates(); ++i )
					{
						std::pair<uint8_t const *, uint8_t const *> consensus = minDG->getCandidate(i);
						logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] candidate " << std::string(consensus.first,consensus.second) << " " << minDG->getCandidateErrorU(MA.begin(),MAo,i) << " weight " << minDG->getCandidateWeight(i) << std::endl;
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
						PV.push_back(PE);
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
								PV.push_back(PE);
								break;
							}
							case libmaus2::lcs::BaseConstants::STEP_DEL:
							{
								PileElement PE(apos++,0,'D');
								PV.push_back(PE);
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

	std::sort(PV.begin(),PV.end());

	if ( producefull && (ita != ite) )
	{
		uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita->aread));

		uint64_t next = 0;
		uint64_t low = 0;
		std::vector < PileElement > NPV;

		while ( low < PV.size() )
		{
			// look for consecutive pile elements
			uint64_t high = low+1;
			while ( high < PV.size() && PV[low].apos == PV[high].apos )
				++high;

			// fill up before elements
			for ( ; static_cast<int64_t>(next) < PV[low].apos; ++next )
				NPV.push_back(PileElement(next,0,::tolower(ua[next])));

			// copy pile elements
			for ( uint64_t i = low; i < high; ++i )
				NPV.push_back(PV[i]);

			next = PV[low].apos+1;
			low = high;
		}

		// fill up until end
		uint64_t const rl = RC.getReadLength(ita->aread);
		for ( ; next < rl; ++next )
			NPV.push_back(PileElement(next,0,::tolower(ua[next])));

		PV = NPV;

		for ( uint64_t i = 1; i < PV.size(); ++i )
			assert ( PV[i-1].apos <= PV[i].apos );
	}

	std::vector< std::pair<uint64_t,uint64_t> > PVI;

	#if 0
	for ( uint64_t i = 0; i < PV.size(); ++i )
		std::cerr << PV[i].apos << "," << PV[i].apre << std::endl;
	#endif

	uint64_t il = 0;
	while ( il < PV.size() )
	{
		uint64_t ih = il+1;

		while (
			ih != PV.size()
			&&
			(PV[ih].apos - PV[ih-1].apos) <= 1
		)
			++ih;

		uint64_t const first = PV[il].apos;
		uint64_t const last = PV[ih-1].apos;

		if ( last-first >= 100 )
		{
			// std::cerr << "first=" << first << " last=" << last << std::endl;
			PVI.push_back(std::pair<uint64_t,uint64_t>(il,ih));
		}

		il = ih;
	}

	for ( uint64_t z = 0; z < PVI.size(); ++z )
	{
		std::pair<uint64_t,uint64_t> const P = PVI[z];

		uint64_t const first = PV[P.first].apos;
		uint64_t const last = PV[P.second-1].apos;

		std::vector<char> CO;

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
				CO.push_back(C[0].second);
		}

		std::reverse(CO.begin(),CO.end());

		std::string const SCO(CO.begin(),CO.end());

		if ( producefull || SCO.size() >= minlen )
		{
			out << '>' << (aid+1) << '/' << wellcounter++ << '/' << first << '_' << first + SCO.size() << " A=[" << first << "," << last << "]" << "\n";
			char const * zp = SCO.c_str();
			char const * ze = zp + SCO.size();
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

	#if 0
	bool Vsufsub = true;
	for ( uint64_t i = 1; i < Vsuf.size(); ++i )
		if ( Vsuf[i-1]+1 != Vsuf[i] )
			Vsufsub = false;
	#endif

	// err << "read id " << ita[0].aread << " " << AS << " failcount=" << failcount << " insufcount=" << insufcount << " unbroken=" << Vsufsub << " " << nnpres << " pvdense=" << pvdense << std::endl;

	// return traces
	for ( std::map<uint64_t,trace_type::shared_ptr_type>::const_iterator m_ita =
		Mtraces.begin(); m_ita != Mtraces.end(); ++m_ita )
		traceFreeList.put(m_ita->second);

	#if defined(HANDLE_DEBUG)
	// return BAM trace for ground truth
	traceFreeList.put(Pbamtrace);
	#endif

	// DG.printSize(err);

	err << "[V] read id " << ita[0].aread << " time " << handlertc << std::endl;
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
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	uint64_t const verbose = arg.uniqueArgPresent("V") ? arg.getUnsignedNumericArg<uint64_t>("V") : getDefaultVerbose();
	uint64_t const maxalign = arg.uniqueArgPresent("d") ? arg.getUnsignedNumericArg<uint64_t>("d") : getDefaultMaxAlign();
	uint64_t const keepeprof = arg.uniqueArgPresent("keepeprof") ? arg.getUnsignedNumericArg<uint64_t>("keepeprof") : 1;
	uint64_t const minwindowcov = arg.uniqueArgPresent("m") ? arg.getUnsignedNumericArg<uint64_t>("m") : getDefaultMinWindowCoverage();
	uint64_t const eminrate = arg.uniqueArgPresent("e") ? arg.getUnsignedNumericArg<uint64_t>("e") : getDefaultMinWindowError();
	uint64_t const minlen = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultMinLen();

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
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
	PDB->computeTrimVector();

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

	//libmaus2::fastx::ReadContainer const RC(readsfn);
	// std::cerr << "opening " << dbptr->fn << std::endl;

	// load read/read alignments
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(lasfn);

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

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp > & RH = *AAH[tid];
			RH.clear();
			libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

			libmaus2::lcs::Aligner & al = *(ANP[tid]);

			// std::vector<libmaus2::dazzler::align::Overlap> VOVL;
			libmaus2::dazzler::align::Overlap OVL;

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
			uint64_t f = 0;
			while ( pdec->getNextOverlap(OVL) )
			{
				uint64_t const score =
					static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

				if ( RH.full() )
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
				std::map<uint64_t,uint64_t> & copyFreq = *(Vcopyfreq[tid]);

				uint64_t lusable = 0, lunusable = 0;
				DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
				DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);

				// if ( VOVL[0].aread == 137 )
				{
					#if 0
					char const * f = RDC.getForwardRead(VOVL[0].aread);
					uint64_t const l = RDC.getReadLength(VOVL[0].aread);
					std::cerr << std::string(f,f+l) << std::endl;
					#endif

					handleIndelEstimate<8>(
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
						copyFreq,
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
				GASlock.unlock();
			}
		}

		std::cerr << "usable=" << usable << " unusable=" << unusable << std::endl;

		#if 0
		copy_freq_type copyFreq;
		for ( uint64_t i = 0; i < Vcopyfreq.size(); ++i )
		{
			copy_freq_type & lcopyfreq = *(Vcopyfreq[i]);
			for ( copy_freq_type::const_iterator ita = lcopyfreq.begin(); ita != lcopyfreq.end(); ++ita )
				copyFreq[ita->first] += ita->second;
		}

		for ( copy_freq_type::const_iterator ita = copyFreq.begin(); ita != copyFreq.end(); ++ita )
		{
			std::cerr << "copyfreq\t" << ita->first << "\t" << ita->second << std::endl;
		}
		#endif

		#if defined(HANDLE_INDEL_ESTIMATE_DEBUG)
		std::cerr << "TAS=" << TAS << std::endl;
		#endif

		GAS.serialise(eproffntmp);
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

	std::vector<uint64_t> DGKV;
	for ( uint64_t k = kmersizelow; k <= kmersizehigh; ++k )
		DGKV.push_back(k);
	std::sort(DGKV.begin(),DGKV.end(),std::greater<uint64_t>());

	std::cerr << "[V] using kmer range [" << kmersizelow << "," << kmersizehigh << "]" << std::endl;

	typedef DebruijnGraphInterface::unique_ptr_type graph_ptr_type;
	typedef libmaus2::autoarray::AutoArray < graph_ptr_type > graph_array_type;
	typedef graph_array_type::unique_ptr_type graph_array_ptr_type;
	libmaus2::autoarray::AutoArray < graph_array_ptr_type > AADG(numthreads);

	for ( uint64_t i = 0; i < AADG.size(); ++i )
	{
		graph_array_ptr_type tptr(new graph_array_type(DGKV.size()));
		AADG[i] = UNIQUE_PTR_MOVE(tptr);
		graph_array_type & ADG = *(AADG[i]);
		assert ( ADG.size() == DGKV.size() );

		for ( uint64_t j = 0; j < ADG.size(); ++j )
		{
			switch ( DGKV[j] )
			{
				case 3:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<3>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 4:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<4>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 5:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<5>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 6:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<6>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 7:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<7>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 8:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<8>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 9:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<9>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 10:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<10>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 11:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<11>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 12:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<12>);
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				default:
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "k-mer size " << DGKV[j] << " is not compiled in" << std::endl;
					lme.finish();
					throw lme;
					break;
				}
			}
		}
	}

	std::cerr << "[V] graph structures set up" << std::endl;

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

		graph_array_type & ADG = *(AADG[tid]);

		// score,id on RO
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t>, PairFirstGreaterComp > & RH = *AAH[tid];
		RH.clear();
		libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

		libmaus2::dazzler::align::Overlap OVL;
		libmaus2::dazzler::align::AlignmentFileDecoder::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileAt(lasfn,z,z+1,*Adalindex[tid]));
		// libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
		uint64_t f = 0;
		while ( pdec->getNextOverlap(OVL) )
		{
			uint64_t const score =
				static_cast<uint64_t>(ldexp((static_cast<double>(OVL.path.diffs) / static_cast<double>(OVL.path.aepos - OVL.path.abpos)),30));

			if ( RH.full() )
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

		std::ostringstream logostr;
		std::ostringstream outstr;

		//libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
		// DLOS << "[S] handling aread " << Vovl[low].aread << std::endl;
		if ( f )
			logostr << "[S] handling aread " << RO[0].aread << " depth " << f << std::endl;

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
					DecodedReadContainer RDC(readDataFreeList,readDecoderFreeList);
					DecodedReadContainer RDC2(readDataFreeList,readDecoderFreeList2);
					assert ( tid < ANP.size() );
					assert ( ANP[tid] );
					libmaus2::lcs::Aligner & al = *(ANP[tid]);

					handle(
						#if defined(CONSENSUS_SERIAL)
						std::cout,
						std::cerr,
						#else
						outstr,
						logostr,
						#endif
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
