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
	return 4;
}

static double getDefaultMinWindowError()
{
	return std::numeric_limits<double>::max();
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

/**
 * dot product
 **/
struct DotProduct
{
	// first significant index
	uint64_t firstsign;
	// vector of coefficients
	std::vector<double> V;

	DotProduct() {}
	DotProduct(uint64_t const rfirstsign, std::vector<double> const & rV)
	: firstsign(rfirstsign), V(rV)
	{

	}

	std::ostream & printData(std::ostream & out) const
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i] )
				out << i+firstsign << "\t" << V[i] << "\n";
		return out;
	}

	uint64_t size() const
	{
		return firstsign + V.size();
	}

	double operator[](uint64_t const i) const
	{
		if ( i < firstsign )
			return 0.0;

		uint64_t j = i - firstsign;

		if ( j < V.size() )
			return V[j];
		else
			return 0.0;
	}

	void normaliseValue(uint64_t const i, double const div)
	{
		if ( i >= firstsign )
		{
			uint64_t const j = i-firstsign;

			if ( j < V.size() )
				V[j] /= div;
		}
	}

	// compute product
	double dotproduct(std::vector < double > const & O) const
	{
		double s = 0;

		// iterate over length of V
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			// compute corresponding index j on O
			uint64_t const j = firstsign + i;

			// if j is in range for O
			if ( j < O.size() )
				s += V[i] * O[j];
			// j is too large, stop
			else
				break;
		}

		return s;
	}

	// normalise the vector (make dot product between the vector and itself 1)
	void normalise()
	{
		double s = 0.0;
		for ( uint64_t i = 0; i < V.size(); ++i )
			s += V[i]*V[i];
		double const c = std::sqrt(1.0/s);
		for ( uint64_t i = 0; i < V.size(); ++i )
			V[i] *= c;
	}
};

/**
 * output operator for dot product
 **/
std::ostream & operator<<(std::ostream & out, DotProduct const & DP)
{
	out << "DotProduct(firstsign=" << DP.firstsign << ",";

	for ( uint64_t j = 0; j < DP.V.size(); ++j )
		out << DP.V[j] << ";";
	out << ")";

	return out;
}

struct OffsetLikely
{
	std::vector<DotProduct> DP;
	std::vector < double > dsum;
	std::vector<DotProduct> DPnorm;
	// read position [] is relevant for reference positions in given interval
	std::vector< std::pair<uint64_t,uint64_t> > Vsupport;
	std::vector<DotProduct> DPnormSquare;

	// get supporting ref position lower bound for read position i
	uint64_t getSupportLow(int64_t const i) const
	{
		return expect_true(i < static_cast<int64_t>(Vsupport.size())) ? Vsupport[i].first : DPnorm.size();
	}

	uint64_t getSupportHigh(int64_t const i) const
	{
		return expect_true(i < static_cast<int64_t>(Vsupport.size())) ? Vsupport[i].second : DPnorm.size();
	}

	uint64_t size() const
	{
		return DP.size();
	}

	void push_back(DotProduct const & D)
	{
		DP.push_back(D);
	}

	DotProduct const & operator[](uint64_t const i) const
	{
		return DP[i];
	}

	void setup()
	{
		dsum.resize(0);

		// maximum size of any dot product (max read position supported + 1)
		uint64_t maxsize = 0;
		for ( uint64_t i = 0; i < size(); ++i )
			maxsize = std::max(maxsize,DP[i].size());

		for ( uint64_t i = 0; i < maxsize; ++i )
		{
			double sum = 0.0;

			for ( uint64_t j = 0; j < size(); ++j )
				sum += (*this)[j][i];

			dsum.push_back(sum);
		}

		DPnorm = DP;
		for ( uint64_t i = 0; i < DPnorm.size(); ++i )
			for ( uint64_t j = 0; j < maxsize; ++j )
				DPnorm[i].normaliseValue(j,dsum[j]);

		uint64_t j = 0, k = 0;
		for ( uint64_t i = 0; i < maxsize; ++i )
		{
			// update lower end
			while ( j < DPnorm.size() && i >= DPnorm[j].firstsign + DPnorm[j].V.size() )
				++j;
			// update upper end
			while ( k < DPnorm.size() && DPnorm[k].firstsign <= i )
				++k;

			Vsupport.push_back(std::pair<uint64_t,uint64_t>(j,k));
		}

		DPnormSquare = DP;
		for ( uint64_t i = 0; i < DPnormSquare.size(); ++i )
			DPnormSquare[i].normalise();
	}

	void plotData(std::string const & prefix) const
	{
		std::vector<std::string> Vfn;
		std::string const gplfn = prefix + ".gpl";
		libmaus2::util::TempFileRemovalContainer::addTempFile(gplfn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type POSI(
			new libmaus2::aio::OutputStreamInstance(gplfn)
		);

		*POSI << "set terminal pdf\n";
		*POSI << "set xlabel \"position in read\"\n";
		*POSI << "set ylabel \"probability\"\n";
		*POSI << "plot [0 to 55] ";
		for ( uint64_t i = 0; i < DP.size(); ++i )
		{
			std::ostringstream fnostr;
			fnostr << prefix << "_" << i << ".gpl";
			std::string const fn = fnostr.str();
			Vfn.push_back(fn);
			libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
			libmaus2::aio::OutputStreamInstance OSI(fn);
			DP[i].printData(OSI);

			if ( i > 0 )
				*POSI << ",";

			*POSI << "\"" << fn << "\" smooth bezier title \"\"";
		}
		*POSI << "\n";
		POSI->flush();
		POSI.reset();

		std::ostringstream comstr;
		comstr << "gnuplot <" << gplfn << " >" << prefix << ".pdf";
		std::string const com = comstr.str();

		int const r = system(com.c_str());
		if ( r != EXIT_SUCCESS )
		{

		}
	}

	void plotDataSingle(std::string const & prefix, uint64_t const z) const
	{
		std::vector<std::string> Vfn;
		std::string const gplfn = prefix + ".gpl";
		libmaus2::util::TempFileRemovalContainer::addTempFile(gplfn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type POSI(
			new libmaus2::aio::OutputStreamInstance(gplfn)
		);

		*POSI << "set terminal pdf\n";
		*POSI << "set xlabel \"position in read\"\n";
		*POSI << "set ylabel \"probability\"\n";
		*POSI << "plot [0 to 30] ";
		bool first = true;
		for ( uint64_t i = z; i < DP.size() && i < z+1; ++i )
		{
			std::ostringstream fnostr;
			fnostr << prefix << "_" << i << ".gpl";
			std::string const fn = fnostr.str();
			Vfn.push_back(fn);
			libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
			libmaus2::aio::OutputStreamInstance OSI(fn);
			DP[i].printData(OSI);

			if ( ! first )
				*POSI << ",";

			*POSI << "\"" << fn << "\" smooth bezier title \"position " << z << " in true sequence\"";
			first = false;
		}
		*POSI << "\n";
		POSI->flush();
		POSI.reset();

		std::ostringstream comstr;
		comstr << "gnuplot <" << gplfn << " >" << prefix << ".pdf";
		std::string const com = comstr.str();

		int const r = system(com.c_str());
		if ( r != EXIT_SUCCESS )
		{

		}
	}
};

struct Node
{
	uint64_t v;
	uint64_t spo;
	uint64_t freq;
	uint64_t numsucc;
	uint64_t numsuccactive;
	uint64_t feaspos;
	uint64_t cfeaspos;
	uint64_t numfeaspos;
	uint64_t numcfeaspos;
	uint64_t pfostart;
	uint64_t pfosize;
	uint64_t cpfostart;
	uint64_t cpfosize;
	uint64_t plow;
	uint64_t phigh;
	uint64_t cplow;
	uint64_t cphigh;

	Node() {}
	Node(uint64_t const rv) : v(rv) {}
	Node(
		uint64_t const rv,
		uint64_t const rspo,
		uint64_t const rfreq,
		uint64_t const rpfostart,
		uint64_t const rpfosize,
		uint64_t const rcpfostart,
		uint64_t const rcpfosize,
		uint64_t const rplow,
		uint64_t const rphigh,
		uint64_t const rcplow,
		uint64_t const rcphigh
	)
	: v(rv), spo(rspo), freq(rfreq), numsucc(0), numsuccactive(0), feaspos(0), cfeaspos(0), numfeaspos(0), numcfeaspos(0),
	  pfostart(rpfostart), pfosize(rpfosize), cpfostart(rcpfostart), cpfosize(rcpfosize), plow(rplow), phigh(rphigh), cplow(rcplow), cphigh(rcphigh) {}
};

std::ostream & operator<<(std::ostream & out, Node const & N)
{
	return out << "Node(v=" << N.v
		<< ",spo=" << N.spo
		<< ",freq=" << N.freq
		<< ",numsucc=" << N.numsucc
		<< ",numsuccactive=" << N.numsuccactive
		<< ",feaspos=" << N.feaspos
		<< ",numfeaspos=" << N.numfeaspos
		<< ",pfostart=" << N.pfostart
		<< ",pfosize=" << N.pfosize
		<< ",plow=" << N.plow
		<< ",phigh=" << N.phigh
		<< ")";
}

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

	uint32_t hash() const
	{
		uint64_t A[2] = { first, ext };
		return libmaus2::hashing::EvaHash::hash2(reinterpret_cast<uint32_t const *>(&A[0]),4);
	}

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

struct Links
{
	public:
	uint64_t A[4];
	uint64_t p;

	public:
	void reset()
	{
		p = 0;
	}

	void push(uint64_t const sym, uint64_t const freq)
	{
		if ( freq )
			A[p++] = (freq << 8) | sym;
	}

	void setSize(uint64_t const rp)
	{
		p = rp;
	}

	Links()
	{
		reset();
	}

	void sort()
	{
		if ( p <= 1 )
			return;

		std::sort(&A[0],&A[p],std::greater<uint64_t>());
	}

	uint64_t size() const
	{
		return p;
	}

	uint64_t getFreq(uint64_t const i) const
	{
		return A[i] >> 8;
	}

	uint64_t getSym(uint64_t const i) const
	{
		return A[i] & 0xFF;
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

struct DebruijnGraphBase
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

struct DebruijnGraphInterface
{
	typedef DebruijnGraphInterface this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual ~DebruijnGraphInterface() {}

	virtual void setup(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o) = 0;
	virtual void filterFreq(uint64_t const f) = 0;
	// static double getDefaultComputeFeasibleKmerPositionsThres() { return 1e-3; }
	virtual void computeFeasibleKmerPositions(OffsetLikely const & offsetLikely, double const thres) = 0;
	virtual void getLevelSuccessors(unsigned int const s) = 0;
	virtual void setupNodes() = 0;
	virtual void setupAddHeap() = 0;
	virtual bool addNextFromHeap(std::ostream * logstr) = 0;
	virtual bool traverse(int64_t const lmin, int64_t const lmax, std::pair< uint8_t const *, uint64_t> const * MA, uint64_t const MAo,
		uint64_t const maxfrontpath, uint64_t const maxfullpath) = 0;
	virtual std::pair<uint64_t,double> checkCandidates(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o) = 0;
	virtual std::pair<uint8_t const *, uint8_t const *> getCandidate(uint64_t const i) const = 0;
	virtual uint64_t getNumCandidates() const = 0;
	virtual double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	) = 0;
	virtual double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o,
		uint64_t const id
	) = 0;
	virtual double getCandidateWeight(uint64_t const i) const = 0;
	virtual Node const * getNodeVirtual(uint64_t const v) const = 0;
	virtual bool isEdgeActiveVirtual(uint64_t const from, uint64_t const to) const = 0;
	virtual void getActiveSuccessorsVirtual(uint64_t const v, Links & L) const = 0;
	virtual void getSuccessorsVirtual(uint64_t const v, Links & L) const = 0;
	virtual uint64_t getKmerSize() const = 0;
	virtual std::string printNode(Node const & S) const = 0;
	virtual std::ostream & print(std::ostream & out) const = 0;
};

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

		CH.clear();
		while ( ! CDH.empty() )
			CH.pushBump(CDH.pop());
		while ( ! CH.empty() )
		{
			ConsensusCandidate CC = CH.pop();

			#if ! defined(CANDIDATE_AVOID_ALIGN)
			double const canderr = getSimpleCandidateError(MA,MAo,
				std::pair<uint8_t const *, uint8_t const *>(Acons.begin()+CC.o,Acons.begin()+CC.o+CC.l)
			);

			CC.error = canderr;
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

	double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o,
		uint64_t const id
	)
	{
		std::pair<uint8_t const *, uint8_t const *> Ucand = getCandidate(id);
		return getCandidateError(I,o,Ucand);
	}

	// setup using list of (string,length) pairs
	std::pair<uint64_t,double> checkCandidates(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o)
	{
		if ( getNumCandidates() )
			return std::pair<uint64_t,double>(0,getCandidateError(I,o,0));
		else
			return std::pair<uint64_t,double>(0,std::numeric_limits<double>::max());
	}
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
	for ( uint64_t z = 0; z < nintv; ++z )
	{
		// get trace from free list
		trace_type::shared_ptr_type Ptrace = traceFreeList.get();

		// update maximum aepos
		if ( ita[z].path.aepos > static_cast<int64_t>(maxaepos) )
			maxaepos = ita[z].path.aepos;

		// get pointers to data
		uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread));
		uint8_t const * ub = reinterpret_cast<uint8_t const *>(ita[z].isInverse() ? RC2.getReverseComplementRead(ita[z].bread) : RC2.getForwardRead(ita[z].bread));

		// compute the trace operations
		ita[z].computeTrace(ua,ub,tspace,*Ptrace,NP);

		// store trace object
		Mtraces[z] = Ptrace;
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

		if ( MAo >= 5 )
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
	double const eminrate,
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

	#if ! defined(WINDOWALIGNAVOID)
	bool const windowstracepointaligned = ((windowsize % tspace) == 0) && ((advancesize % tspace) == 0);
	#else
	bool const windowstracepointaligned = false;
	#endif

	// number of overlaps with A read
	uint64_t const nintv = ite-ita;

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
	uint64_t maxaepos = 0;
	for ( uint64_t z = 0; z < nintv; ++z )
	{
		if ( ita[z].path.aepos > static_cast<int64_t>(maxaepos) )
			maxaepos = ita[z].path.aepos;

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

	uint64_t const ylimit = (maxaepos + advancesize >= windowsize) ? ((maxaepos + advancesize - windowsize) / advancesize) : 0;
	assert ( ylimit * advancesize + windowsize > maxaepos );

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

	uint64_t z = 0;
	for ( uint64_t y = 0; y < ylimit; ++y )
	{
		bool printlogstr = false;
		std::ostringstream logstr;

		// start of window on A
		uint64_t const astart = y * advancesize;
		uint64_t const aend = astart + windowsize;
		assert ( aend <= maxaepos );

		logstr << std::string(20,'=') << " aread=" << aid << " window y=" << y << "/" << ylimit << " [" << astart << "," << aend << ")" << std::endl;

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

				// sanity check for aoff
				assert ( (ita[z].path.abpos + aoff) % advancesize == 0 );

				uint8_t const * ua = reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)) + astart;
				assert ( (!windowstracepointaligned) || ((((ua - reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread))) % tspace) == 0)));
				assert (((ua - reinterpret_cast<uint8_t const *>(RC.getForwardRead(ita[z].aread)))) == static_cast< ::std::ptrdiff_t > (y*advancesize));

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
		while ( (! (E.empty())) && E.top().first <= aend )
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
				std::pair<uint64_t,uint64_t> const advadv = libmaus2::lcs::AlignmentTraceContainer::advanceA(AE.ta,AE.te,advancesize);
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
			AE.ua += advancesize;
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
			double minrate = eminrate;
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
							std::pair<uint64_t,double> const MR = DG.checkCandidates(MA.begin(),MAo);

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
						logstr << "[filterfreq=" << minDG->getKmerSize() << "," << filterfreq << "] candidate " << std::string(consensus.first,consensus.second) << " " << minDG->getCandidateError(MA.begin(),MAo,i) << " weight " << minDG->getCandidateWeight(i) << std::endl;
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
		uint64_t const refadv = frontsoftclip >= advancesize ? 0 : advancesize - frontsoftclip;
		// advancement in trace operations
		std::pair<uint64_t,uint64_t> const advrefadv = libmaus2::lcs::AlignmentTraceContainer::advanceB(rta,rte,refadv);
		// advancement in string length
		std::pair<uint64_t,uint64_t> const slrefadv = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(rta,rta+advrefadv.second);
		// update front soft clipping
		frontsoftclip -= std::min(frontsoftclip,advancesize);
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

OffsetLikely computeOffsetLikely(
	// maximum length considered
	uint64_t const maxl,
	// I prob
	double const p_i,
	// D prob
	double const p_d
)
{
	// std::vector < DotProduct > VD;
	OffsetLikely VD;

	// compute probability bins for number of insertions at/before a single position
	// double const p_i = 0.1;
	double const q_i = 1.0 - p_i;
	double f_i = q_i;

	std::vector < double > P_I;
	while ( f_i >= 1e-7 )
	{
		P_I.push_back(f_i);
		f_i *= p_i;
	}

	// convolution accumulation vector for insertions
	std::vector < double > C_I(1,1.0);

	// double const p_d = 0.02;

	// iterate over length
	for ( uint64_t l = 0; l <= maxl; ++l )
	{
		// std::cerr << std::string(40,'*') << l << std::string(40,'*') << std::endl;
		//
		// add another round of insertions
		#if 1
		C_I = libmaus2::math::Convolution::convolutionFFTRef(C_I,P_I);
		#else
		C_I = libmaus2::math::Convolution::convolutionFFT(C_I,P_I);
		#endif

		// likelihood of deletion count
		std::vector < libmaus2::math::GmpFloat > V_D = libmaus2::math::Binom::binomVector(p_d, l, 512);
		assert ( V_D.size() );

		// insertion vector with V_D.size()-1 leading zeroes
		std::vector < double > V_I (V_D.size()-1+C_I.size());
		std::copy(
			C_I.begin(),C_I.end(),
			V_I.begin() + (V_D.size()-1)
		);

		// reverse deletion vector
		std::reverse(V_D.begin(),V_D.end());

		#if 0
		for ( uint64_t j = 0; j < V_D.size(); ++j )
			if ( static_cast<double>(V_D[j]) >= 1e-5 )
				std::cerr << "D j=" << j << " V_D=" << V_D[j] << std::endl;
		for ( uint64_t j = 0; j < V_I.size(); ++j )
			if ( V_I[j] >= 1e-5 )
				std::cerr << "I j=" << j << " V_I=" << V_I[j] << std::endl;
		#endif

		// compute convolution
		std::vector < double > F_I = libmaus2::math::Convolution::convolutionFFTRef(V_D,V_I);
		// first significant value found
		bool signfound = false;
		// index of first significant value
		int64_t firstsign = std::numeric_limits<int64_t>::min();
		// value vector
		std::vector < double > VP;

		for ( uint64_t j = 0; j < F_I.size(); ++j )
			if ( F_I[j] >= 1e-5 )
			{
				if ( ! signfound )
				{
					signfound = true;
					firstsign = static_cast<int64_t>(j)-static_cast<int64_t>(l);
				}

				uint64_t const offset =
					static_cast<int64_t>(j)-static_cast<int64_t>(l) - firstsign;

				while ( ! (offset < VP.size()) )
					VP.push_back(0);
				VP[offset] = F_I[j];

				// std::cerr << "F j=" << static_cast<int64_t>(j)-static_cast<int64_t>(2*l) << " abs=" << j-l << " off=" << offset << " F_I=" << F_I[j] << std::endl;
			}
		assert ( firstsign >= 0 );

		VD.push_back(DotProduct(firstsign,VP));
		// std::cerr << "first sign " << firstsign << std::endl;

		// add another round of insertions
		//C_I = libmaus2::math::Convolution::convolutionFFT(C_I,P_I);
	}

	#if 0
	for ( uint64_t i = 0; i < VD.size(); ++i )
		VD[i].normalise();
	#endif

	VD.setup();

	return VD;
}

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
	double const eminrate = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultMinWindowError();
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
	libmaus2::autoarray::AutoArray< libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >::unique_ptr_type > AAH(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >::unique_ptr_type tptr(
			new libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> >(maxinput)
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

			libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & RH = *AAH[tid];
			RH.clear();
			libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > & RO = AAOVL[tid];

			libmaus2::lcs::Aligner & al = *(ANP[tid]);

			// std::vector<libmaus2::dazzler::align::Overlap> VOVL;
			libmaus2::dazzler::align::Overlap OVL;

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type pdec(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(lasfn,z,z+1,memlasindexname));
			uint64_t f = 0;
			while ( pdec->getNextOverlap(OVL) )
			{
				uint64_t const score = OVL.path.aepos - OVL.path.abpos;

				if ( RH.full() )
				{
					if ( score <= RH.top().first )
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

				DLOS << "[S] handling aread " << RO[0].aread << std::endl;

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

		libmaus2::util::FiniteSizeHeap< std::pair<uint64_t,uint64_t> > & RH = *AAH[tid];
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
				if ( score <= RH.top().first )
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
			logostr << "[S] handling aread " << RO[0].aread << std::endl;

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
