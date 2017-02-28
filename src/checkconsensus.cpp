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

#include <complex>
#include <config.h>
#include <libmaus2/aio/ConcatInputStream.hpp>
#include <libmaus2/aio/DebugLineOutputStream.hpp>
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/aio/OutputStreamFactoryContainer.hpp>
#include <libmaus2/aio/TempFileArray.hpp>
#include <libmaus2/autoarray/AutoArray2d.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/BamNumericalIndexDecoder.hpp>
#include <libmaus2/bambam/BamNumericalIndexGenerator.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/db/InqualContainer.hpp>
#include <libmaus2/fastx/KmerFreqInWindowStats.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/graph/POGraph.hpp>
#include <libmaus2/hashing/hash.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NPL.hpp>
#include <libmaus2/lcs/SimdX86GlobalAlignmentY256_8.hpp>
#include <libmaus2/lcs/SuffixArrayLCS.hpp>
#include <libmaus2/lcs/SuffixPrefix.hpp>
#include <libmaus2/math/Convolution.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/math/binom.hpp>
#include <libmaus2/math/ipow.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/rank/ERank222B.hpp>
#include <libmaus2/rmq/QuickDynamicRMQ.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/sv/sv.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/BorderArray.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>
#include <libmaus2/util/ToUpperTable.hpp>
#include <libmaus2/wavelet/WaveletTree.hpp>

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

static uint64_t getDefaultNumPacks()
{
	return 1;
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
	optionMap . push_back ( std::pair < std::string, std::string >("p", formatRHS("number of packages for batchlist",getDefaultNumPacks())));

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

void patternToUpper(libmaus2::fastx::Pattern & pattern)
{
	for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
		pattern.spattern[i] = ::toupper(pattern.spattern[i]);

	pattern.pattern = pattern.spattern.c_str();
}

// load reference text from a FastA file
static std::vector<std::string> loadTextVector(std::string const & textfn, std::vector<std::string> & names)
{
	libmaus2::aio::InputStreamInstance ISI(textfn);
	libmaus2::fastx::StreamFastAReaderWrapper SFARW(ISI);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
	std::vector<std::string> Vtext;

	while ( SFARW.getNextPatternUnlocked(pattern) )
	{
		patternToUpper(pattern);
		Vtext.push_back(pattern.spattern);
		names.push_back(pattern.getShortStringId());
	}

	return Vtext;
}

std::string printDoubleFixed(double const v, uint64_t const w)
{
	std::ostringstream ostr;
	ostr << std::setprecision(w) << std::fixed << v;
	return ostr.str();
}

uint64_t getReadsHigh(std::string const & consfnindex)
{
	assert ( libmaus2::util::GetFileSize::getFileSize(consfnindex) % sizeof(libmaus2::fastx::FastInterval) == 0 );

	uint64_t const numentries = libmaus2::util::GetFileSize::getFileSize(consfnindex) / sizeof(libmaus2::fastx::FastInterval);
	uint64_t numreads = 0;
	if ( numentries )
	{
		libmaus2::aio::InputStreamInstance ISI(consfnindex);
		ISI.seekg((numentries-1) * sizeof(libmaus2::fastx::FastInterval));
		libmaus2::fastx::FastInterval FI(ISI);
		numreads = FI.high;
	}

	return numreads;
}

uint64_t getReadsLow(std::string const & consfnindex)
{
	assert ( libmaus2::util::GetFileSize::getFileSize(consfnindex) % sizeof(libmaus2::fastx::FastInterval) == 0 );

	uint64_t const numentries = libmaus2::util::GetFileSize::getFileSize(consfnindex) / sizeof(libmaus2::fastx::FastInterval);
	uint64_t numreads = 0;
	if ( numentries )
	{
		libmaus2::aio::InputStreamInstance ISI(consfnindex);
		libmaus2::fastx::FastInterval FI(ISI);
		numreads = FI.low;
	}

	return numreads;
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

static uint64_t getBamId(std::string const & name)
{
	std::deque<std::string> const tokens = libmaus2::util::stringFunctions::tokenize<std::string>(name,std::string("/"));
	assert ( tokens.size() >= 2 );
	std::string sid = tokens[1];
	std::istringstream istr(sid);
	uint64_t id;
	istr >> id;
	assert ( istr && istr.peek() == std::istream::traits_type::eof() );
	return id;
}

struct ConsensusEntry
{
	uint64_t id;
	libmaus2::lcs::AlignmentStatistics AS;
	libmaus2::math::IntegerInterval<int64_t> IN;
	uint64_t refid;
	uint64_t refpos;
	uint64_t reflen;

	ConsensusEntry()
	{

	}

	ConsensusEntry(std::istream & in)
	: id(libmaus2::util::NumberSerialisation::deserialiseNumber(in)), AS(in)
	{
		IN.from = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		IN.to = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		refid = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		refpos = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		reflen = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
	}

	ConsensusEntry(
		uint64_t const rid,
		libmaus2::lcs::AlignmentStatistics const & RAS,
		libmaus2::math::IntegerInterval<int64_t> const & RIN,
		uint64_t const rrefid,
		uint64_t const rrefpos,
		uint64_t const rreflen
	)
	: id(rid), AS(RAS), IN(RIN), refid(rrefid), refpos(rrefpos), reflen(rreflen) {}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,id);
		AS.serialise(out);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,IN.from);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,IN.to);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,refid);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,refpos);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,reflen);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = ConsensusEntry(in);
		return in;
	}

	bool operator<(ConsensusEntry const & O) const
	{
		return id < O.id;
	}
};

struct EmapEntry
{
	double v;
	uint64_t i;

	EmapEntry()
	{

	}

	EmapEntry(std::istream & in)
	: v(libmaus2::util::NumberSerialisation::deserialiseDouble(in)),
	  i(libmaus2::util::NumberSerialisation::deserialiseNumber(in)) {}

	EmapEntry(double const rv, uint64_t const ri) : v(rv), i(ri) {}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseDouble(out,v);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,i);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = EmapEntry(in);
		return in;
	}

	bool operator<(EmapEntry const & E) const
	{
		if ( v != E.v )
			return v < E.v;
		else
			return i < E.i;
	}

	bool operator>(EmapEntry const & E) const
	{
		if ( v != E.v )
			return v > E.v;
		else
			return i > E.i;
	}
};

void mergeEmap(std::string const & fn)
{
	std::string const outfn = fn + ".merge";
	libmaus2::util::TempFileRemovalContainer::addTempFile(outfn);

	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(outfn));
	libmaus2::aio::InputStreamInstance::unique_ptr_type ISI(new libmaus2::aio::InputStreamInstance(fn));

	uint64_t c = 0;
	double prev = std::numeric_limits<double>::max();
	EmapEntry E;

	while ( *ISI && ISI->peek() != std::istream::traits_type::eof() )
	{
		E.deserialise(*ISI);

		if ( E.v != prev && c )
		{
			EmapEntry(prev,c).serialise(*OSI);
			c = 0;
		}

		c += E.i;
		prev = E.v;
	}

	if ( c )
		EmapEntry(prev,c).serialise(*OSI);

	OSI->flush();
	OSI.reset();
	ISI.reset();

	libmaus2::aio::OutputStreamFactoryContainer::rename(outfn,fn);
}

void reverseEmap(std::vector<std::string> const & in, std::string const & out, std::string const & tmp)
{
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmp);
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<EmapEntry,std::greater<EmapEntry> >::unique_ptr_type SSBOF(
		new libmaus2::sorting::SerialisingSortingBufferedOutputFile<EmapEntry,std::greater<EmapEntry>>(tmp));
	libmaus2::aio::ConcatInputStream::unique_ptr_type ISI(new libmaus2::aio::ConcatInputStream(in));
	EmapEntry E;

	while ( *ISI && ISI->peek() != std::istream::traits_type::eof() )
	{
		E.deserialise(*ISI);
		SSBOF->put(E);
	}
	ISI.reset();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile< EmapEntry,std::greater<EmapEntry> >::merger_ptr_type Pmerger(SSBOF->getMerger());
	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(out));
	while ( Pmerger->getNext(E) )
		E.serialise(*OSI);

	OSI->flush();
	OSI.reset();

	Pmerger.reset();
	SSBOF.reset();

	libmaus2::aio::FileRemoval::removeFile(tmp);
}

struct MarkedInterval
{
	uint64_t id;
	uint64_t from;
	uint64_t to;

	MarkedInterval() {}
	MarkedInterval(uint64_t const rid, uint64_t const rfrom, uint64_t const rto)
	: id(rid), from(rfrom), to(rto) {}

	MarkedInterval(std::istream & in)
	:
		id(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		from(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		to(libmaus2::util::NumberSerialisation::deserialiseNumber(in))
	{

	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,id);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,from);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,to);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = MarkedInterval(in);
		return in;
	}

	bool operator<(MarkedInterval const & M) const
	{
		if ( id != M.id )
			return id < M.id;
		else if ( from != M.from )
			return from < M.from;
		else
			return to < M.to;
	}
};

struct MarkedIntervalInput
{
	typedef MarkedIntervalInput this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus2::aio::ConcatInputStream ISI;

	MarkedInterval peekslot;
	bool peekslotactive;

	MarkedIntervalInput(std::vector<std::string> const & fn)
	: ISI(fn), peekslotactive(false)
	{

	}

	bool getNext(MarkedInterval & M)
	{
		if ( peekslotactive )
		{
			M = peekslot;
			peekslotactive = false;
			return true;
		}
		else if ( ISI.peek() != std::istream::traits_type::eof() )
		{
			M.deserialise(ISI);
			return true;
		}
		else
		{
			return false;
		}
	}

	bool peekNext(MarkedInterval & M)
	{
		if ( ! peekslotactive )
			peekslotactive = getNext(peekslot);

		if ( peekslotactive )
		{
			M = peekslot;
			return true;
		}
		else
		{
			return false;
		}
	}
};

void mergeMIV(std::string const & fn)
{
	std::string const outfn = fn + ".merge";
	libmaus2::util::TempFileRemovalContainer::addTempFile(outfn);

	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI(new libmaus2::aio::OutputStreamInstance(outfn));
	libmaus2::aio::InputStreamInstance::unique_ptr_type ISI(new libmaus2::aio::InputStreamInstance(fn));

	MarkedInterval E;
	std::vector < MarkedInterval > V;

	while ( *ISI && ISI->peek() != std::istream::traits_type::eof() )
	{
		E.deserialise(*ISI);

		if ( V.size() && V.back().id != E.id )
		{
			uint64_t const id = V.back().id;
			std::vector< libmaus2::math::IntegerInterval<int64_t> > IV;
			for ( uint64_t i = 0; i < V.size(); ++i )
				IV.push_back(libmaus2::math::IntegerInterval<int64_t>(V[i].from,V[i].to));
			IV = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(IV);
			for ( uint64_t i = 0; i < IV.size(); ++i )
				MarkedInterval(id,IV[i].from,IV[i].to).serialise(*OSI);

			V.resize(0);
		}

		V.push_back(E);
	}

	if ( V.size() )
	{
		uint64_t const id = V.back().id;
		std::vector< libmaus2::math::IntegerInterval<int64_t> > IV;
		for ( uint64_t i = 0; i < V.size(); ++i )
			IV.push_back(libmaus2::math::IntegerInterval<int64_t>(V[i].from,V[i].to));
		IV = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(IV);
		for ( uint64_t i = 0; i < IV.size(); ++i )
			MarkedInterval(id,IV[i].from,IV[i].to).serialise(*OSI);
	}

	OSI->flush();
	OSI.reset();
	ISI.reset();

	libmaus2::aio::OutputStreamFactoryContainer::rename(outfn,fn);
}

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

void printResult(
	libmaus2::util::ArgParser const & arg,
	std::string const & bamfn,
	std::vector<std::string> const & TCONStmp,
	std::vector<std::string> const & TEMAPtmp,
	std::vector<std::string> const & TMIVtmp
)
{
	std::string const tmpfilebase = getTmpFileBase(arg);

	libmaus2::lcs::AlignmentStatistics GAS;
	{
		ConsensusEntry CE;
		libmaus2::aio::ConcatInputStream ISI(TCONStmp);
		while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
		{
			CE.deserialise(ISI);
			GAS += CE.AS;
			std::cout << CE.id << "\t" << CE.AS << "\t" << GAS << "\t" << CE.IN << "\t" << CE.refid << ":" << CE.refpos << "," << CE.refpos+CE.reflen << std::endl;
		}
	}

	uint64_t grlen = 0;
	uint64_t gclen = 0;

	MarkedIntervalInput::unique_ptr_type PMII(new MarkedIntervalInput(TMIVtmp));

	libmaus2::bambam::BamDecoder bamdec(bamfn);
	for ( uint64_t i = 0; bamdec.readAlignment(); ++i )
	{
		libmaus2::bambam::BamAlignment const & algn = bamdec.getAlignment();
		uint64_t const rlen = algn.getReferenceLength();
		uint64_t clen = 0;

		std::vector < libmaus2::math::IntegerInterval<int64_t> > V;
		MarkedInterval M;
		while ( PMII->peekNext(M) && M.id == i )
		{
			V.push_back(libmaus2::math::IntegerInterval<int64_t>(M.from,M.to));
			PMII->getNext(M);
		}

		for ( uint64_t j = 0; j < V.size(); ++j )
			clen += V[j].diameter();

		std::cout << "[C]\t" << i << "\t" << rlen << "\t" << clen << "\t" << static_cast<double>(clen)/static_cast<double>(rlen) << std::endl;

		grlen += rlen;
		gclen += clen;
	}

	std::cout << "[G]\t" << grlen << "\t" << gclen << "\t" << static_cast<double>(gclen)/static_cast<double>(grlen) << "\t" << GAS << std::endl;

	std::string const TEMAPrev = tmpfilebase + ".emap_rev";
	libmaus2::util::TempFileRemovalContainer::addTempFile(TEMAPrev);
	std::string const TEMAPrevtmp = tmpfilebase + ".emap_revtmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(TEMAPrevtmp);
	reverseEmap(TEMAPtmp,TEMAPrev, TEMAPrevtmp);

	{
		uint64_t s = 0;

		{
			libmaus2::aio::ConcatInputStream ISI(TEMAPtmp);
			EmapEntry E;
			while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
			{
				E.deserialise(ISI);
				s += E.i;
			}
		}

		{
			uint64_t a = 0;

			libmaus2::aio::InputStreamInstance ISI(TEMAPrev);
			EmapEntry E;
			while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
			{
				E.deserialise(ISI);
				assert ( E.i != 0 );
				a += E.i;
				std::cout << "[EM]\t" << printDoubleFixed(E.v,10) << "\t" << a << "\t" << static_cast<double>(a)/s << std::endl;
			}
		}

		{
			uint64_t a = 0;

			libmaus2::aio::ConcatInputStream ISI(TEMAPtmp);
			EmapEntry E;
			while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
			{
				E.deserialise(ISI);
				assert ( E.i != 0 );
				a += E.i;
				std::cout << "[EP]\t" << printDoubleFixed(E.v,10) << "\t" << a << "\t" << static_cast<double>(a)/s << std::endl;
			}
		}
	}

	libmaus2::aio::FileRemoval::removeFile(TEMAPrev);
	libmaus2::aio::FileRemoval::removeFile(TEMAPrevtmp);
}

static uint64_t const fastaindexmod = 64;

static std::string getFastAIndexFileName(std::string const & consfn)
{
	std::string const consfnindex = consfn + ".checkconsensus.index";
	return consfnindex;
}

static std::string getNumericalIndexName(std::string const & bamfn)
{
	return libmaus2::bambam::BamNumericalIndexBase::getIndexName(bamfn);
}

int generateindex(
	libmaus2::util::ArgParser const & arg,
	std::string const & /* textfn */,
	std::string const & bamfn,
	std::string const & consfn
)
{
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

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

	libmaus2::bambam::BamNumericalIndexGenerator::indexFileCheck(bamfn,64,numthreads,true);

	return EXIT_SUCCESS;
}

struct CheckResultSet
{
	std::string cons;
	std::string emap;
	std::string miv;

	CheckResultSet() {}
	CheckResultSet(
		std::string const & rcons,
		std::string const & remap,
		std::string const & rmiv
	) : cons(rcons), emap(remap), miv(rmiv) {}
	CheckResultSet(std::istream & in)
	:
		cons(libmaus2::util::StringSerialisation::deserialiseString(in)),
		emap(libmaus2::util::StringSerialisation::deserialiseString(in)),
		miv(libmaus2::util::StringSerialisation::deserialiseString(in))
	{

	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::StringSerialisation::serialiseString(out,cons);
		libmaus2::util::StringSerialisation::serialiseString(out,emap);
		libmaus2::util::StringSerialisation::serialiseString(out,miv);
		return out;
	}

	void serialise(std::string const & fn) const
	{
		libmaus2::aio::OutputStreamInstance OSI(fn);
		serialise(OSI);
		OSI.flush();
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = CheckResultSet(in);
		return in;
	}

	void remove()
	{
		libmaus2::aio::FileRemoval::removeFile(cons);
		libmaus2::aio::FileRemoval::removeFile(emap);
		libmaus2::aio::FileRemoval::removeFile(miv);
	}
};

CheckResultSet checkconsensus(
	libmaus2::util::ArgParser const & arg,
	std::string const & textfn,
	std::string const & bamfn,
	std::string const & consfn,
	std::string const & outputprefix,
	uint64_t const ilow = 0,
	uint64_t const ihigh = std::numeric_limits<uint64_t>::max()
)
{
	libmaus2::util::ToUpperTable const toup;

	std::cerr << "[V] loading ref seq...";
	std::vector<std::string> refnames;
	std::vector<std::string> const Vtext = loadTextVector(textfn,refnames);
	std::cerr << "done" << std::endl;

	std::map<std::string,uint64_t> refmap;
	std::map<uint64_t,uint64_t> refremap;
	for ( uint64_t i = 0; i < refnames.size(); ++i )
		refmap[refnames[i]] = i;

	{
		libmaus2::bambam::BamDecoder bamdec(bamfn);
		libmaus2::bambam::BamHeader const & header = bamdec.getHeader();
		//assert ( header.getNumRef() <= Vtext.size() );
		for ( uint64_t i = 0; i < header.getNumRef(); ++i  )
		{
			std::string const refname = header.getRefIDName(i);
			assert ( refmap.find(refname) != refmap.end() );
			refremap[i] = refmap.find(refname)->second;
		}
	}

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

	std::string const tmpfilebase = getTmpFileBase(arg);

	libmaus2::aio::TempFileArray TEMAP("emap",numthreads);
	libmaus2::aio::TempFileArray TMIV("miv",numthreads);
	libmaus2::aio::TempFileArray TCONS("cons",numthreads);

	std::string const consfnindex = getFastAIndexFileName(consfn);
	std::string const bamindexfn = getNumericalIndexName(bamfn);
	libmaus2::bambam::BamNumericalIndexDecoder bamindex(bamindexfn);

	uint64_t const rlow = std::max(ilow,getReadsLow(consfnindex));
	uint64_t const rhigh = std::min(ihigh,getReadsHigh(consfnindex));
	uint64_t const numreads = rhigh-rlow;
	//std::cerr << "numreads=" << numreads << std::endl;

	uint64_t const tnumpacks = numthreads * 4;
	uint64_t const readsperpack = (numreads + tnumpacks - 1)/tnumpacks;
	uint64_t const numpacks = readsperpack ? ((numreads + readsperpack - 1)/readsperpack) : 0;

	//std::cerr << "numpacks=" << numpacks << std::endl;

	libmaus2::lcs::AlignmentStatistics TGAS;
	libmaus2::parallel::PosixSpinLock TGASlock;

	int volatile failed = 0;
	libmaus2::parallel::PosixSpinLock failedlock;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t pack = 0; pack < numpacks; ++pack )
	{
		try
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			// std::ostream & out = *(Aout[tid]);

			uint64_t const packlow = rlow + pack * readsperpack;
			uint64_t const packhigh = std::min(packlow+readsperpack,rhigh);
			assert ( packhigh > packlow );

			ReadAccessor RA(consfn,consfnindex,fastaindexmod);

			uint64_t const idlow = RA.getReadId(packlow);
			libmaus2::lz::BgzfInflateFile::unique_ptr_type tptr(bamindex.getStreamAt(bamfn,idlow));

			uint64_t const b = packlow / fastaindexmod;
			libmaus2::aio::InputStreamInstance IISI(consfnindex);
			IISI.seekg(b * sizeof(libmaus2::fastx::FastInterval));
			libmaus2::fastx::FastInterval FI(IISI);
			libmaus2::aio::InputStreamInstance ISI(consfn);
			ISI.seekg(FI.fileoffset);
			libmaus2::fastx::StreamFastAReaderWrapper SFAR(ISI);
			uint64_t o = packlow - b * fastaindexmod;
			for ( uint64_t j = 0; j < o; ++j )
				SFAR.skipPattern();
			libmaus2::fastx::FastAReader::pattern_type P;

			#if 0
			bool const ok = SFAR.getNextPatternUnlocked(P);
			assert ( ok );
			patternToUpper(P);
			#endif

			libmaus2::bambam::BamAlignment algn;
			libmaus2::bambam::BamAlignmentDecoder::readAlignmentGz(*tptr,algn);

			for ( uint64_t i = packlow; i < packhigh; ++i )
			{
				bool const ok = SFAR.getNextPatternUnlocked(P);
				assert ( ok );
				patternToUpper(P);

				// std::cerr << "new name " << P.sid << std::endl;

				uint64_t const id = ReadAccessor::getReadId(P.sid);
				assert ( getBamId(algn.getName()) <= id );

				while ( getBamId(algn.getName()) < id )
				{
					bool const ok = libmaus2::bambam::BamAlignmentDecoder::readAlignmentGz(*tptr,algn);
					assert ( ok );
				}


				// std::cerr << P.sid << " " << getBamId(algn.getName()) << std::endl;

				assert ( getBamId(algn.getName()) == id );

				// libmaus2::bambam::BamAlignment const & algn = *(Vbam[id]);
				uint64_t const refid = refremap.find(algn.getRefID())->second;
				std::string const & text = Vtext [ refid ];

				// part of text covered by read
				bool const bamrev = algn.isReverse();
				uint64_t const refstretchpos = algn.getPos()-algn.getFrontDel();
				//std::cerr << "id=" << id << " pos=" << refstretchpos << " refid=" << refid << " len " << text.size() << std::endl;
				std::string const psreftext = toup(text.substr(refstretchpos,algn.getReferenceLength()));
				//std::cerr << "got it" << std::endl;
				// reference text matched by A read
				std::string const sreftext = bamrev ? libmaus2::fastx::reverseComplementUnmapped(psreftext) : psreftext;

				// std::cerr << id << std::endl;

				uint8_t const * uref = reinterpret_cast<uint8_t const *>(sreftext.c_str());
				uint8_t const * urefe = uref + sreftext.size();

				uint8_t const * uSCO = reinterpret_cast<uint8_t const *>(P.spattern.c_str());
				uint8_t const * uSCOe = uSCO + P.spattern.size();

				libmaus2::lcs::SuffixArrayLCS::LCSResult const lcsres =
					libmaus2::lcs::SuffixArrayLCS::lcsmin(std::string(uref,urefe),std::string(uSCO,uSCOe));
				libmaus2::lcs::NNP nnp;
				libmaus2::lcs::NNPTraceContainer nnptrace;
				libmaus2::lcs::NNPAlignResult nnpres = nnp.align(uref,urefe,lcsres.maxpos_a,uSCO,uSCOe,lcsres.maxpos_b,nnptrace);
				//nnptrace.printTraceLines(err,uref + nnpres.abpos,uSCO + nnpres.bbpos);

				if ( nnpres.bepos != P.spattern.size() && nnpres.aepos != sreftext.size() )
				{
					assert ( nnpres.bepos < P.spattern.size() );

					libmaus2::lcs::NPL npl;
					npl.np(
						uref+nnpres.abpos,
						urefe,
						uSCO+nnpres.bbpos,
						uSCOe
					);

					std::pair<uint64_t,uint64_t> const P = npl.getTraceContainer().getStringLengthUsed();

					nnpres.aepos = nnpres.abpos + P.first;
					nnpres.bepos = nnpres.bbpos + P.second;
				}
				if ( nnpres.bbpos != 0 && nnpres.abpos != 0 )
				{
					assert ( nnpres.bbpos > 0 );

					libmaus2::lcs::NPL npl;
					npl.np(
						::std::reverse_iterator<uint8_t const *>(uref + nnpres.aepos),
						::std::reverse_iterator<uint8_t const *>(uref),
						::std::reverse_iterator<uint8_t const *>(uSCO + nnpres.bepos),
						::std::reverse_iterator<uint8_t const *>(uSCO)
					);

					std::pair<uint64_t,uint64_t> const P = npl.getTraceContainer().getStringLengthUsed();

					nnpres.abpos = nnpres.aepos - P.first;
					nnpres.bbpos = nnpres.bepos - P.second;
				}

				libmaus2::lcs::NP np;
				np.np(
					uref + nnpres.abpos,
					uref + nnpres.aepos,
					uSCO + nnpres.bbpos,
					uSCO + nnpres.bepos
				);

				libmaus2::lcs::AlignmentTraceContainer const & ATC = np.getTraceContainer();

				int64_t const orefpos = refstretchpos + ( bamrev ? (sreftext.size() - nnpres.aepos) : nnpres.abpos );
				int64_t const oreflen = nnpres.aepos - nnpres.abpos;

				#if 0
				libmaus2::lcs::AlignmentTraceContainer ATC;
				nnptrace.computeTrace(ATC);
				#endif

				std::ostringstream algnstr;
				libmaus2::lcs::AlignmentPrint::printAlignmentLines(
					algnstr,
					uref+nnpres.abpos,
					nnpres.aepos-nnpres.abpos,
					uSCO+nnpres.bbpos,
					nnpres.bepos-nnpres.bbpos,
					80,
					ATC.ta,
					ATC.te
				);

				libmaus2::lcs::AlignmentStatistics const AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(ATC.ta,ATC.te);

				EmapEntry(AS.getErrorRate(),1).serialise(TEMAP[tid]);

				libmaus2::math::IntegerInterval<int64_t> IN(nnpres.abpos,static_cast<int64_t>(nnpres.aepos)-1);

				ConsensusEntry(id,AS,IN,refid,orefpos,oreflen).serialise(TCONS[tid]);

				{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				TGAS += AS;
				std::cerr << "[A]\t" << id << "\t" << AS << "\t" << "R(" << refid << ":" << orefpos << "," << orefpos+oreflen << ")" << "\t" << IN << "\t" << TGAS << std::endl;
				}

				MarkedInterval(id,IN.from,IN.to).serialise(TMIV[tid]);
			}
		}
		catch(std::exception const & ex)
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << ex.what();
			failedlock.lock();
			failed = 1;
			failedlock.unlock();
		}
	}

	if ( failed )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] failed in OMP loop" << std::endl;
		lme.finish();
		throw lme;
	}

	TCONS.flush();
	std::string const TCONStmp = tmpfilebase + ".cons_tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(TCONStmp);
	std::less<ConsensusEntry> consorder;
	TCONS.merge< ConsensusEntry, std::less<ConsensusEntry> >(TCONStmp,consorder);

	TEMAP.flush();
	std::string const TEMAPtmp = tmpfilebase + ".emap_tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(TEMAPtmp);
	std::less<EmapEntry> emaporder;
	TEMAP.merge< EmapEntry,std::less<EmapEntry> >(TEMAPtmp,emaporder);
	mergeEmap(TEMAPtmp);

	TMIV.flush();
	std::string const TMIVtmp = tmpfilebase + ".miv_tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(TMIVtmp);
	std::less<MarkedInterval> MIVorder;
	TMIV.merge < MarkedInterval, std::less<MarkedInterval> >(TMIVtmp,MIVorder);
	mergeMIV(TMIVtmp);

	std::string const consout = outputprefix + ".cons";
	std::string const emapout = outputprefix + ".emap";
	std::string const mivout = outputprefix + ".miv";

	libmaus2::aio::OutputStreamFactoryContainer::rename(TCONStmp,consout);
	libmaus2::aio::OutputStreamFactoryContainer::rename(TEMAPtmp,emapout);
	libmaus2::aio::OutputStreamFactoryContainer::rename(TMIVtmp,mivout);

	#if 0
	printResult(arg,bamfn,TCONStmp,TEMAPtmp,TMIVtmp);
	libmaus2::aio::FileRemoval::removeFile(TCONStmp);
	libmaus2::aio::FileRemoval::removeFile(TEMAPtmp);
	libmaus2::aio::FileRemoval::removeFile(TMIVtmp);
	#endif

	return CheckResultSet(consout,emapout,mivout);
}

uint64_t parseInt(std::string const & s)
{
	std::istringstream istr(s);

	uint64_t u;
	istr >> u;

	if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] failed to parse " << s << std::endl;
		lme.finish();
		throw lme;
	}

	return u;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 4 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] <index|check|batchlist|batchprocess|batchmerge|cleanup> ref.fasta reads.bam reads_cons.fasta\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			if ( arg[0] == "index" )
			{
				return generateindex(arg,arg[1],arg[2],arg[3]);
			}
			else if ( arg[0] == "check" )
			{
				std::string const tmpfilebase = getTmpFileBase(arg);
				std::ostringstream outprefixstr;
				outprefixstr << tmpfilebase << "_output";
				std::string const outprefix = outprefixstr.str();
				CheckResultSet CRS = checkconsensus(arg,arg[1],arg[2],arg[3],outprefix);

				printResult(arg,arg[2],std::vector<std::string>(1,CRS.cons),std::vector<std::string>(1,CRS.emap),std::vector<std::string>(1,CRS.miv));

				CRS.remove();

				return EXIT_SUCCESS;
			}
			else if ( arg[0] == "batchlist" )
			{
				std::string const consfn = arg[3];
				std::string const consfnindex = getFastAIndexFileName(consfn);
				uint64_t const rlow = getReadsLow(consfnindex);
				uint64_t const rhigh = getReadsHigh(consfnindex);
				uint64_t const range = rhigh - rlow;

				uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

				uint64_t const tnumpacks = arg.uniqueArgPresent("p") ? arg.getUnsignedNumericArg<uint64_t>("p") : getDefaultNumPacks();
				uint64_t const packsize = (range + tnumpacks - 1) / tnumpacks;
				uint64_t const numpacks = packsize ? ((range + packsize - 1)/packsize) : 0;

				std::string const tmp = getTmpFileBase(arg);

				for ( uint64_t i = 0; i < numpacks; ++i )
				{
					uint64_t const low = rlow + i*packsize;
					uint64_t const high = std::min(low+packsize,rhigh);

					std::ostringstream tmpstr;
					tmpstr << tmp << "_" << std::setw(6) << std::setfill('0') << i;
					std::string const subtmp = tmpstr.str();

					std::cout << arg.progname << " -T" << subtmp << " -t" << numthreads << " -I" << low << "," << high << " batchprocess "
						<< arg[1] << " " << arg[2] << " " << arg[3] << " " << consfn << "_check_" << std::setw(6) << std::setfill('0') << i << std::setw(0) << std::endl;
				}

				std::cout << arg.progname << " -T" << tmp << " -t" << numthreads << " -p" << tnumpacks << " batchmerge "
					<< arg[1] << " " << arg[2] << " " << arg[3] << std::endl;
				std::cout << arg.progname << " -T" << tmp << " -t" << numthreads << " -p" << tnumpacks << " cleanup "
					<< arg[1] << " " << arg[2] << " " << arg[3] << std::endl;
			}
			else if ( arg[0] == "batchprocess" )
			{
				if ( ! arg.uniqueArgPresent("I") )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] missing argument I" << std::endl;
					lme.finish();
					throw lme;
				}

				std::string Is = arg["I"];
				if ( Is.find(',') == std::string::npos )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "[E] argument I malformed" << std::endl;
					lme.finish();
					throw lme;
				}

				uint64_t const I0 = parseInt(Is.substr(0,Is.find(',')));
				uint64_t const I1 = parseInt(Is.substr(Is.find(',')+1));

				CheckResultSet const CRS = checkconsensus(arg,arg[1],arg[2],arg[3],arg[4],I0,I1);
				CRS.serialise(arg[4]);
			}
			else if ( arg[0] == "batchmerge" )
			{
				std::string const consfn = arg[3];
				std::string const consfnindex = getFastAIndexFileName(consfn);
				uint64_t const rlow = getReadsLow(consfnindex);
				uint64_t const rhigh = getReadsHigh(consfnindex);
				uint64_t const range = rhigh - rlow;

				uint64_t const tnumpacks = arg.uniqueArgPresent("p") ? arg.getUnsignedNumericArg<uint64_t>("p") : 1;
				uint64_t const packsize = (range + tnumpacks - 1) / tnumpacks;
				uint64_t const numpacks = packsize ? ((range + packsize - 1)/packsize) : 0;

				std::vector < CheckResultSet > V;
				std::vector < std::string > Vcons;
				std::vector < std::string > Vemap;
				std::vector < std::string > Vmiv;

				for ( uint64_t i = 0; i < numpacks; ++i )
				{
					std::ostringstream fnstr;
					fnstr << consfn << "_check_" << std::setw(6) << std::setfill('0') << i;
					libmaus2::aio::InputStreamInstance ISI(fnstr.str());
					CheckResultSet const CRS(ISI);
					V.push_back(CRS);

					Vcons.push_back(CRS.cons);
					Vemap.push_back(CRS.emap);
					Vmiv.push_back(CRS.miv);
				}

				std::string const emapmergefn = consfn + "_check.emap";
				libmaus2::sorting::SerialisingSortingBufferedOutputFile<EmapEntry>::reduce(Vemap,emapmergefn);

				printResult(arg,arg[2],Vcons,std::vector<std::string>(1,emapmergefn),Vmiv);

				libmaus2::aio::FileRemoval::removeFile(emapmergefn);
			}
			else if ( arg[0] == "cleanup" )
			{
				std::string const consfn = arg[3];
				std::string const consfnindex = getFastAIndexFileName(consfn);
				uint64_t const rlow = getReadsLow(consfnindex);
				uint64_t const rhigh = getReadsHigh(consfnindex);
				uint64_t const range = rhigh - rlow;

				uint64_t const tnumpacks = arg.uniqueArgPresent("p") ? arg.getUnsignedNumericArg<uint64_t>("p") : 1;
				uint64_t const packsize = (range + tnumpacks - 1) / tnumpacks;
				uint64_t const numpacks = packsize ? ((range + packsize - 1)/packsize) : 0;

				std::vector < CheckResultSet > V;
				std::vector < std::string > remlist;
				for ( uint64_t i = 0; i < numpacks; ++i )
				{
					std::ostringstream fnstr;
					fnstr << consfn << "_check_" << std::setw(6) << std::setfill('0') << i;
					libmaus2::aio::InputStreamInstance ISI(fnstr.str());
					CheckResultSet const CRS(ISI);
					V.push_back(CRS);

					remlist.push_back(CRS.cons);
					remlist.push_back(CRS.emap);
					remlist.push_back(CRS.miv);
					remlist.push_back(fnstr.str());
				}

				for ( uint64_t i = 0; i < remlist.size(); ++i )
					libmaus2::aio::FileRemoval::removeFile(remlist[i]);
			}
			else
			{
				std::cerr << "[E] unknown subcommand " << arg[0] << std::endl;
				return EXIT_FAILURE;
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
