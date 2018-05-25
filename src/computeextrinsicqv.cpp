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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NNPCorL.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <config.h>

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

int64_t getId(libmaus2::fastx::FastAReader::pattern_type const & pb)
{
	std::string sid = pb.getShortStringId();

	if ( sid.find('/') != std::string::npos )
	{
		sid = sid.substr(0,sid.find('/'));

		std::istringstream istr(sid);
		int64_t id;
		istr >> id;

		if ( istr && istr.peek() == std::istream::traits_type::eof() )
			return id-1;
		else
			return -1;
	}
	else
	{
		return -1;
	}
}

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}


int64_t getDefaultTSpace()
{
	return 100;
}

static std::string helpMessage(libmaus2::util::ArgParser const & /* arg */)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
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

struct IndexEntry
{
	uint64_t offset;
	uint64_t count;

	IndexEntry()
	{}

	IndexEntry(uint64_t const roffset, uint64_t const rcount) : offset(roffset), count(rcount) {}
	IndexEntry(std::istream & in) { deserialise(in); }

	std::istream & deserialise(std::istream & in)
	{
		offset = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		count = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		return in;
	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,offset);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,count);
		return out;
	}
};

static std::vector<IndexEntry> indexConsFile(std::string const & fn, uint64_t const nr)
{
	std::vector<IndexEntry> V(nr);
	libmaus2::fastx::FastAReader FA(fn);
	libmaus2::fastx::FastAReader::pattern_type pattern;

	int64_t prevread = std::numeric_limits<int64_t>::min();
	uint64_t prevcount = 0;
	uint64_t prevo = 0;
	uint64_t next = 0;

	while ( FA.foundnextmarker )
	{
		uint64_t const o = FA.getC()-1;

		bool const ok = FA.getNextPatternUnlocked(pattern);
		assert ( ok );

		int64_t const id = getId(pattern);
		assert ( id >= 0 );

		if ( id != prevread )
		{
			if ( prevread >= 0 )
			{
				while ( static_cast<int64_t>(next) < prevread )
				{
					V[next] = IndexEntry(0,0);
					++next;
				}
				// std::cerr << "o=" << prevo << " id=" << prevread << " count=" << prevcount << std::endl;
				V[next] = IndexEntry(prevo,prevcount);

				++next;
			}

			prevo = o;
			prevread = id;
			prevcount = 0;
		}

		++prevcount;
	}

	if ( prevread >= 0 )
	{
		while ( static_cast<int64_t>(next) < prevread )
		{
			V[next] = IndexEntry(0,0);
			++next;
		}
		// std::cerr << "o=" << prevo << " id=" << prevread << " count=" << prevcount << std::endl;
		V[next] = IndexEntry(prevo,prevcount);
		++next;
	}

	while ( next < nr )
	{
		V[next] = IndexEntry(0,0);
		++next;
	}

	return V;
}

struct StringId
{
	uint64_t id;
	uint64_t len;
	libmaus2::autoarray::AutoArray<uint8_t> A;

	StringId()
	{

	}

	template<typename iterator>
	StringId(uint64_t const rid, uint64_t const rlen, iterator rA)
	: id(rid), len(rlen), A(len,false)
	{
		std::copy(rA,rA+rlen,A.begin());
	}

	StringId(std::istream & in)
	{
		deserialise(in);
	}

	std::istream & deserialise(std::istream & in)
	{
		id = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		len = libmaus2::util::NumberSerialisation::deserialiseNumber(in);

		A.resize(len);
		in.read(reinterpret_cast<char *>(A.begin()),len);
		assert ( in.gcount() == static_cast<int64_t>(len) );

		return in;
	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,id);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,len);
		out.write(reinterpret_cast<char const *>(A.begin()),len);
		return out;
	}

	bool operator<(StringId const & O) const
	{
		return id < O.id;
	}
};

int computequality(libmaus2::util::ArgParser const & arg)
{
	std::string const db0name = arg[1];

	std::cerr << "[V] loading data for " << db0name << " to memory...";
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type Pdb0data(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(db0name)
	);
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet const * db0data = Pdb0data.get();
	std::cerr << "done." << std::endl;

	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB0(
		new libmaus2::dazzler::db::DatabaseFile(db0data->getDBURL())
	);
	libmaus2::dazzler::db::DatabaseFile & DB = *PDB0;
	DB.computeTrimVector();
	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	std::vector<IndexEntry> const Vindex = indexConsFile(arg[0],DB.size());

	if ( arg.uniqueArgPresent("indexonly") )
		return EXIT_SUCCESS;

	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	std::string const tmpfilebase = getTmpFileBase(arg);

	std::vector < std::string > Vtmp(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::aio::OutputStreamInstance::unique_ptr_type > Aout(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		std::ostringstream fnostr;
		fnostr << tmpfilebase << "_" << i << "_datatmp";
		std::string const fn = fnostr.str();
		Vtmp[i] = fn;
		libmaus2::util::TempFileRemovalContainer::addTempFile(fn);

		libmaus2::aio::OutputStreamInstance::unique_ptr_type tptr(
			new libmaus2::aio::OutputStreamInstance(fn)
		);

		Aout[i] = UNIQUE_PTR_MOVE(tptr);
	}

	int64_t const tspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : getDefaultTSpace();

	// int64_t first = -1;

	int64_t Icnt = 0;
	int64_t Idiv = 1;
	if ( arg.uniqueArgPresent("J") )
	{

		std::string const Js = arg["J"];
		std::istringstream istr(Js);
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

		istr >> Idiv;

		if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Js << std::endl;
			lme.finish();
			throw lme;
		}

		if ( !Idiv )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] denominator of J argument cannot be zero" << std::endl;
			lme.finish();
			throw lme;
		}
	}

	uint64_t const n = DB.size();
	uint64_t const readsperpack = (n + Idiv - 1) / Idiv;

	uint64_t const alow = std::min(Icnt * readsperpack,n);
	uint64_t const ahigh = std::min(alow + readsperpack,n);

	libmaus2::autoarray::AutoArray<uint64_t> annosize(ahigh-alow);
	std::fill(annosize.begin(),annosize.end(),0ull);

	std::cerr << "[V] processing [" << alow << "," << ahigh << ")" << std::endl;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t aid = alow; aid < ahigh; ++aid )
	{
		uint64_t const l = RL[aid];
		// number of trace point intervals
		uint64_t const nt = (l + tspace - 1)/tspace;
		annosize[aid - alow] = nt;
		std::vector<uint8_t> V(nt,std::numeric_limits<uint8_t>::max());

		libmaus2::fastx::FastAReader::pattern_type pb;
		libmaus2::lcs::NNPCorL np;

		libmaus2::aio::InputStreamInstance ISI(arg[0]);
		ISI.seekg(Vindex[aid].offset);
		libmaus2::fastx::StreamFastAReaderWrapper SFC(ISI);

		std::string const afull = DB[aid];

		for ( uint64_t z = 0; z < Vindex[aid].count; ++z )
		{
			bool const ok = SFC.getNextPatternUnlocked(pb);
			assert ( ok );
			assert ( getId(pb) == static_cast<int64_t>(aid) );
			// FC.getNext(pb);

			assert ( pb.sid.find("A=[") != std::string::npos );
			std::string sid = pb.sid;
			sid = sid.substr(sid.find("A=[") + strlen("A=["));
			assert ( sid.find("]") != std::string::npos );
			sid = sid.substr(0,sid.find("]"));

			std::istringstream istr(sid);
			int64_t from, to;

			istr >> from;
			assert ( istr.peek() == ',' );
			istr.get();
			istr >> to;
			to += 1;
			to = std::min(to,static_cast<int64_t>(RL[aid]));
			assert ( istr.peek() == std::istream::traits_type::eof() );

			std::string const asub = afull.substr(from,to-from);
			std::string const bsub = pb.spattern;

			np.np(asub.begin(),asub.end(),bsub.begin(),bsub.end());

			#if 0
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				asub.begin(),
				asub.size(),
				bsub.begin(),
				bsub.size(),
				80,
				np.ta,
				np.te
			);
			#endif

			if ( from % tspace != 0 )
			{
				std::pair<uint64_t,uint64_t> const P = libmaus2::lcs::AlignmentTraceContainer::advanceA(np.ta,np.te,tspace - (from % tspace));

				from += P.first;
				np.ta += P.second;
			}

			assert ( from == to || (from % tspace == 0) );

			while ( to-from >= tspace )
			{
				assert ( from % tspace == 0 );

				std::pair<uint64_t,uint64_t> const P = libmaus2::lcs::AlignmentTraceContainer::advanceA(np.ta,np.te,tspace);

				bool const ok = static_cast<int64_t>(P.first) == tspace;
				if ( ! ok )
				{
					std::cerr << "Failure for " << aid << " tspace=" << tspace << " P.first=" << P.first << " from=" << from << " to=" << to << " RL=" << RL[aid] << std::endl;
				}
				assert ( ok );

				libmaus2::lcs::AlignmentStatistics AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(np.ta,np.ta+P.second);

				uint64_t const e = std::floor(AS.getErrorRate() * std::numeric_limits<uint8_t>::max() + 0.5);
				assert ( e <= std::numeric_limits<uint8_t>::max() );

				V [ from / tspace ] = std::min(V[from/tspace],static_cast<uint8_t>(e));

				from += tspace;
				np.ta += P.second;
			}

			if ( (from % tspace == 0) && (to > from) && (to == static_cast<int64_t>(l)) )
			{
				uint64_t const d = to-from;
				assert ( static_cast<int64_t>(d) < tspace );

				std::pair<uint64_t,uint64_t> const P = libmaus2::lcs::AlignmentTraceContainer::advanceA(np.ta,np.te,d);
				assert ( static_cast<int64_t>(P.first) == static_cast<int64_t>(d) );

				libmaus2::lcs::AlignmentStatistics AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(np.ta,np.ta+P.second);

				uint64_t const e = std::floor(AS.getErrorRate() * std::numeric_limits<uint8_t>::max() + 0.5);
				assert ( e <= std::numeric_limits<uint8_t>::max() );

				V [ from / tspace ] = std::min(V[from/tspace],static_cast<uint8_t>(e));

				from += d;
				np.ta += P.second;
			}
		}

		#if 0
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			// std::cerr << aid << " " << i << " " << (double)V[i]/255.0 << std::endl;
			dataOSI.put(V[i]);
		}
		#endif

		StringId SI(aid,V.size(),V.begin());

		#if defined(_OPENMP)
		uint64_t const tid = omp_get_thread_num();
		#else
		uint64_t const tid = 0;
		#endif

		SI.serialise(*(Aout[tid]));

		if ( (aid-alow) % 1024 == 0 )
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << "[V] " << (aid-alow) << std::endl;
		}
	}

	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		Aout[i]->flush();
		Aout[i].reset();
	}

	std::cerr << "[V] " << (ahigh-alow) << std::endl;

	PDB0.reset();
	libmaus2::dazzler::db::DatabaseFile ondiskDB(db0name);

	std::string const annofn = ondiskDB.getBlockTrackAnnoFileName("exqual",(Idiv > 1) ? (Icnt+1) : 0);
	std::string const datafn = ondiskDB.getBlockTrackDataFileName("exqual",(Idiv > 1) ? (Icnt+1) : 0);

	std::ostringstream fnostr;
	fnostr << tmpfilebase << "_mergetmp";
	std::string const mergetmpfn = fnostr.str();
	libmaus2::util::TempFileRemovalContainer::addTempFile(mergetmpfn);
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<StringId>::reduce(Vtmp,mergetmpfn);

	for ( uint64_t i = 0; i < numthreads; ++i )
		libmaus2::aio::FileRemoval::removeFile(Vtmp[i]);

	libmaus2::aio::InputStreamInstance dataISI(mergetmpfn);
	libmaus2::aio::OutputStreamInstance dataOSI(datafn);

	while ( dataISI && (dataISI.peek() != std::istream::traits_type::eof()) )
	{
		StringId SI(dataISI);
		for ( uint64_t i = 0; i < SI.len; ++i )
			dataOSI.put(SI.A[i]);
	}

	uint64_t const p = dataOSI.tellp();
	assert ( p == std::accumulate(annosize.begin(),annosize.end(),0ull) );

	dataOSI.flush();

	libmaus2::aio::OutputStreamInstance annoOSI(annofn);

	// write inqual anno file
	uint64_t annooff = 0;
	libmaus2::dazzler::db::OutputBase::putLittleEndianInteger4(annoOSI,annosize.size() /* tracklen */,annooff);
	libmaus2::dazzler::db::OutputBase::putLittleEndianInteger4(annoOSI,8 /* size of pointer */,annooff);
	uint64_t s = 0;
	for ( uint64_t i = 0; i < annosize.size(); ++i )
	{
		libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(annoOSI,s,annooff);
		s += annosize[i];
	}
	libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(annoOSI,s,annooff);
	annoOSI.flush();

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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] reads_cons.fasta reads.db\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return computequality(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
