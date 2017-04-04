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
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <config.h>

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


int computequality(libmaus2::util::ArgParser const & arg)
{
	libmaus2::fastx::FastaPeeker FA(arg[0]);
	libmaus2::fastx::FastaPeeker FC(arg[1]);
	std::string const dbname = arg[2];
	libmaus2::dazzler::db::DatabaseFile const DB(dbname);
	uint64_t const part = 3 < arg.size() ? arg.getParsedRestArg<uint64_t>(3) : 0;

	int64_t const tspace = arg.uniqueArgPresent("tspace") ? arg.getUnsignedNumericArg<uint64_t>("tspace") : getDefaultTSpace();

	std::string const annofn = DB.getBlockTrackAnnoFileName("exqual",part);
	std::string const datafn = DB.getBlockTrackDataFileName("exqual",part);

	libmaus2::aio::OutputStreamInstance dataOSI(datafn);

	libmaus2::fastx::FastAReader::pattern_type pa;
	libmaus2::fastx::FastAReader::pattern_type pb;

	// int64_t const tspace = 100;

	libmaus2::lcs::NP np;

	int64_t first = -1;
	std::vector<uint64_t> annosize;

	while ( FA.getNext(pa) )
	{
		std::string sid = pa.getShortStringId();
		if ( sid.find('/') != std::string::npos )
		{
			sid = sid.substr(sid.find('/')+1);

			if ( sid.find('/') != std::string::npos )
				sid = sid.substr(0,sid.find('/'));
			else
				continue;
		}
		else
		{
			continue;
		}

		std::istringstream istr(sid);
		int64_t id;
		istr >> id;

		if ( first < 0 )
			first = id;

		uint64_t const aid = id - first;

		while ( !(aid < annosize.size()) )
			annosize.push_back(0);

		assert ( aid < annosize.size() );

		uint64_t const l = pa.spattern.size();
		uint64_t const nt = (l + tspace - 1)/tspace;
		annosize[aid] = nt;
		std::vector<uint8_t> V(nt,std::numeric_limits<uint8_t>::max());

		while ( FC.peekNext(pb) && getId(pb) == id )
		{
			FC.getNext(pb);

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
			assert ( istr.peek() == std::istream::traits_type::eof() );

			std::string const asub = pa.spattern.substr(from,to-from+1);
			std::string const bsub = pb.spattern;

			np.np(asub.begin(),asub.end(),bsub.begin(),bsub.end());

			/*
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
			*/

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
				assert ( static_cast<int64_t>(P.first) == tspace );

				libmaus2::lcs::AlignmentStatistics AS = libmaus2::lcs::AlignmentTraceContainer::getAlignmentStatistics(np.ta,np.ta+P.second);

				uint64_t const e = std::floor(AS.getErrorRate() * std::numeric_limits<uint8_t>::max() + 0.5);
				assert ( e <= std::numeric_limits<uint8_t>::max() );

				V [ from / tspace ] = std::min(V[from/tspace],static_cast<uint8_t>(e));

				from += tspace;
				np.ta += P.second;
			}
		}

		for ( uint64_t i = 0; i < V.size(); ++i )
			dataOSI.put(V[i]);

		std::cerr << "[V] " << id << std::endl;
	}

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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] reads.fasta reads_cons.fasta reads.db [blockid]\n";
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
