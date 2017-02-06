/*
    consensus
    Copyright (C) 2016 German Tischler

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
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/bambam/StrCmpNum.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>

struct FastAEntry
{
	std::string id;
	std::string data;

	FastAEntry() {}
	FastAEntry(std::string const & rid, std::string const & rdata) : id(rid), data(rdata) {}
	FastAEntry(std::istream & in)
	:
		id(libmaus2::util::StringSerialisation::deserialiseString(in)),
		data(libmaus2::util::StringSerialisation::deserialiseString(in))
	{}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::StringSerialisation::serialiseString(out,id);
		libmaus2::util::StringSerialisation::serialiseString(out,data);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		*this = FastAEntry(in);
		return in;
	}

	bool operator<(FastAEntry const & O) const
	{
		return libmaus2::bambam::StrCmpNum::strcmpnum(id.c_str(),O.id.c_str()) < 0;
	}
};

int sortfasta(libmaus2::util::ArgParser const & arg)
{
	libmaus2::fastx::FastAReader::pattern_type fullpat;
	libmaus2::fastx::StreamFastAReaderWrapper fullread(std::cin);
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	std::string const sortfile = tmpfilebase + ".sort";
	libmaus2::util::TempFileRemovalContainer::addTempFile(sortfile);
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<FastAEntry> SBO(sortfile);
	while ( fullread.getNextPatternUnlocked(fullpat) )
		SBO.put(FastAEntry(fullpat.sid,fullpat.spattern));
	libmaus2::sorting::SerialisingSortingBufferedOutputFile<FastAEntry>::merger_ptr_type Pmerger(SBO.getMerger());
	FastAEntry O;
	while ( Pmerger->getNext(O) )
	{
		std::cout << '>' << O.id << "\n";
		std::cout << O.data << "\n";
	}

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


int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] <in.fasta >out.fasta\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return sortfasta(arg);
		}


	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
