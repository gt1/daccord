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
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <config.h>

uint64_t getId(libmaus2::fastx::FastAReader::pattern_type const & pat)
{
	std::string sid = pat.getShortStringId();
	assert ( sid.find('_') != std::string::npos );
	sid = sid.substr(sid.find('_')+1);
	std::istringstream istr(sid);
	uint64_t id;
	istr >> id;
	assert ( istr && istr.peek() == std::istream::traits_type::eof() );
	return id;
}

int fillfasta(libmaus2::util::ArgParser const & arg)
{
	std::string const fullfn = arg[0];
	std::string const partfn = arg[1];

	libmaus2::fastx::FastAReader::pattern_type fullpat;
	libmaus2::fastx::FastAReader::pattern_type partpat;

	libmaus2::fastx::FastaPeeker fullpeek(fullfn);
	libmaus2::fastx::FastaPeeker partpeek(partfn);

	while ( fullpeek.peekNext(fullpat) && partpeek.peekNext(partpat) )
	{
		uint64_t const fullid = getId(fullpat);
		uint64_t const partid = getId(partpat);

		if ( fullid < partid )
		{
			fullpeek.getNext(fullpat);

			for ( uint64_t i = 0; i < fullpat.spattern.size(); ++i )
				fullpat.spattern[i] = ::tolower(fullpat.spattern[i]);
			fullpat.pattern = fullpat.spattern.c_str();

			fullpat.printMultiLine(std::cout,80);
		}
		else if ( partid < fullid )
		{
			partpeek.getNext(partpat);
			assert ( false );
		}
		else
		{
			assert ( partid == fullid );
			fullpeek.getNext(fullpat);
			partpeek.getNext(partpat);
			partpat.printMultiLine(std::cout,80);
		}
	}

	assert ( ! partpeek.peekNext(partpat) );

	while ( fullpeek.getNext(fullpat) )
	{
		for ( uint64_t i = 0; i < fullpat.spattern.size(); ++i )
			fullpat.spattern[i] = ::tolower(fullpat.spattern[i]);
		fullpat.pattern = fullpat.spattern.c_str();

		fullpat.printMultiLine(std::cout,80);
	}

	return EXIT_SUCCESS;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " uncorrected.fasta corrected_mapped.fasta\n";
			return EXIT_SUCCESS;
		}
		else
		{
			return fillfasta(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
