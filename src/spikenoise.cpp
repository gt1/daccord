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
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/random/DNABaseNoiseSpiker.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <sys/time.h>
#include <config.h>

static double getDefaultSubstRate()
{
	return 1;
}

static double getDefaultInsRate()
{
	return 0.00;
}

static double getDefaultDelRate()
{
	return 0.00;
}

static double getDefaultErrorRate()
{
	return 0.01;
}

static double getDefaultErrorRateStdDev()
{
	return 0;
}

int spikenoise(libmaus2::util::ArgParser const & arg)
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	libmaus2::random::Random::setup(static_cast<time_t>((static_cast<uint64_t>(tv.tv_sec) ^ static_cast<uint64_t>(tv.tv_usec))));

	double const substrate = arg.uniqueArgPresent("s") ? arg.getParsedArg<double>("s") : getDefaultSubstRate();
	double const insrate = arg.uniqueArgPresent("i") ? arg.getParsedArg<double>("i") : getDefaultInsRate();
	double const delrate = arg.uniqueArgPresent("d") ? arg.getParsedArg<double>("d") : getDefaultDelRate();
	double const erate = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultErrorRate();
	double const eratestddev = arg.uniqueArgPresent("stddev") ? arg.getParsedArg<double>("stddev") : getDefaultErrorRateStdDev();
	bool const omitoriginal = arg.uniqueArgPresent("omitoriginal");

	libmaus2::fastx::StreamFastAReaderWrapper SFAR(std::cin);
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

	// read sequence
	while ( SFAR.getNextPatternUnlocked(pattern) )
	{
		if ( ! omitoriginal )
		{
			// output unmodified sequence
			std::cout << pattern;
		}

		// compute modified sequence + operations applied
		std::pair<std::string,std::string> const P = libmaus2::random::DNABaseNoiseSpiker::modifyAndComment(pattern.spattern,substrate,insrate,delrate,0.0 /* homopol */,erate,eratestddev);

		// set data in pattern object
		pattern.sid = P.second;
		pattern.spattern = P.first;
		pattern.pattern = pattern.spattern.c_str();

		// output modified sequence
		std::cout << pattern;
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

static std::string helpMessage(libmaus2::util::ArgParser const & /* arg */)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("s", formatRHS("substitution rate",getDefaultSubstRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("i", formatRHS("insertion rate",getDefaultInsRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("d", formatRHS("deletion rate",getDefaultDelRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("e", formatRHS("error rate",getDefaultErrorRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("stddev", formatRHS("error rate standard deviation",getDefaultErrorRateStdDev())));

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
			return spikenoise(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
