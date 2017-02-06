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

#include <libmaus2/dazzler/align/TrueOverlap.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/bambam/parallel/FragmentAlignmentBufferFragment.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

int bamidrename(libmaus2::util::ArgParser const & /* arg */, libmaus2::util::ArgInfo const & arginfo)
{
	libmaus2::bambam::BamDecoder dec(std::cin);
	libmaus2::bambam::BamAlignment & algn = dec.getAlignment();
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(libmaus2::bambam::BamBlockWriterBaseFactory::construct(dec.getHeader(),arginfo, 0));
	uint64_t id = 0;

	while ( dec.readAlignment() )
	{
		std::ostringstream ostr;

		ostr << "L0/" << id++ << "/0_" << algn.getLseq();

		std::string const newname = ostr.str();

		algn.replaceName(newname.begin(),newname.size());

		Pwriter->writeAlignment(algn);
	}

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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] <in.bam >out.bam\n";
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::util::ArgInfo const arginfo(argc,argv);
			return bamidrename(arg,arginfo);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
