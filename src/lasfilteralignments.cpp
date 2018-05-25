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

#include <libmaus2/dazzler/align/OverlapInfoIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/OverlapProperCheck.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/db/InqualContainer.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/dazzler/align/SortingOverlapOutputBuffer.hpp>

static double getDefaultTermVal()
{
	return 0.35;
}

struct NamedInterval : std::pair<int64_t,int64_t>
{
	std::string name;

	NamedInterval() {}
	NamedInterval(
		std::pair<int64_t,int64_t> const & rI,
		std::string const & rname
	) : std::pair<int64_t,int64_t>(rI), name(rname) {}
};


int lasfilteralignments(libmaus2::util::ArgParser const & arg)
{
	double const termval = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultTermVal();
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	std::string const outfn = arg[0];
	std::string const symkillfn = outfn + ".sym";
	std::string const dbfn = arg[1];

	std::vector<std::string> Vinfn;
	for ( uint64_t i = 2; i < arg.size(); ++i )
		Vinfn.push_back(arg[i]);

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vinfn);

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	if ( DB.part != 0 )
	{
		std::cerr << "Partial databases are not supported." << std::endl;
		return EXIT_FAILURE;
	}

	DB.computeTrimVector();

	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	libmaus2::dazzler::db::Track::unique_ptr_type Ptrack(DB.readTrack("inqual",0));
	libmaus2::dazzler::align::OverlapProperCheck OPC(RL,*Ptrack,termval);
	libmaus2::aio::OutputStreamInstance::unique_ptr_type Psymkill(new libmaus2::aio::OutputStreamInstance(symkillfn));

	std::vector < NamedInterval > VNI;

	libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AW(new libmaus2::dazzler::align::AlignmentWriter(outfn,tspace,false /* index */));

	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		std::string const infn = Vinfn[i];
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));

		libmaus2::dazzler::align::Overlap OVL;

		int64_t prevread = -1;
		uint64_t kept = 0;
		uint64_t removed = 0;

		for ( uint64_t c = 0 ; Plas->getNextOverlap(OVL) ; ++c )
		{
			if ( OVL.aread != prevread )
			{
				if ( prevread >= 0 )
				{
					if ( prevread > OVL.aread )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "[E] file " << infn << " is unsorted" << std::endl;
						lme.finish();
						throw lme;
					}
				}

				if ( kept || removed )
					std::cerr << "[V] " << prevread << " kept " <<kept << " removed " << removed << std::endl;
				kept = 0;
				removed = 0;
				prevread = OVL.aread;
			}

			bool const keep = OPC(OVL,tspace).proper;

			if ( keep )
			{
				AW->put(OVL);
				++kept;
			}
			else
			{
				libmaus2::dazzler::align::OverlapInfo info = OVL.getHeader().getInfo().swapped();

				if ( OVL.isInverse() )
				{
					uint64_t const alen = (*(OPC.RL))[info.aread >> 1];
					uint64_t const blen = (*(OPC.RL))[info.bread >> 1];
					info = info.inverse(alen,blen);
				}

				assert ( (info.aread & 1) == 0 );

				info.serialise(*Psymkill);

				++removed;
			}
		}

		if ( kept || removed )
			std::cerr << "[V] " << prevread << " kept " <<kept << " removed " << removed << std::endl;
	}

	AW.reset();

	// sort output
	libmaus2::dazzler::align::SortingOverlapOutputBuffer<
		libmaus2::dazzler::align::OverlapFullComparator
	>::sort(outfn,tmpfilebase);

	Psymkill->flush();
	Psymkill.reset();

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::dazzler::align::OverlapInfo>::sort(symkillfn,16*1024*1024);
	libmaus2::dazzler::align::OverlapInfoIndexer::createInfoIndex(symkillfn,DB.size());
	libmaus2::dazzler::align::OverlapIndexer::constructIndex(outfn);

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
			std::cerr << "usage: " << arg.progname << " [options] out.las in.db in1.las ...\n";
			std::cerr << std::endl;
			std::cerr << "optional parameters:" << std::endl << std::endl;
			std::cerr << " -e: error threshold for proper alignment termination (default: " << getDefaultTermVal() << ")" << std::endl;
			return EXIT_SUCCESS;
		}
		else
		{
			return lasfilteralignments(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
