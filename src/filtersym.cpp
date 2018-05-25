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
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>

void filterBySym(libmaus2::util::ArgParser const & arg, std::string const outfn, std::string const symkillmerge)
{
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(outfn);
	bool const verbose = arg.uniqueArgPresent("verbose");

	int64_t firstaread = -1;

	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Ain(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(outfn));
		libmaus2::dazzler::align::Overlap OVL;
		if ( Ain->peekNextOverlap(OVL) )
			firstaread = OVL.aread;
	}

	if ( firstaread < 0 )
		return;

	uint64_t p = 0;

	{
		std::string const indexfn = symkillmerge + ".index";
		libmaus2::aio::InputStreamInstance ISI(indexfn);
		ISI.clear();
		ISI.seekg(firstaread * sizeof(uint64_t));
		p = libmaus2::util::NumberSerialisation::deserialiseNumber(ISI);
	}

	libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Ain(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(outfn));
	libmaus2::aio::InputStreamInstance ISIsymkill(symkillmerge);
	ISIsymkill.clear();
	ISIsymkill.seekg(p);
	libmaus2::aio::SerialisedPeeker<libmaus2::dazzler::align::OverlapInfo> AP(ISIsymkill);
	libmaus2::dazzler::align::Overlap OVL;
	libmaus2::dazzler::align::OverlapInfo info;

	std::string const tmpoutfn = outfn + ".tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(tmpoutfn);
	libmaus2::dazzler::align::AlignmentWriter::unique_ptr_type AWtmp(new libmaus2::dazzler::align::AlignmentWriter(tmpoutfn,tspace,false /* index */));

	libmaus2::dazzler::align::OverlapInfo OVLprev;
	bool OVLprevvalid = false;
	libmaus2::dazzler::align::OverlapInfo infoprev;
	bool infoprevvalid = false;

	while ( Ain->peekNextOverlap(OVL) && AP.peekNext(info) )
	{
		libmaus2::dazzler::align::OverlapInfo const OVLinfo = OVL.getHeader().getInfo();

		bool const ovlok = (!OVLprevvalid) || (OVLprev < OVLinfo);
		bool const infook = (!infoprevvalid) || (infoprev < info);

		if ( ! ovlok )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] error in overlap stream " << OVLprev << " " << OVLinfo << std::endl;
			lme.finish();
			throw lme;
		}
		if ( ! infook )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] error in filter stream " << infoprev << " " << info << std::endl;
			lme.finish();
			throw lme;
		}

		OVLprev = OVLinfo;
		OVLprevvalid = true;
		infoprev = info;
		infoprevvalid = true;

		if ( OVLinfo < info )
		{
			while ( Ain->peekNextOverlap(OVL) && OVL.getHeader().getInfo() == OVLprev )
			{
				Ain->getNextOverlap(OVL);
				AWtmp->put(OVL);
			}
			infoprevvalid = false;
		}
		else if ( info < OVLinfo )
		{
			while ( AP.peekNext(info) && info == infoprev )
			{
				AP.getNext(info);
			}
			OVLprevvalid = false;
		}
		else
		{
			assert ( info == OVLinfo );

			while ( Ain->peekNextOverlap(OVL) && OVL.getHeader().getInfo() == OVLprev )
				Ain->getNextOverlap(OVL);
			while ( AP.peekNext(info) && info == infoprev )
				AP.getNext(info);
		}

	}

	while ( Ain->getNextOverlap(OVL) )
	{
		libmaus2::dazzler::align::OverlapInfo const OVLinfo = OVL.getHeader().getInfo();
		bool const ovlok = (!OVLprevvalid) || (OVLprev < OVLinfo);

		if ( ! ovlok )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] error in overlap stream " << OVLprev << " " << OVLinfo << std::endl;
			lme.finish();
			throw lme;
		}

		OVLprev = OVLinfo;
		OVLprevvalid = true;

		AWtmp->put(OVL);
	}

	AWtmp.reset();

	libmaus2::aio::OutputStreamFactoryContainer::rename(tmpoutfn,outfn);
	libmaus2::dazzler::align::OverlapIndexer::constructIndex(outfn);
}

int filtersym(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const symfn = arg[1];
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	filterBySym(arg,outfn,symfn);

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
			std::cerr << "usage: " << arg.progname << " [options] in.las in.sym\n";
			std::cerr << std::endl;
			std::cerr << "optional parameters:" << std::endl << std::endl;
			return EXIT_SUCCESS;
		}
		else
		{
			return filtersym(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
