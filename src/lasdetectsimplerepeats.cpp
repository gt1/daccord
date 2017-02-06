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

#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/OverlapProperCheck.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>

static double getDefaultTermVal()
{
	return 0.35;
}

static uint64_t getDefaultDepthThres()
{
	return 10;
}

struct Handler
{
	std::set<uint64_t> active;
	libmaus2::util::FiniteSizeHeap< std::pair<int64_t,uint64_t> > H;

	Handler() : H(1024) {}

	void handleVector(
		std::ostream & out,
		int64_t const prevread,
		std::vector< std::pair<int64_t,int64_t> > & V,
		uint64_t thres
	)
	{
		if ( V.size() )
		{
			// sort repeat evidence
			std::sort(V.begin(),V.end());

			std::vector < std::pair<int64_t,int64_t> > IV;

			int64_t istart = 0;
			for ( uint64_t i = 0; i < V.size(); ++i )
			{
				// remove inactive repeat regions before start of V[i]
				while ( !H.empty() && H.top().first <= V[i].first )
				{
					std::pair<int64_t,uint64_t> P = H.pop();

					if ( P.first > istart )
					{
						if ( active.size() >= thres )
						{
							// std::cerr << prevread << " " << istart << "," << P.first << "," << active.size() << std::endl;

							if ( IV.size() && IV.back().second == istart )
								IV.back().second = P.first;
							else
								IV.push_back(std::pair<int64_t,int64_t>(istart,P.first));
						}
					}
					istart = P.first;

					active.erase(P.second);
				}

				if ( active.size() && V[i].first > istart )
				{
					if ( active.size() >= thres )
					{
						// std::cerr << prevread << " " << istart << "," << V[i].first << "," << active.size() << std::endl;

						if ( IV.size() && IV.back().second == istart )
							IV.back().second = V[i].first;
						else
							IV.push_back(std::pair<int64_t,int64_t>(istart,V[i].first));
					}
				}

				H.pushBump(std::pair<int64_t,uint64_t>(V[i].second,i));
				active.insert(i);

				istart = V[i].first;
			}

			while ( !H.empty() )
			{
				std::pair<int64_t,uint64_t> P = H.pop();

				if ( P.first > istart )
				{
					if ( active.size() >= thres )
					{
						// std::cerr << prevread << " " << istart << "," << P.first << "," << active.size() << std::endl;

						if ( IV.size() && IV.back().second == istart )
							IV.back().second = P.first;
						else
							IV.push_back(std::pair<int64_t,int64_t>(istart,P.first));
					}
				}
				istart = P.first;

				active.erase(P.second);
			}

			// write repeat regions as triples (readid,from,to)
			for ( uint64_t i = 0; i < IV.size(); ++i )
			{
				libmaus2::util::NumberSerialisation::serialiseNumber(out,prevread);
				libmaus2::util::NumberSerialisation::serialiseNumber(out,IV[i].first);
				libmaus2::util::NumberSerialisation::serialiseNumber(out,IV[i].second);
				// std::cerr << prevread << " " << IV[i].first << "," << IV[i].second << std::endl;
			}

			V.resize(0);
		}
	}
};


int lasdetectsimplerepeats(libmaus2::util::ArgParser const & arg)
{
	double const termval = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultTermVal();
	uint64_t const dthres = arg.uniqueArgPresent("d") ? arg.getParsedArg<uint64_t>("d") : getDefaultDepthThres();

	std::string const dbfn = arg[0];
	// std::string const graphfn = arg[1];

	std::vector<std::string> Vinfn;
	for ( uint64_t i = 1; i < arg.size(); ++i )
		Vinfn.push_back(arg[i]);

	libmaus2::dazzler::db::DatabaseFile DB(dbfn);
	if ( DB.part != 0 )
	{
		std::cerr << "Partial databases are not supported." << std::endl;
		return EXIT_FAILURE;
	}

	// trim database
	DB.computeTrimVector();

	// get read length vector
	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	// load inqual track
	libmaus2::dazzler::db::Track::unique_ptr_type Ptrack(DB.readTrack("inqual",0));
	// overlap checker class
	libmaus2::dazzler::align::OverlapProperCheck OPC(RL,*Ptrack,termval);
	// repeat handler
	Handler H;

	// output stream
	std::ostream & out = std::cout;

	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		std::string const infn = Vinfn[i];

		// open LAS file
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type
			Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(infn));
		// get tspace parameter
		int64_t const tspace = Plas->Palgn->tspace;

		// overlap
		libmaus2::dazzler::align::Overlap OVL;

		// repeat evidence vector
		std::vector< std::pair<int64_t,int64_t> > V;

		// previous read id
		int64_t prevread = -1;

		// read overlaps from LAS file
		for ( uint64_t c = 0 ; Plas->getNextOverlap(OVL) ; ++c )
		{
			// if this is a new A read id
			if ( OVL.aread != prevread )
			{
				// handle data for previous A read
				H.handleVector(out,prevread,V,dthres);
				if ( prevread != -1 )
					std::cerr << "[V] " << prevread << std::endl;
				prevread = OVL.aread;
			}

			// get if this is a proper overlap
			libmaus2::dazzler::align::OverlapProperCheck::OverlapProperCheckInfo const proper = OPC(OVL,tspace);

			// if improper on left and improper on right
			if ( ! proper.termleft && ! proper.termright )
			{
				// std::cerr << "repeat " << OVL.aread << " " << OVL.path.abpos << "," << OVL.path.aepos << std::endl;

				// put in evidence vector
				V.push_back(std::pair<uint64_t,uint64_t>(OVL.path.abpos,OVL.path.aepos));
			}
		}

		// handle last read
		H.handleVector(out,prevread,V,dthres);

		if ( prevread != -1 )
			std::cerr << "[V] " << prevread << std::endl;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] in.db in.las\n";
			std::cerr << std::endl;
			std::cerr << "optional parameters:" << std::endl << std::endl;
			std::cerr << " -e: error threshold for proper alignment termination (default: " << getDefaultTermVal() << ")" << std::endl;
			std::cerr << " -d: depth threshold for repeat detection (default: " << getDefaultDepthThres() << ")" << std::endl;
			return EXIT_SUCCESS;
		}
		else
		{
			return lasdetectsimplerepeats(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
