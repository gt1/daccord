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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>
#include <libmaus2/util/ArgParser.hpp>

void handleOverlaps(
	// error histogram
	libmaus2::autoarray::AutoArray < std::pair < uint64_t, double > > & AH,
	uint64_t const oh,
	// aread id
	int64_t const cura,
	// trace point distance
	int64_t const tspace,
	// read length array
	std::vector<uint64_t> const & RL,
	uint64_t const d,
	libmaus2::autoarray::AutoArray<uint64_t> & annosize,
	std::ostream & odata,
	std::pair<uint64_t,uint64_t> const & blockintv,
	uint64_t const * spanS
)
{
	std::sort(AH.begin(),AH.begin()+oh);

	uint64_t const alen = RL[cura];
	uint64_t const ablocks = (alen + tspace - 1) / tspace;

	// extract all blocks with at least d/4 other blocks mapping there
	std::vector < double > D( ablocks , 1.0 );
	for ( uint64_t i = 0; i < ablocks; ++i )
		assert ( D[i] == 1.0 );

	uint64_t const dthres = std::max(static_cast<uint64_t>(d/4),static_cast<uint64_t>(1));
	uint64_t ilow = 0;
	while ( ilow < oh )
	{
		uint64_t const off = AH[ilow].first;
		uint64_t ihigh = ilow+1;
		while ( ihigh < oh && AH[ihigh].first == off )
			++ihigh;

		uint64_t const depth = ihigh-ilow;

		if ( depth >= dthres && spanS[off] >= dthres )
		{
			double s = 0.0;
			for ( uint64_t i = 0; i < dthres; ++i )
				s += AH[ilow+i].second;
			double const erate = s / dthres;

			bool const ok = erate >= 0 && erate <= 1.0;

			if ( ! ok )
			{
				std::cerr << "s=" << s << std::endl;
				std::cerr << "dthres=" << dthres << std::endl;
				assert ( ok );
			}


			// store error rate
			D[off] = erate;
		}

		#if 0
		if ( depth >= dthres && spanS[off] < dthres )
		{
			std::cerr << "[V] cura=" << cura << " off=" << off << " depth=" << depth << " spanS[off]=" << spanS[off] << std::endl;
		}
		#endif

		ilow = ihigh;
	}

	#if 0
	if ( cura == 10797 || cura == 17917 )
	for ( uint64_t i = 0; i < D.size(); ++i )
		std::cerr << cura << "(" << i << "," << D[i] << ")" << std::endl;
	//std::cerr << std::endl;
	#endif

	// implicit quality values
	annosize [ cura - blockintv.first ] = D.size();
	uint64_t const maxval = std::numeric_limits<unsigned char>::max();
	for ( uint64_t i = 0; i < D.size(); ++i )
	{
		#if 0
		if ( D[i] > 0.4 )
			std::cerr << "aread=" << cura << " block=" << i << " erate=" << D[i] << std::endl;
		#endif

		uint64_t const e = ::std::floor(D[i] * 255.0 + 0.5);
		//std::cerr << D[i] << " " << e << std::endl;
		odata.put(std::min(e,maxval));
	}
}

int main(int argc, char *argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.size() < 2 )
		{
			std::cerr << "usage: " << argv[0] << " -d<depth> <reads.db> <alignments.las>" << std::endl;
			return EXIT_FAILURE;
		}

		std::string const dbfn = arg[0];

		if ( ! arg.uniqueArgPresent("d") )
		{
			std::cerr << "[V] argument -d required" << std::endl;
			return EXIT_FAILURE;
		}

		uint64_t const d = arg.getUnsignedNumericArg<uint64_t>("d");

		libmaus2::dazzler::db::DatabaseFile DB(dbfn);
		DB.computeTrimVector();

		if ( DB.part != 0 )
		{
			std::cerr << "Partial databases are not supported." << std::endl;
			return EXIT_FAILURE;
		}

		// read all meta data
		std::vector<uint64_t> RL;
		DB.getAllReadLengths(RL);

		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(
			std::vector<std::string>(arg.restargs.begin()+1,arg.restargs.end())
		);

		for ( uint64_t a = 1; a < arg.size(); ++a )
		{
			std::string const aligns = arg[a];

			std::cerr << "[V] processing " << aligns << std::endl;

			int64_t blockid = 0;
			if ( aligns.find('.') != std::string::npos )
			{
				std::string const suffix = aligns.substr(aligns.find_first_of('.')+1);
				if ( suffix.size() && isdigit(suffix[0]) )
				{
					char const * p = suffix.c_str();
					char const * pe = p;

					while ( *pe && isdigit(*pe) )
					{
						char const d = *pe;
						int64_t const digit = d - '0';
						blockid *= 10;
						blockid += digit;
						++pe;
					}
				}
			}

			std::pair<uint64_t,uint64_t> const blockintv = DB.getTrimmedBlockInterval(blockid);

			std::string const annofilename = DB.getBlockTrackAnnoFileName("inqual",blockid);
			std::string const datafilename = DB.getBlockTrackDataFileName("inqual",blockid);

			libmaus2::autoarray::AutoArray<uint64_t> annosize(blockintv.second-blockintv.first,true);

			libmaus2::aio::OutputStreamInstance annofile(annofilename);
			libmaus2::aio::OutputStreamInstance datafile(datafilename);

			// open alignment file
			libmaus2::aio::InputStream::unique_ptr_type Palgnfile(libmaus2::aio::InputStreamFactoryContainer::constructUnique(aligns));
			libmaus2::dazzler::align::AlignmentFile algn(*Palgnfile);

			libmaus2::dazzler::align::Overlap OVL;

			// current a-read id
			int64_t preva = std::numeric_limits<int64_t>::min();
			int64_t cura = -1;
			int64_t nextexptd = blockintv.first;

			uint64_t o = 0;
			libmaus2::autoarray::AutoArray < std::pair < uint64_t, double > > A;

			int64_t lp = 0;

			libmaus2::autoarray::AutoArray < bool > spanA;
			libmaus2::autoarray::AutoArray < uint64_t > spanS;

			// get next overlap
			while ( algn.getNextOverlap(*Palgnfile,OVL) )
			{
				if ( OVL.aread < preva )
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "ids of a reads are not increasing: OVL.aread=" << OVL.aread << " preva=" << preva << std::endl;
					lme.finish();
					throw lme;
				}
				preva = OVL.aread;

				if ( OVL.aread != cura && cura >= 0 )
				{
					while ( nextexptd < cura )
					{
						// std::cerr << "[V] filling empty " << nextexptd << std::endl;
						handleOverlaps(A,0,nextexptd++,tspace,RL,d,annosize,datafile,blockintv,spanS.begin());
					}

					assert ( nextexptd == cura );
					handleOverlaps(A,o,cura,tspace,RL,d,annosize,datafile,blockintv,spanS.begin());

					#if 0
					if ( cura == 10797 || cura == 17917 )
					for ( uint64_t i = 0; i < spanS.size(); ++i )
					{
						std::cerr << "cura=" << cura << " i=" << i << " spanS[i]=" << spanS[i] << std::endl;
					}
					#endif

					nextexptd++;

					if ( cura/1024 != lp/1024 )
					{
						lp = cura;
						std::cerr << "[V] " << lp << std::endl;
					}
				}

				if ( OVL.aread != cura )
				{
					o = 0;
					int64_t const numnewblocks = (RL[OVL.aread] + tspace - 1)/tspace;
					spanA.resize(numnewblocks);
					std::fill(spanA.begin(),spanA.end(),0ull);
					spanS.resize(numnewblocks);
					std::fill(spanS.begin(),spanS.end(),0ull);
				}

				cura = OVL.aread;
				o = OVL.fillErrorHistogram(tspace,A,o,RL[OVL.aread]);
				OVL.fillSpanHistogram(tspace,RL[OVL.aread],0.3 /* ethres */,1 /* bthres */,spanA,spanS);
			}

			if ( o )
			{
				while ( nextexptd < cura )
				{
					// std::cerr << "[V] filling empty " << nextexptd << std::endl;
					handleOverlaps(A,0,nextexptd++,tspace,RL,d,annosize,datafile,blockintv,spanS.begin());
				}

				assert ( nextexptd == cura );
				handleOverlaps(A,o,cura,tspace,RL,d,annosize,datafile,blockintv,spanS.begin());
				o = 0;
				nextexptd++;
			}

			while ( nextexptd < static_cast<int64_t>(blockintv.second) )
			{
				// std::cerr << "[V] filling empty " << nextexptd << std::endl;
				handleOverlaps(A,0,nextexptd++,tspace,RL,d,annosize,datafile,blockintv,spanS.begin());
			}
			assert ( nextexptd == static_cast<int64_t>(blockintv.second) );

			// flush and close inqual data file
			datafile.flush();

			// write inqual anno file
			uint64_t annooff = 0;
			libmaus2::dazzler::db::OutputBase::putLittleEndianInteger4(annofile,annosize.size() /* tracklen */,annooff);
			libmaus2::dazzler::db::OutputBase::putLittleEndianInteger4(annofile,8 /* size of pointer */,annooff);
			uint64_t s = 0;
			for ( uint64_t i = 0; i < annosize.size(); ++i )
			{
				libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(annofile,s,annooff);
				s += annosize[i];
			}
			libmaus2::dazzler::db::OutputBase::putLittleEndianInteger8(annofile,s,annooff);
			annofile.flush();
		}

		return EXIT_SUCCESS;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
