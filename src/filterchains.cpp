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
#include <config.h>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/timing/RealTimeClock.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <ChainSet.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/FiniteSizeHeap.hpp>

uint64_t getDefaultMinLength()
{
	return 4000;
}

void patternToUpper(libmaus2::fastx::Pattern & pattern)
{
	for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
		pattern.spattern[i] = ::toupper(pattern.spattern[i]);

	pattern.pattern = pattern.spattern.c_str();
}

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

static std::string getFastAIndexFileName(std::string const & consfn)
{
	std::string const consfnindex = consfn + ".phaser.index";
	return consfnindex;
}


int generateindex(
	libmaus2::util::ArgParser const & arg,
	std::string const & consfn,
	uint64_t const fastaindexmod
)
{

	std::string const consfnindex = getFastAIndexFileName(consfn);
	std::string const tmpfilebase = getTmpFileBase(arg);

	if (
		! libmaus2::util::GetFileSize::fileExists(consfnindex)
		||
		libmaus2::util::GetFileSize::isOlder(consfnindex,consfn)
	)
	{
		std::string const tmpfn = tmpfilebase + ".faindextmp";
		libmaus2::fastx::FastAReader::enumerateOffsets(consfn,tmpfn,fastaindexmod);
		libmaus2::aio::OutputStreamFactoryContainer::rename(tmpfn,consfnindex);
	}

	return EXIT_SUCCESS;
}

struct ReadAccessor
{
	std::string const consfn;
	std::string const consfnindex;
	uint64_t const mod;

	ReadAccessor(std::string const & rconsfn, std::string const & rconsfnindex, uint64_t const rmod)
	: consfn(rconsfn), consfnindex(rconsfnindex), mod(rmod)
	{

	}

	libmaus2::fastx::FastAReader::pattern_type operator[](uint64_t const i)
	{
		uint64_t const b = i / mod;
		libmaus2::aio::InputStreamInstance IISI(consfnindex);
		IISI.seekg(b * sizeof(libmaus2::fastx::FastInterval));
		libmaus2::fastx::FastInterval FI(IISI);
		libmaus2::fastx::FastAReader FARE(consfn,FI);
		uint64_t o = i - b * mod;
		for ( uint64_t j = 0; j < o; ++j )
			FARE.skipPattern();
		libmaus2::fastx::FastAReader::pattern_type P;
		bool const ok = FARE.getNextPatternUnlocked(P);

		if ( ! ok )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "ReadAccessor::operator[]: failed to get read number " << i << std::endl;
			lme.finish();
			throw lme;
		}

		patternToUpper(P);
		return P;
	}

	std::string getName(uint64_t const i)
	{
		libmaus2::fastx::FastAReader::pattern_type const P = (*this)[i];
		return P.sid;
	}

	static uint64_t getReadId(std::string name)
	{
		if ( name.find('/') != std::string::npos )
			name = name.substr(0,name.find('/'));

		std::istringstream istr(name);
		uint64_t id;
		istr >> id;

		assert ( istr && istr.peek() == std::istream::traits_type::eof() );

		assert ( id );

		return id-1;
	}

	uint64_t getReadId(uint64_t const i)
	{
		std::string name = getName(i);
		return getReadId(name);
	}
};

uint64_t getId(std::string sid, libmaus2::fastx::FastAReader::pattern_type const & pattern)
{
	if ( sid.find('/') == std::string::npos )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] unparsable read name (no slash) " << sid << std::endl;
		lme.finish();
		throw lme;
	}

	sid = sid.substr(0,sid.find('/'));

	for ( uint64_t i = 0; i < sid.size(); ++i )
		if ( ! ::isdigit(sid[i]) )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unparsable read name (string " << sid << " before slash not numerical) " << pattern.sid << std::endl;
			lme.finish();
			throw lme;
		}

	std::istringstream istr(sid);
	uint64_t id;
	istr >> id;

	assert ( istr.peek() == std::istream::traits_type::eof() );

	assert ( id );

	return id - 1;
}

std::vector<uint64_t> countConsensus(std::string const & consin)
{
	libmaus2::fastx::FastAReader FA(consin);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	std::vector<uint64_t> V;

	while ( FA.getNextPatternUnlocked(pattern) )
	{
		// >2/24/500_15948 A=[500,17349]

		std::string sid = pattern.sid;

		uint64_t const id = getId(sid,pattern);

		while ( ! (id<V.size()) )
			V.push_back(0);
		assert ( id < V.size() );

		V[id]++;
	}

	V.push_back(0);

	::libmaus2::util::PrefixSums::prefixSums(V.begin(),V.end());

	return V;
}

struct ReadFragment
{
	uint64_t from;
	uint64_t to;
	std::string fragment;

	ReadFragment() {}
	ReadFragment(uint64_t const rfrom, uint64_t const rto, std::string const & rfragment)
	: from(rfrom), to(rto), fragment(rfragment)
	{

	}
};

struct FragmentContainer
{
	int64_t cur;
	std::vector<uint64_t> const & VCNT;
	std::string const & consin;
	std::string const & consfnindex;
	uint64_t const fastaindexmod;
	std::vector < libmaus2::fastx::FastAReader::pattern_type > Vpat;
	std::vector < ReadFragment > fragments;
	std::vector < libmaus2::math::IntegerInterval<int64_t> > VIV;

	FragmentContainer(
		std::vector<uint64_t> const & rVCNT,
		std::string const & rconsin,
		std::string const & rconsfnindex,
		uint64_t const rfastaindexmod
	) : cur(-1), VCNT(rVCNT), consin(rconsin), consfnindex(rconsfnindex), fastaindexmod(rfastaindexmod) {}

	void load(int64_t const aread)
	{
		if ( aread != cur && aread + 1 < static_cast<int64_t>(VCNT.size()) )
		{
			// std::cerr << "z=" << z << " " << VCNT[z] << " " << VCNT[z+1] << std::endl;

			uint64_t const low = VCNT[aread];
			uint64_t const high = VCNT[aread+1];
			ReadAccessor RA(consin, consfnindex, fastaindexmod);

			Vpat.resize(0);
			fragments.resize(0);
			VIV.resize(0);
			for ( uint64_t i = low; i < high; ++i )
			{
				Vpat.push_back(RA[i]);

				assert (
					static_cast<int64_t>(getId(
						Vpat.back().sid,
						Vpat.back()
					))
					==
					aread
				);

				// std::cerr << "got " << Vpat.back().sid << " for " << z << std::endl;
			}

			for ( uint64_t i = 0; i < Vpat.size(); ++i )
			{
				std::string sid = Vpat[i].sid;

				assert ( sid.find(" A=[") != std::string::npos );
				sid = sid.substr(sid.find(" A=[") + strlen(" A=["));

				std::istringstream istr(sid);
				uint64_t first = 0;
				istr >> first;
				assert ( istr && istr.peek() == ',' );
				istr.get();
				uint64_t last = 0;
				istr >> last;
				assert ( istr && istr.peek() == ']' );

				std::string const SCO = Vpat[i].spattern;

				ReadFragment RF(first,last+1,SCO);

				fragments.push_back(RF);

				VIV.push_back(libmaus2::math::IntegerInterval<int64_t>(fragments.back().from,fragments.back().to-1));
			}

			cur = aread;
		}
	}
};

struct OverlapPosComparator
{
	bool operator()(libmaus2::dazzler::align::Overlap const & lhs, libmaus2::dazzler::align::Overlap const & rhs) const
	{
		return lhs.path.abpos < rhs.path.abpos;
	}
};

struct DepthLine
{
	uint64_t from;
	uint64_t to;
	int64_t d;

	DepthLine() {}
	DepthLine(uint64_t const rfrom, uint64_t const rto, int64_t const rd) : from(rfrom), to(rto), d(rd) {}
};

std::ostream & operator<<(std::ostream & out, DepthLine const & D)
{
	return out << "DepthLine(" << D.from << "," << D.to << "," << D.d << ")";
}

void handle(libmaus2::dazzler::align::Overlap * A, uint64_t const f, libmaus2::dazzler::align::AlignmentWriter & AW, int64_t const /* dthres */)
{
	std::sort(A,A+f,libmaus2::dazzler::align::OverlapFullComparator());
	for ( uint64_t i = 0; i < f; ++i )
		AW.put(A[i]);
}

int filterchains(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const consin = arg[1];
	int64_t const dthres = 10;

	uint64_t const fastaindexmod = 1;
	generateindex(arg,consin,fastaindexmod);
	std::vector<uint64_t> VCNT = countConsensus(consin);
	std::string const consfnindex = getFastAIndexFileName(consin);

	std::vector<std::string> Vin;
	for ( uint64_t i = 2; i < arg.size(); ++i )
		Vin.push_back(arg[i]);
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vin);
	libmaus2::dazzler::align::AlignmentWriter AW(outfn,tspace);
	std::pair<int64_t,int64_t> prevaread(-1,-1);
	libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > A;
	libmaus2::autoarray::AutoArray < libmaus2::dazzler::align::Overlap > B;
	uint64_t fb = 0;
	int64_t idb = -1;
	uint64_t const minlength = arg.uniqueArgPresent("l") ? arg.getUnsignedNumericArg<uint64_t>("l") : getDefaultMinLength();

	FragmentContainer FC(VCNT,consin,consfnindex,fastaindexmod);

	for ( uint64_t z = 0; z < Vin.size(); ++z )
	{
		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type AFR(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(Vin[z]));

		libmaus2::dazzler::align::Overlap refOVL, OVL;
		while ( AFR->peekNextOverlap(refOVL) )
		{
			assert ( (refOVL.aread > prevaread.first) || (refOVL.aread == prevaread.first && refOVL.bread > prevaread.second)  );

			uint64_t f = 0;
			while ( AFR->peekNextOverlap(OVL) && OVL.aread == refOVL.aread && OVL.bread == refOVL.bread )
			{
				AFR->getNextOverlap(OVL);
				A.push(f,OVL);
			}

			prevaread = std::pair<int64_t,int64_t>(refOVL.aread,refOVL.bread);

			if ( prevaread.first != idb )
			{
				if ( fb )
					handle(B.begin(),fb,AW,dthres);
				idb = prevaread.first;
				fb = 0;
			}

			FC.load(refOVL.aread);
			ChainSet CS(A.begin(),f);

			for ( uint64_t chainid = 0; chainid < CS.size(); ++chainid )
			{
				uint64_t s = 0;

				for ( uint64_t chainsubid = 0; chainsubid < CS.size(chainid); ++chainsubid )
				{
					uint64_t const i = CS(chainid,chainsubid);
					libmaus2::math::IntegerInterval<int64_t> I(A[i].path.abpos,A[i].path.aepos-1);

					for ( uint64_t i = 0; i < FC.VIV.size(); ++i )
						s += I.intersection(FC.VIV[i]).diameter();
				}

				// std::cerr << prevaread.first << "," << prevaread.second << " s=" << s << std::endl;

				if ( s >= minlength )
				{
					for ( uint64_t chainsubid = 0; chainsubid < CS.size(chainid); ++chainsubid )
					{
						uint64_t const i = CS(chainid,chainsubid);
						B.push(fb,A[i]);
						// AW.put(A[i]);
					}
				}
			}

			// std::cerr << prevaread.first << "," << prevaread.second << " " << f << std::endl;
		}
	};

	if ( fb )
		handle(B.begin(),fb,AW,dthres);

	return EXIT_SUCCESS;
}


template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

/*
 parameters:
 */
static std::string helpMessage(libmaus2::util::ArgParser const & /* arg */)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	// optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("l", formatRHS("minimum chain length",getDefaultMinLength())));

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
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 1 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las cons.fasta in1.las ...\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			libmaus2::timing::RealTimeClock rtc;
			rtc.start();

			int r = EXIT_FAILURE;

			r = filterchains(arg);

			std::cerr << "[V] processing time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

			return r;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << "[E] " << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
