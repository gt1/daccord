/*
    daccord
    Copyright (C) 2018 German Tischler-Höhle

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
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/lcs/NPLinMem.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/parallel/NumCpus.hpp>

int64_t getRawId(libmaus2::fastx::FastAReader::pattern_type const & P)
{
	std::string s = P.getShortStringId();
	std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(s,std::string("/"));

	if ( 1 < tokens.size() )
	{
		std::istringstream istr(tokens[1]);
		int64_t i;
		istr >> i;

		if ( istr && istr.peek() == std::istream::traits_type::eof() )
			return i;
	}

	std::cerr << "[W] unable to parse " << s << std::endl;

	return -1;
}

int64_t getConsId(libmaus2::fastx::FastAReader::pattern_type const & P)
{
	std::string s = P.getShortStringId();
	std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(s,std::string("/"));

	if ( 0 < tokens.size() )
	{
		std::istringstream istr(tokens[0]);
		int64_t i;
		istr >> i;

		if ( istr && istr.peek() == std::istream::traits_type::eof() )
			return i-1;
	}

	std::cerr << "[W] unable to parse " << s << std::endl;

	return -1;
}

std::vector < std::string > stok(std::string const & s)
{
	std::vector < std::string > V;
	uint64_t low = 0;

	while ( low < s.size() )
	{
		while ( low < s.size() && ::isspace(s[low]) )
			++low;

		uint64_t start = low;
		while ( low < s.size() && !::isspace(s[low]) )
			++low;

		if ( low > start )
			V.push_back(s.substr(start,low-start));
	}

	return V;
}

libmaus2::math::IntegerInterval<int64_t> getInterval(libmaus2::fastx::FastAReader::pattern_type const & P)
{
	std::vector < std::string > VT = stok(P.sid);

	if ( 1 < VT.size() )
	{
		std::string A = VT[1];

		if ( A.size() >= 3 && A[0] == 'A' && A[1] == '=' && A[2] == '[' && A[A.size()-1] == ']' )
		{
			A = A.substr(3);
			A = A.substr(0,A.size()-1);

			std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize(A,std::string(","));

			if ( tokens.size() == 2 )
			{
				std::istringstream astr(tokens[0]);
				std::istringstream bstr(tokens[1]);
				int64_t ia;
				int64_t ib;
				astr >> ia;
				bstr >> ib;

				if (
					astr && astr.peek() == std::istream::traits_type::eof() &&
					bstr && bstr.peek() == std::istream::traits_type::eof()
				)
				{
					return libmaus2::math::IntegerInterval<int64_t>(ia,ib);
				}
				else
				{
					return libmaus2::math::IntegerInterval<int64_t>::empty();
				}
			}
			else
			{
				return libmaus2::math::IntegerInterval<int64_t>::empty();
			}
		}
		else
		{
			return libmaus2::math::IntegerInterval<int64_t>::empty();
		}
	}
	else
	{
		return libmaus2::math::IntegerInterval<int64_t>::empty();
	}
	// A=[0,6552]
}

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const outfn = arg[0];
		std::string const dbfn = arg[1];
		std::string const rawfast = arg[2];
		std::string const consfast = arg[3];
		std::string const lasin = arg[4];

		// number of threads
		uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

		libmaus2::dazzler::db::DatabaseFile DB(dbfn);
		DB.computeTrimVector();

		int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(lasin);

		libmaus2::fastx::FastaPeeker::unique_ptr_type praw(new libmaus2::fastx::FastaPeeker(rawfast));
		libmaus2::fastx::FastaPeeker::unique_ptr_type pcons(new libmaus2::fastx::FastaPeeker(consfast));

		libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type plasin(
			libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(lasin)
		);

		libmaus2::dazzler::align::AlignmentWriter AW(outfn,tspace);

		libmaus2::fastx::FastAReader::pattern_type pa, pb;
		while (
			praw->peekNext(pa)
			&&
			pcons->peekNext(pb)
		)
		{
			if ( getRawId(pa) < getConsId(pb) )
			{
				praw->getNext(pa);
			}
			else if ( getConsId(pb) < getRawId(pa) )
			{
				pcons->getNext(pb);
			}
			else
			{
				int64_t const id = getRawId(pa);

				libmaus2::dazzler::align::Overlap OVL;
				while ( plasin->peekNextOverlap(OVL) && OVL.aread < id )
					plasin->getNextOverlap(OVL);

				std::vector < libmaus2::dazzler::align::Overlap > VOVL;
				while ( plasin->peekNextOverlap(OVL) && OVL.aread == id )
				{
					plasin->getNextOverlap(OVL);
					VOVL.push_back(OVL);
				}


				std::vector < libmaus2::fastx::FastAReader::pattern_type > VRAW;
				std::vector < libmaus2::fastx::FastAReader::pattern_type > VCONS;
				std::vector < libmaus2::math::IntegerInterval<int64_t> > VI;

				while ( praw->peekNext(pa) && getRawId(pa) == id )
				{
					praw->getNext(pa);
					VRAW.push_back(pa);
				}
				while ( pcons->peekNext(pb) && getConsId(pb) == id )
				{
					pcons->getNext(pb);
					VCONS.push_back(pb);
					VI.push_back(getInterval(pb));
				}

				assert ( VRAW.size() == 1 );

				int64_t const l = VRAW[0].spattern.size();

				libmaus2::math::IntegerInterval<int64_t> IF(0,l-1);

				std::vector< libmaus2::math::IntegerInterval<int64_t> > VC = libmaus2::math::IntegerInterval<int64_t>::difference(IF,VI);

				struct MarkedInterval
				{
					bool raw;
					libmaus2::math::IntegerInterval<int64_t> I;
					uint64_t id;

					MarkedInterval() {}
					MarkedInterval(bool const rraw, libmaus2::math::IntegerInterval<int64_t> const &  rI, uint64_t const rid) : raw(rraw), I(rI), id(rid) {}

					bool operator<(MarkedInterval const & M) const
					{
						return I < M.I;
					}
				};

				std::vector < MarkedInterval > VR;
				for ( uint64_t i = 0; i < VI.size(); ++i )
					VR.push_back(MarkedInterval(false,VI[i],i));
				for ( uint64_t i = 0; i < VC.size(); ++i )
					VR.push_back(MarkedInterval(true,VC[i],i));

				std::sort(VR.begin(),VR.end());

				std::ostringstream basestr;
				for ( uint64_t i = 0; i < VR.size(); ++i )
				{
					if ( VR[i].raw )
						basestr << VRAW[0].spattern.substr(VR[i].I.from,VR[i].I.diameter());
					else
						basestr << VCONS[VR[i].id].spattern;
				}

				std::string const sraw = VRAW[0].spattern;
				std::string const spatch = basestr.str();

				libmaus2::lcs::NPLinMem np;
				// A raw
				// B patch
				np.np(sraw.begin(),sraw.end(),spatch.begin(),spatch.end());

				std::cerr << id << " " << np.getAlignmentStatistics() << " " << VOVL.size() << std::endl;

				std::cout << ">" << VRAW[0].sid << "\n";
				uint64_t p = 0;
				while ( p < spatch.size() )
				{
					uint64_t const cols = 80;
					uint64_t const l = spatch.size();
					uint64_t const rest = l - p;
					uint64_t const toprint = std::min(rest,cols);
					std::cout << spatch.substr(p,toprint) << "\n";
					p += toprint;
				}

				#if defined(_OPENMP)
				#pragma omp parallel for num_threads(numthreads) schedule(dynamic,1)
				#endif
				for ( uint64_t i = 0; i < VOVL.size(); ++i )
				{
					libmaus2::dazzler::align::Overlap & OVL = VOVL[i];

					int64_t abpos = OVL.path.abpos;
					int64_t aepos = OVL.path.aepos;

					std::pair<uint64_t,uint64_t> const adva = libmaus2::lcs::AlignmentTraceContainer::advanceMaxA(np.ta,np.te,abpos);
					std::pair<int64_t,int64_t> SLA = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(np.ta,np.ta + adva.second);
					std::pair<uint64_t,uint64_t> const advb = libmaus2::lcs::AlignmentTraceContainer::advanceMaxA(np.ta,np.te,aepos);
					std::pair<int64_t,int64_t> SLB = libmaus2::lcs::AlignmentTraceContainer::getStringLengthUsed(np.ta,np.ta + advb.second);

					OVL.path.abpos = SLA.second;
					OVL.path.aepos = SLB.second;

					libmaus2::lcs::NPLinMem nploc;
					std::string const b = DB.decodeRead(OVL.bread,OVL.isInverse());
					nploc.np(
						spatch.begin() + OVL.path.abpos,
						spatch.begin() + OVL.path.aepos,
						b.begin() + OVL.path.bbpos,
						b.begin() + OVL.path.bepos
					);


					libmaus2::dazzler::align::Overlap const NOVL = libmaus2::dazzler::align::Overlap::computeOverlap(
						OVL.flags,
						OVL.aread, // ref
						OVL.bread, // read
						OVL.path.abpos,
						OVL.path.aepos,
						OVL.path.bbpos,
						OVL.path.bepos,
						tspace,
						nploc
					);

					// std::cerr << nploc.getAlignmentStatistics() << std::endl;

					VOVL[i] = NOVL;
				}

				{
					libmaus2::lcs::NPLinMem nploc;
					std::string const b = DB.decodeRead(id,false);
					nploc.np(
						spatch.begin(),
						spatch.end(),
						b.begin(),
						b.end()
					);
					libmaus2::dazzler::align::Overlap const NOVL = libmaus2::dazzler::align::Overlap::computeOverlap(
						libmaus2::dazzler::align::Overlap::getPrimaryFlag(),
						id, // ref
						id, // read
						0,
						spatch.size(),
						0,
						b.size(),
						tspace,
						nploc
					);

					AW.put(NOVL);
				}

				for ( uint64_t i = 0; i < VOVL.size(); ++i )
					AW.put(VOVL[i]);
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
