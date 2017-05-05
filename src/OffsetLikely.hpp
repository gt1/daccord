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
#if ! defined(OFFSETLIKELY_HPP)
#define OFFSETLIKELY_HPP

#include <DotProduct.hpp>
#include <libmaus2/util/TempFileRemovalContainer.hpp>

struct OffsetLikely
{
	std::vector<DotProduct> DP;
	std::vector < double > dsum;
	std::vector<DotProduct> DPnorm;
	// read position [] is relevant for reference positions in given interval
	std::vector< std::pair<uint64_t,uint64_t> > Vsupport;
	std::vector<DotProduct> DPnormSquare;

	// get supporting ref position lower bound for read position i
	uint64_t getSupportLow(int64_t const i) const
	{
		return expect_true(i < static_cast<int64_t>(Vsupport.size())) ? Vsupport[i].first : DPnorm.size();
	}

	uint64_t getSupportHigh(int64_t const i) const
	{
		return expect_true(i < static_cast<int64_t>(Vsupport.size())) ? Vsupport[i].second : DPnorm.size();
	}

	uint64_t size() const
	{
		return DP.size();
	}

	void push_back(DotProduct const & D)
	{
		DP.push_back(D);
	}

	DotProduct const & operator[](uint64_t const i) const
	{
		return DP[i];
	}

	void setup()
	{
		dsum.resize(0);

		// maximum size of any dot product (max read position supported + 1)
		uint64_t maxsize = 0;
		for ( uint64_t i = 0; i < size(); ++i )
			maxsize = std::max(maxsize,DP[i].size());

		for ( uint64_t i = 0; i < maxsize; ++i )
		{
			double sum = 0.0;

			for ( uint64_t j = 0; j < size(); ++j )
				sum += (*this)[j][i];

			dsum.push_back(sum);
		}

		DPnorm = DP;
		for ( uint64_t i = 0; i < DPnorm.size(); ++i )
			for ( uint64_t j = 0; j < maxsize; ++j )
				DPnorm[i].normaliseValue(j,dsum[j]);

		uint64_t j = 0, k = 0;
		for ( uint64_t i = 0; i < maxsize; ++i )
		{
			// update lower end
			while ( j < DPnorm.size() && i >= DPnorm[j].firstsign + DPnorm[j].V.size() )
				++j;
			// update upper end
			while ( k < DPnorm.size() && DPnorm[k].firstsign <= i )
				++k;

			Vsupport.push_back(std::pair<uint64_t,uint64_t>(j,k));
		}

		DPnormSquare = DP;
		for ( uint64_t i = 0; i < DPnormSquare.size(); ++i )
			DPnormSquare[i].normalise();
	}

	void plotData(std::string const & prefix) const
	{
		std::vector<std::string> Vfn;
		std::string const gplfn = prefix + ".gpl";
		libmaus2::util::TempFileRemovalContainer::addTempFile(gplfn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type POSI(
			new libmaus2::aio::OutputStreamInstance(gplfn)
		);

		*POSI << "set terminal pdf\n";
		*POSI << "set xlabel \"position in read\"\n";
		*POSI << "set ylabel \"probability\"\n";
		*POSI << "plot [0 to 55] ";
		for ( uint64_t i = 0; i < DP.size(); ++i )
		{
			std::ostringstream fnostr;
			fnostr << prefix << "_" << i << ".gpl";
			std::string const fn = fnostr.str();
			Vfn.push_back(fn);
			libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
			libmaus2::aio::OutputStreamInstance OSI(fn);
			DP[i].printData(OSI);

			if ( i > 0 )
				*POSI << ",";

			*POSI << "\"" << fn << "\" smooth bezier title \"\"";
		}
		*POSI << "\n";
		POSI->flush();
		POSI.reset();

		std::ostringstream comstr;
		comstr << "gnuplot <" << gplfn << " >" << prefix << ".pdf";
		std::string const com = comstr.str();

		int const r = system(com.c_str());
		if ( r != EXIT_SUCCESS )
		{

		}
	}

	void plotDataSingle(std::string const & prefix, uint64_t const z) const
	{
		std::vector<std::string> Vfn;
		std::string const gplfn = prefix + ".gpl";
		libmaus2::util::TempFileRemovalContainer::addTempFile(gplfn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type POSI(
			new libmaus2::aio::OutputStreamInstance(gplfn)
		);

		*POSI << "set terminal pdf\n";
		*POSI << "set xlabel \"position in read\"\n";
		*POSI << "set ylabel \"probability\"\n";
		*POSI << "plot [0 to 30] ";
		bool first = true;
		for ( uint64_t i = z; i < DP.size() && i < z+1; ++i )
		{
			std::ostringstream fnostr;
			fnostr << prefix << "_" << i << ".gpl";
			std::string const fn = fnostr.str();
			Vfn.push_back(fn);
			libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
			libmaus2::aio::OutputStreamInstance OSI(fn);
			DP[i].printData(OSI);

			if ( ! first )
				*POSI << ",";

			*POSI << "\"" << fn << "\" smooth bezier title \"position " << z << " in true sequence\"";
			first = false;
		}
		*POSI << "\n";
		POSI->flush();
		POSI.reset();

		std::ostringstream comstr;
		comstr << "gnuplot <" << gplfn << " >" << prefix << ".pdf";
		std::string const com = comstr.str();

		int const r = system(com.c_str());
		if ( r != EXIT_SUCCESS )
		{

		}
	}
};
#endif
