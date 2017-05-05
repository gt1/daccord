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
#if ! defined(DOTPRODUCT_HPP)
#define DOTPRODUCT_HPP

#include <libmaus2/types/types.hpp>
#include <vector>
#include <ostream>
#include <cmath>

/**
 * dot product
 **/
struct DotProduct
{
	// first significant index
	uint64_t firstsign;
	// vector of coefficients
	std::vector<double> V;

	DotProduct() {}
	DotProduct(uint64_t const rfirstsign, std::vector<double> const & rV)
	: firstsign(rfirstsign), V(rV)
	{

	}

	std::ostream & printData(std::ostream & out) const
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i] )
				out << i+firstsign << "\t" << V[i] << "\n";
		return out;
	}

	uint64_t size() const
	{
		return firstsign + V.size();
	}

	double operator[](uint64_t const i) const
	{
		if ( i < firstsign )
			return 0.0;

		uint64_t j = i - firstsign;

		if ( j < V.size() )
			return V[j];
		else
			return 0.0;
	}

	void normaliseValue(uint64_t const i, double const div)
	{
		if ( i >= firstsign )
		{
			uint64_t const j = i-firstsign;

			if ( j < V.size() )
				V[j] /= div;
		}
	}

	// compute product
	double dotproduct(std::vector < double > const & O) const
	{
		double s = 0;

		// iterate over length of V
		for ( uint64_t i = 0; i < V.size(); ++i )
		{
			// compute corresponding index j on O
			uint64_t const j = firstsign + i;

			// if j is in range for O
			if ( j < O.size() )
				s += V[i] * O[j];
			// j is too large, stop
			else
				break;
		}

		return s;
	}

	// normalise the vector (make dot product between the vector and itself 1)
	void normalise()
	{
		double s = 0.0;
		for ( uint64_t i = 0; i < V.size(); ++i )
			s += V[i]*V[i];
		double const c = std::sqrt(1.0/s);
		for ( uint64_t i = 0; i < V.size(); ++i )
			V[i] *= c;
	}
};

/**
 * output operator for dot product
 **/
inline std::ostream & operator<<(std::ostream & out, DotProduct const & DP)
{
	out << "DotProduct(firstsign=" << DP.firstsign << ",";

	for ( uint64_t j = 0; j < DP.V.size(); ++j )
		out << DP.V[j] << ";";
	out << ")";

	return out;
}
#endif

