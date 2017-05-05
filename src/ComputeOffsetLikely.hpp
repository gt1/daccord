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
#if ! defined(COMPUTEOFFSETLIKELY_HPP)
#define COMPUTEOFFSETLIKELY_HPP

#include <OffsetLikely.hpp>
#include <libmaus2/math/Convolution.hpp>
#include <libmaus2/math/GmpFloat.hpp>
#include <libmaus2/math/binom.hpp>

OffsetLikely computeOffsetLikely(
	// maximum length considered
	uint64_t const maxl,
	// I prob
	double const p_i,
	// D prob
	double const p_d
)
{
	// std::vector < DotProduct > VD;
	OffsetLikely VD;

	// compute probability bins for number of insertions at/before a single position
	// double const p_i = 0.1;
	double const q_i = 1.0 - p_i;
	double f_i = q_i;

	std::vector < double > P_I;
	while ( f_i >= 1e-7 )
	{
		P_I.push_back(f_i);
		f_i *= p_i;
	}

	// convolution accumulation vector for insertions
	std::vector < double > C_I(1,1.0);

	// double const p_d = 0.02;

	// iterate over length
	for ( uint64_t l = 0; l <= maxl; ++l )
	{
		// std::cerr << std::string(40,'*') << l << std::string(40,'*') << std::endl;
		//
		// add another round of insertions
		#if 1
		C_I = libmaus2::math::Convolution::convolutionFFTRef(C_I,P_I);
		#else
		C_I = libmaus2::math::Convolution::convolutionFFT(C_I,P_I);
		#endif

		// likelihood of deletion count
		std::vector < libmaus2::math::GmpFloat > V_D = libmaus2::math::Binom::binomVector(p_d, l, 512);
		assert ( V_D.size() );

		// insertion vector with V_D.size()-1 leading zeroes
		std::vector < double > V_I (V_D.size()-1+C_I.size());
		std::copy(
			C_I.begin(),C_I.end(),
			V_I.begin() + (V_D.size()-1)
		);

		// reverse deletion vector
		std::reverse(V_D.begin(),V_D.end());

		#if 0
		for ( uint64_t j = 0; j < V_D.size(); ++j )
			if ( static_cast<double>(V_D[j]) >= 1e-5 )
				std::cerr << "D j=" << j << " V_D=" << V_D[j] << std::endl;
		for ( uint64_t j = 0; j < V_I.size(); ++j )
			if ( V_I[j] >= 1e-5 )
				std::cerr << "I j=" << j << " V_I=" << V_I[j] << std::endl;
		#endif

		// compute convolution
		std::vector < double > F_I = libmaus2::math::Convolution::convolutionFFTRef(V_D,V_I);
		// first significant value found
		bool signfound = false;
		// index of first significant value
		int64_t firstsign = std::numeric_limits<int64_t>::min();
		// value vector
		std::vector < double > VP;

		for ( uint64_t j = 0; j < F_I.size(); ++j )
			if ( F_I[j] >= 1e-5 )
			{
				if ( ! signfound )
				{
					signfound = true;
					firstsign = static_cast<int64_t>(j)-static_cast<int64_t>(l);
				}

				uint64_t const offset =
					static_cast<int64_t>(j)-static_cast<int64_t>(l) - firstsign;

				while ( ! (offset < VP.size()) )
					VP.push_back(0);
				VP[offset] = F_I[j];

				// std::cerr << "F j=" << static_cast<int64_t>(j)-static_cast<int64_t>(2*l) << " abs=" << j-l << " off=" << offset << " F_I=" << F_I[j] << std::endl;
			}
		assert ( firstsign >= 0 );

		VD.push_back(DotProduct(firstsign,VP));
		// std::cerr << "first sign " << firstsign << std::endl;

		// add another round of insertions
		//C_I = libmaus2::math::Convolution::convolutionFFT(C_I,P_I);
	}

	#if 0
	for ( uint64_t i = 0; i < VD.size(); ++i )
		VD[i].normalise();
	#endif

	VD.setup();

	return VD;
}
#endif
