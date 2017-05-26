/*
    libmaus2
    Copyright (C) 2009-2013 German Tischler
    Copyright (C) 2011-2013 Genome Research Limited

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
#if ! defined(HETTHRESHOLD_HPP)
#define HETTHRESHOLD_HPP

#include <libmaus2/random/Poisson.hpp>
#include <libmaus2/math/binom.hpp>

struct HetThreshold
{
	static uint64_t getHetThreshold(uint64_t const d, double const e, double const thres = 0.995)
	{
		double const p = 1-e;
		std::vector<double> VPO = libmaus2::random::Poisson(d).getVector(1e-8);

		double cs = 1.0;
		uint64_t detd = 1;
		for ( uint64_t dp = 0; true; ++dp )
		{
			double s = 0.0;

			for ( uint64_t i = dp; i < VPO.size(); ++i )
			{
				double const c =
					VPO[i] *
					libmaus2::math::gpow(p,dp) *
					libmaus2::math::gpow(1-p,i-dp) *
					libmaus2::math::Binom::binomialCoefficientInteger(dp,i);

				s += c;
			}

			if ( cs - s < thres )
			{
				detd = dp;
				break;
			}
			else
			{
				cs -= s;
			}
		}

		return detd;
	}
};
#endif
