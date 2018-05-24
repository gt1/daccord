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
#if ! defined(DEBRUIJNGRAPHCONTAINER_HPP)
#define DEBRUIJNGRAPHCONTAINER_HPP

#include <DebruijnGraph.hpp>

struct DebruijnGraphContainer
{
	typedef DebruijnGraphContainer this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	libmaus2::autoarray::AutoArray < DebruijnGraphInterface::unique_ptr_type > ADG;

	DebruijnGraphContainer(
		double const est_cor,
		uint64_t const kmersizelow,
		uint64_t const kmersizehigh,
		std::map < uint64_t, KmerLimit::shared_ptr_type > const & MKL
	) : ADG( kmersizehigh-kmersizelow + 1)
	{
		for ( uint64_t k = kmersizelow; k <= kmersizehigh; ++k )
		{
			uint64_t const j = k - kmersizelow;

			switch ( k )
			{
				case 3:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<3>(est_cor,*(MKL.find(3)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 4:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<4>(est_cor,*(MKL.find(4)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 5:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<5>(est_cor,*(MKL.find(5)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 6:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<6>(est_cor,*(MKL.find(6)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 7:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<7>(est_cor,*(MKL.find(7)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 8:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<8>(est_cor,*(MKL.find(8)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 9:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<9>(est_cor,*(MKL.find(9)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 10:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<10>(est_cor,*(MKL.find(10)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 11:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<11>(est_cor,*(MKL.find(11)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				case 12:
				{
					DebruijnGraphInterface::unique_ptr_type tptr(new DebruijnGraph<12>(est_cor,*(MKL.find(12)->second)));
					ADG[j] = UNIQUE_PTR_MOVE(tptr);
					break;
				}
				default:
				{
					libmaus2::exception::LibMausException lme;
					lme.getStream() << "k-mer size " << k << " is not compiled in" << std::endl;
					lme.finish();
					throw lme;
					break;
				}
			}
		}
	}
};
#endif
