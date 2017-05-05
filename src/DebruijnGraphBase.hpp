/*
    libmaus2
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
#if ! defined(DEBRUIJNGRAPHBASE_HPP)
#define DEBRUIJNGRAPHBASE_HPP

#include <libmaus2/lcs/Aligner.hpp>
#include <libmaus2/lcs/AlignerFactory.hpp>

struct DebruijnGraphBase
{
	static libmaus2::lcs::Aligner::unique_ptr_type getAligner()
	{
		std::set<libmaus2::lcs::AlignerFactory::aligner_type> const S = libmaus2::lcs::AlignerFactory::getSupportedAligners();

		libmaus2::lcs::AlignerFactory::aligner_type T[] = {
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_y256_8,
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_x128_8,
			libmaus2::lcs::AlignerFactory::libmaus2_lcs_AlignerFactory_NP
		};

		for ( uint64_t i = 0; i < sizeof(T)/sizeof(T[0]); ++i )
			if ( S.find(T[i]) != S.end() )
			{
				libmaus2::lcs::Aligner::unique_ptr_type tptr(libmaus2::lcs::AlignerFactory::construct(T[i]));
				return UNIQUE_PTR_MOVE(tptr);
			}

		libmaus2::exception::LibMausException lme;
		lme.getStream() << "DebruijnGraph::getAligner: no suitable class found" << std::endl;
		lme.finish();
		throw lme;
	}
};
#endif
