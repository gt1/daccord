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
#if ! defined(ACTIVEELEMENT_HPP)
#define ACTIVEELEMENT_HPP

#include <libmaus2/lcs/AlignmentTraceContainer.hpp>

/**
 * active sequence in multiple alignment
 **/
struct ActiveElement
{
	// base data for a sequence
	uint8_t const * ua;
	// base data for b sequence
	uint8_t const * ub;
	// current trace pointer between a and b sequence
	libmaus2::lcs::AlignmentTraceContainer::step_type const * ta;
	libmaus2::lcs::AlignmentTraceContainer::step_type const * te;
	// offset on b sequence
	uint64_t uboff;
	// error rate of alignment
	double erate;

	ActiveElement() {}
	ActiveElement(
		uint8_t const * rua,
		uint8_t const * rub,
		libmaus2::lcs::AlignmentTraceContainer::step_type const * rta,
		libmaus2::lcs::AlignmentTraceContainer::step_type const * rte,
		uint64_t const ruboff,
		double const rerate
	) : ua(rua), ub(rub), ta(rta), te(rte), uboff(ruboff), erate(rerate) {}
};
#endif
