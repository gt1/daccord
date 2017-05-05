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
#if ! defined(LINKS_HPP)
#define LINKS_HPP

#include <libmaus2/types/types.hpp>

struct Links
{
	public:
	uint64_t A[4];
	uint64_t p;

	public:
	void reset()
	{
		p = 0;
	}

	void push(uint64_t const sym, uint64_t const freq)
	{
		if ( freq )
			A[p++] = (freq << 8) | sym;
	}

	void setSize(uint64_t const rp)
	{
		p = rp;
	}

	Links()
	{
		reset();
	}

	void sort()
	{
		if ( p <= 1 )
			return;

		std::sort(&A[0],&A[p],std::greater<uint64_t>());
	}

	uint64_t size() const
	{
		return p;
	}

	uint64_t getFreq(uint64_t const i) const
	{
		return A[i] >> 8;
	}

	uint64_t getSym(uint64_t const i) const
	{
		return A[i] & 0xFF;
	}
};
#endif
