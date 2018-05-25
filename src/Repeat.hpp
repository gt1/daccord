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
#if ! defined(REPEAT_HPP)
#define REPEAT_HPP

#include <libmaus2/autoarray/AutoArray.hpp>
#include <libmaus2/util/NumberSerialisation.hpp>

struct Repeat
{
	uint64_t id;
	uint64_t abpos;
	uint64_t aepos;

	Repeat(uint64_t const rid = 0) : id(rid), abpos(0), aepos(0)
	{}

	Repeat(uint64_t const rid,
		uint64_t const rabpos,
		uint64_t const raepos
	) : id(rid), abpos(rabpos), aepos(raepos) {}

	Repeat(std::istream & in)
	:
		id(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		abpos(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		aepos(libmaus2::util::NumberSerialisation::deserialiseNumber(in))
	{}

	bool operator<(Repeat const & R) const
	{
		if ( id != R.id )
			return id < R.id;
		else if ( abpos != R.abpos )
			return abpos < R.abpos;
		else
			return aepos < R.aepos;
	}

	std::ostream & serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,id);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,abpos);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,aepos);
		return out;
	}

	std::istream & deserialise(std::istream & in)
	{
		id    = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		abpos = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		aepos = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		return in;
	}

	static libmaus2::autoarray::AutoArray<Repeat> loadArray(std::istream & in)
	{
		uint64_t const recsize = 3*sizeof(uint64_t);
		in.clear();
		in.seekg(0,std::ios::end);
		uint64_t const fs = in.tellg();
		in.clear();
		in.seekg(0,std::ios::beg);
		assert ( fs % recsize == 0 );
		uint64_t const n = fs / recsize;
		libmaus2::autoarray::AutoArray<Repeat> A(n,false);
		for ( uint64_t i = 0; i < n; ++i )
			A[i] = Repeat(in);
		return A;
	}

	static libmaus2::autoarray::AutoArray<Repeat> loadArray(std::string const & fn)
	{
		libmaus2::aio::InputStreamInstance ISI(fn);
		return loadArray(ISI);
	}
};
#endif
