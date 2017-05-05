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
#if ! defined(NODE_HPP)
#define NODE_HPP

struct Node
{
	uint64_t v;
	uint64_t spo;
	uint64_t freq;
	uint64_t numsucc;
	uint64_t numsuccactive;
	uint64_t feaspos;
	uint64_t cfeaspos;
	uint64_t numfeaspos;
	uint64_t numcfeaspos;
	uint64_t pfostart;
	uint64_t pfosize;
	uint64_t cpfostart;
	uint64_t cpfosize;
	uint64_t plow;
	uint64_t phigh;
	uint64_t cplow;
	uint64_t cphigh;

	Node() {}
	Node(uint64_t const rv) : v(rv) {}
	Node(
		uint64_t const rv,
		uint64_t const rspo,
		uint64_t const rfreq,
		uint64_t const rpfostart,
		uint64_t const rpfosize,
		uint64_t const rcpfostart,
		uint64_t const rcpfosize,
		uint64_t const rplow,
		uint64_t const rphigh,
		uint64_t const rcplow,
		uint64_t const rcphigh
	)
	: v(rv), spo(rspo), freq(rfreq), numsucc(0), numsuccactive(0), feaspos(0), cfeaspos(0), numfeaspos(0), numcfeaspos(0),
	  pfostart(rpfostart), pfosize(rpfosize), cpfostart(rcpfostart), cpfosize(rcpfosize), plow(rplow), phigh(rphigh), cplow(rcplow), cphigh(rcphigh) {}
};

inline std::ostream & operator<<(std::ostream & out, Node const & N)
{
	return out << "Node(v=" << N.v
		<< ",spo=" << N.spo
		<< ",freq=" << N.freq
		<< ",numsucc=" << N.numsucc
		<< ",numsuccactive=" << N.numsuccactive
		<< ",feaspos=" << N.feaspos
		<< ",numfeaspos=" << N.numfeaspos
		<< ",pfostart=" << N.pfostart
		<< ",pfosize=" << N.pfosize
		<< ",plow=" << N.plow
		<< ",phigh=" << N.phigh
		<< ")";
}
#endif
