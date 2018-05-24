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
#if ! defined(DEBRUIJNGRAPHINTERFACE_HPP)
#define DEBRUIJNGRAPHINTERFACE_HPP

struct DebruijnGraphInterface
{
	typedef DebruijnGraphInterface this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	virtual ~DebruijnGraphInterface() {}

	virtual void setup(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o) = 0;
	virtual void filterFreq(uint64_t const f, uint64_t const no) = 0;
	// static double getDefaultComputeFeasibleKmerPositionsThres() { return 1e-3; }
	virtual void computeFeasibleKmerPositions(OffsetLikely const & offsetLikely, double const thres) = 0;
	virtual void getLevelSuccessors(unsigned int const s) = 0;
	virtual void setupNodes() = 0;
	virtual void setupAddHeap(uint64_t const o) = 0;
	virtual bool addNextFromHeap(std::ostream * logstr) = 0;
	virtual bool traverse(int64_t const lmin, int64_t const lmax, std::pair< uint8_t const *, uint64_t> const * MA, uint64_t const MAo,
		uint64_t const maxfrontpath, uint64_t const maxfullpath) = 0;
	virtual std::pair<uint64_t,double> checkCandidates(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o) = 0;
	virtual std::pair<uint64_t,uint64_t> checkCandidatesU(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o) = 0;
	virtual std::pair<uint8_t const *, uint8_t const *> getCandidate(uint64_t const i) const = 0;
	virtual uint64_t getNumCandidates() const = 0;
	virtual double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I,
		uint64_t const o,
		std::pair<uint8_t const *, uint8_t const *> Ucand
	) = 0;
	virtual double getCandidateError(
		std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o,
		uint64_t const id
	) = 0;
	virtual double getCandidateWeight(uint64_t const i) const = 0;
	virtual Node const * getNodeVirtual(uint64_t const v) const = 0;
	virtual bool isEdgeActiveVirtual(uint64_t const from, uint64_t const to) const = 0;
	virtual void getActiveSuccessorsVirtual(uint64_t const v, Links & L) const = 0;
	virtual void getSuccessorsVirtual(uint64_t const v, Links & L) const = 0;
	virtual uint64_t getKmerSize() const = 0;
	virtual std::string printNode(Node const & S) const = 0;
	virtual std::ostream & print(std::ostream & out) const = 0;

	virtual uint64_t getCandidateErrorU(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o, std::pair<uint8_t const *, uint8_t const *> Ucand) = 0;
	virtual uint64_t getCandidateErrorU(std::pair< uint8_t const *, uint64_t> const * I, uint64_t const o, uint64_t const id) = 0;
	virtual uint64_t getActiveEdgeWeightVirtual(uint64_t const from, uint64_t const to) const = 0;
	virtual void printStretches(uint64_t const first, std::ostream & ostr) const = 0;
	virtual void toDot(std::ostream & out) const = 0;
};
#endif
