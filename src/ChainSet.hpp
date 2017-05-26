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
#if ! defined(CHAINSET_HPP)
#define CHAINSET_HPP

#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/util/PrefixSums.hpp>

struct ChainSet
{
	uint64_t const f;
	std::vector<int64_t> chainprev;
	std::vector<uint64_t> childcount;
	std::vector<uint64_t> childoffsets;
	std::vector<int64_t> childlinks;
	std::vector<int64_t> childdiagdif;
	std::vector<uint64_t> chains;
	std::vector<uint64_t> chainlengths;
	uint64_t numchains;

	ChainSet(libmaus2::dazzler::align::Overlap const * RO, uint64_t const rf)
	: f(rf), chainprev(f,-1), childcount(f+1,0), childoffsets(f+1,0), childlinks(), childdiagdif(f,0), chains(), chainlengths(), numchains(0)
	{
		uint64_t b_low = 0;

		while ( b_low < f )
		{
			uint64_t b_high = b_low+1;

			while ( b_high < f && RO[b_high].bread == RO[b_low].bread && RO[b_high].isInverse() == RO[b_low].isInverse() )
				++b_high;

			// std::cerr << "bread=" << RO[b_low].bread << " cnt=" << b_high-b_low << std::endl;

			// std::sort(RO+b_low,RO+b_high,OverlapPosComparator());

			for ( uint64_t i = b_low; i < b_high; ++i )
			{
				int64_t mindif = std::numeric_limits<int64_t>::max();

				for ( uint64_t j = b_low; j < i; ++j )
				{
					if ( RO[j].path.aepos <= RO[i].path.abpos )
					{
						#if 0
						int64_t const diagend = RO[j].path.aepos - RO[j].path.bepos;
						int64_t const diagstart = RO[i].path.abpos - RO[i].path.bbpos;
						#endif
						int64_t const difa = RO[i].path.abpos - RO[j].path.aepos;
						int64_t const difb = RO[i].path.bbpos - RO[j].path.bepos;
						assert ( difa >= 0 );

						if ( difb >= 0 )
						{
							int64_t const diagdif = difa+difb; // std::abs(diagstart - diagend);

							if ( diagdif < mindif )
							{
								// std::cerr << "check " << RO[j] << " < " << RO[i] << " diagdif=" << diagdif << std::endl;
								mindif = diagdif;
								chainprev[i] = j;
								childdiagdif[i] = diagdif;
							}
						}
					}
				}
			}

			b_low = b_high;
		}

		for ( uint64_t i = 0; i < f; ++i )
			if ( chainprev[i] >= 0 )
				assert ( RO[i].bread == RO[chainprev[i]].bread );

		for ( uint64_t i = 0; i < f; ++i )
			if ( chainprev[i] >= 0 )
				childcount[chainprev[i]]++;

		childoffsets = childcount;
		libmaus2::util::PrefixSums::prefixSums(childoffsets.begin(),childoffsets.end());

		childlinks.resize(childoffsets.back(),-1);

		for ( uint64_t i = 0; i < f; ++i )
			if ( chainprev[i] >= 0 )
			{
				uint64_t const parent = chainprev[i];
				uint64_t const offset = childoffsets[ parent ]++;
				assert ( childlinks[offset] == -1 );
				childlinks [ offset ] = i;
			}

		// restore childoffsets
		childoffsets = childcount;
		libmaus2::util::PrefixSums::prefixSums(childoffsets.begin(),childoffsets.end());

		for ( uint64_t i = 0; i < f; ++i )
			for ( uint64_t j = 0; j < childcount[i]; ++j )
				assert ( RO[i].bread == RO[childlinks[childoffsets[i]+j]].bread );

		// sort childlinks by increasing order of diagdif
		for ( uint64_t i = 0; i < f; ++i )
			if ( childcount[i] )
			{
				std::vector < std::pair<int64_t,int64_t> > V(childcount[i]);

				for ( uint64_t j = 0; j < childcount[i]; ++j )
				{
					int64_t const to = childlinks[childoffsets[i]+j];
					int64_t const diagdif = childdiagdif[to];

					V.at(j) = std::pair<int64_t,int64_t>(diagdif,to);
				}

				std::sort(V.begin(),V.end());

				for ( uint64_t j = 0; j < childcount[i]; ++j )
					childlinks[childoffsets[i]+j] = V[j].second;
			}

		#if 0
		for ( uint64_t i = 0; i < f; ++i )
			for ( uint64_t j = 0; j < childcount[i]; ++j )
			{
				int64_t const from = i;
				int64_t const to = childlinks[childoffsets[i]+j];
				std::cerr << from << " " << to << " " << childdiagdif[to] << " "
					<< RO[from].getHeader() << " -> " << RO[to].getHeader() << std::endl;
			}
		#endif

		std::set<uint64_t> sunused;
		for ( uint64_t i = 0; i < f; ++i )
			sunused.insert(i);

		while ( sunused.size() )
		{
			int64_t cur = *(sunused.begin());
			uint64_t len = 0;

			while ( cur >= 0 )
			{
				len += 1;
				chains.push_back(cur);
				sunused.erase(cur);

				if ( childcount[cur] )
				{
					int64_t next = childlinks[childoffsets[cur]];
					childoffsets[cur]++;
					childcount[cur]--;
					cur = next;
				}
				else
				{
					cur = -1;
				}
			}

			chainlengths.push_back(len);
		}
		numchains = chainlengths.size();
		chainlengths.push_back(0);

		libmaus2::util::PrefixSums::prefixSums(chainlengths.begin(),chainlengths.end());
	}

	uint64_t size() const
	{
		return numchains;
	}

	uint64_t size(uint64_t const i) const
	{
		assert ( i < size() );
		return chainlengths[i+1] - chainlengths[i];
	}

	uint64_t operator()(uint64_t const i, uint64_t const j) const
	{
		return chains [ chainlengths[i] + j ];
	}
};
#endif
