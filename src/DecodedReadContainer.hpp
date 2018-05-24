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
#if ! defined(DECODEDREADCONTAINER_HPP)
#define DECODEDREADCONTAINER_HPP

#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>

struct ReadData
{
	typedef ReadData this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::autoarray::AutoArray<char> A;
	uint64_t l;

	ReadData()
	{
	}
};

struct ReadDataTypeInfo
{
	typedef ReadData element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		return pointer_type();
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct ReadDataAllocator
{
	typedef ReadData element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct ReadDecoder
{
	typedef ReadDecoder this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::dazzler::db::DatabaseFile const & DB;
	libmaus2::aio::InputStream::unique_ptr_type bpsfile;
	libmaus2::aio::InputStream::unique_ptr_type idxfile;

	ReadDecoder(libmaus2::dazzler::db::DatabaseFile const & rDB)
	: DB(rDB), bpsfile(DB.openBaseStream()), idxfile(DB.openIndexStream())
	{

	}

	size_t decodeRead(uint64_t const id, libmaus2::autoarray::AutoArray<char> & A) const
	{
		return DB.decodeReadAndRC(*idxfile,*bpsfile,id,A);
	}
};

struct ReadDecoderTypeInfo
{
	typedef ReadDecoder element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		return pointer_type();
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct ReadDecoderAllocator
{
	typedef ReadDecoder element_type;
	typedef element_type::shared_ptr_type pointer_type;

	libmaus2::dazzler::db::DatabaseFile const * DB;

	ReadDecoderAllocator(libmaus2::dazzler::db::DatabaseFile const * rDB)
	: DB(rDB)
	{

	}

	pointer_type operator()() const
	{
		return pointer_type(new element_type(*DB));
	}
};

struct DecodedReadContainer
{
	std::map<uint64_t,uint64_t> M;
	std::vector<ReadData::shared_ptr_type> V;

	libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & readDataFreeList;
	libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & readDecoderFreeList;

	ReadDecoder::shared_ptr_type decoder;

	DecodedReadContainer(
		libmaus2::parallel::LockedGrowingFreeList<ReadData,ReadDataAllocator,ReadDataTypeInfo> & rreadDataFreeList,
		libmaus2::parallel::LockedGrowingFreeList<ReadDecoder,ReadDecoderAllocator,ReadDecoderTypeInfo> & rreadDecoderFreeList
	) : readDataFreeList(rreadDataFreeList), readDecoderFreeList(rreadDecoderFreeList), decoder(readDecoderFreeList.get())
	{
	}

	~DecodedReadContainer()
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i] )
				readDataFreeList.put(V[i]);
		readDecoderFreeList.put(decoder);
	}

	void erase(uint64_t const id)
	{
		std::map<uint64_t,uint64_t>::iterator it = M.find(id);

		if ( it != M.end() )
		{
			ReadData::shared_ptr_type R = V[it->second];
			readDataFreeList.put(R);
			V[it->second] = ReadData::shared_ptr_type();
			M.erase(it);
		}
	}

	uint64_t ensurePresent(uint64_t const id)
	{
		std::map<uint64_t,uint64_t>::iterator it = M.find(id);

		if ( it == M.end() )
		{
			ReadData::shared_ptr_type R = readDataFreeList.get();
			R->l = decoder->decodeRead(id,R->A);

			uint64_t const vid = V.size();
			V.push_back(R);
			M[id] = vid;
		}

		it = M.find(id);
		assert ( it != M.end() );

		return it->second;
	}

	char const * getForwardRead(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->A.begin();
	}

	char const * getReverseComplementRead(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->A.begin() + R->l;
	}

	uint64_t getReadLength(uint64_t const id)
	{
		uint64_t const vid = ensurePresent(id);
		ReadData::shared_ptr_type const & R = V [ vid ];
		return R->l;
	}
};
#endif
