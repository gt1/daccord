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
#include <config.h>
#include <libmaus2/util/LineBuffer.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/DigitTable.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NPL.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/geometry/RangeSet.hpp>
#include <queue>
#include <libmaus2/parallel/NumCpus.hpp>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

struct Match
{
	int64_t off;
	int64_t l;

	Match()
	{

	}
	Match(int64_t const roff, int64_t const rl)
	: off(roff), l(rl) {}

	bool operator<(Match const & M) const
	{
		return off < M.off;
	}

	int64_t getFrom() const
	{
		return off;
	}

	int64_t getTo() const
	{
		return off + l;
	}

	libmaus2::math::IntegerInterval<int64_t> getInterval() const
	{
		return libmaus2::math::IntegerInterval<int64_t>(getFrom(),getTo()-1);
	}
};

struct Kmer
{
	uint64_t v;
	uint64_t p;

	Kmer() {}
	Kmer(uint64_t const rv, uint64_t const rp) : v(rv), p(rp) {}

	bool operator<(Kmer const & K) const
	{
		if ( v != K.v )
			return v < K.v;
		else
			return p < K.p;
	}
};

static uint64_t getKmers(std::string const & s, uint64_t const k, libmaus2::autoarray::AutoArray<Kmer> & K)
{
	uint64_t const n = s.size();
	uint64_t const numk = (n >= k) ? (n-k+1) : 0;
	uint64_t const kmask = libmaus2::math::lowbits(2*(k-1));
	assert ( k );

	uint64_t v = 0;
	for ( uint64_t i = 0; i < k-1; ++i )
	{
		v <<= 2;
		v |= libmaus2::fastx::mapChar(s[i]);
	}

	uint64_t o = 0;
	for ( uint64_t i = 0; i < numk ; ++i )
	{
		v &= kmask;
		v <<= 2;
		v |= libmaus2::fastx::mapChar(s[i+k-1]);
		K.push(o,Kmer(v,i));
	}

	std::sort(K.begin(),K.begin()+o);

	return o;
}

struct HeapTodo
{
	Kmer const * ka;
	Kmer const * kb_c;
	Kmer const * kb_e;

	HeapTodo()
	{

	}

	HeapTodo(
		Kmer const * rka,
		Kmer const * rkb_c,
		Kmer const * rkb_e)
	: ka(rka), kb_c(rkb_c), kb_e(rkb_e)
	{

	}

	bool hasNext() const
	{
		return kb_c+1 != kb_e;
	}

	HeapTodo getNext()
	{
		HeapTodo T = *this;
		++T.kb_c;
		return T;
	}

	int64_t getDiag() const
	{
		return static_cast<int64_t>(ka->p) - static_cast<int64_t>(kb_c->p);
	}

	bool operator<(HeapTodo const & O) const
	{
		return getDiag() < O.getDiag();
	}

	bool operator>(HeapTodo const & O) const
	{
		return getDiag() > O.getDiag();
	}
};

struct KmerRange
{
	uint64_t v;
	uint64_t aindex;
	uint64_t ap;
	uint64_t bfrom;
	uint64_t bto;

	KmerRange() {}
	KmerRange(
		uint64_t const rv,
		uint64_t const raindex,
		uint64_t const rap,
		uint64_t const rbfrom,
		uint64_t const rbto
	) : v(rv), aindex(raindex), ap(rap), bfrom(rbfrom), bto(rbto)
	{

	}
};

std::ostream & operator<<(std::ostream & out, KmerRange const & K)
{
	out << "KmerRange(v=" << K.v << ",aindex=" << K.aindex << ",ap=" << K.ap << ",bfrom=" << K.bfrom << ",bto=" << K.bto << ")";
	return out;
}

struct KmerRangePosComparator
{
	bool operator()(KmerRange const & A, KmerRange const & B) const
	{
		return A.ap < B.ap;
	}
};

struct SeedAndExtend
{
	typedef SeedAndExtend this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	int64_t minlen;
	int64_t w_err;
	int64_t w_back;
	unsigned int k;
	int64_t tspace;

	libmaus2::autoarray::AutoArray<Kmer> KA;
	libmaus2::autoarray::AutoArray<Kmer> KB;
	libmaus2::autoarray::AutoArray<KmerRange> Arange;
	libmaus2::lcs::NNP nnp;
	libmaus2::lcs::NNPTraceContainer nnptracecontainer;
	libmaus2::autoarray::AutoArray<uint64_t> Ascore;
	libmaus2::autoarray::AutoArray<int64_t> Alast;
	libmaus2::lcs::AlignmentTraceContainer ATC;

	std::string ra;
	uint64_t oa;
	int64_t aid;

	SeedAndExtend(int64_t const rminlen, int64_t const rw_err, int64_t const rw_back, unsigned int const rk, int64_t const rtspace)
	: minlen(rminlen), w_err(rw_err), w_back(rw_back), k(rk), tspace(rtspace), nnp(w_err,w_back), oa(0), aid(-1)
	{}

	void setupA(std::string const & sa, uint64_t const raid)
	{
		ra = sa;
		oa = getKmers(ra,k,KA);
		aid = raid;
	}

	bool seedAndExtend(std::string const & sb, int64_t const bid, libmaus2::dazzler::align::Overlap & ROVL)
	{
		int64_t const minscore = 3*k;

		std::string rb = sb;

		// compute k-mers in read a and b
		uint64_t const ob = getKmers(rb,k,KB);

		uint64_t ia = 0;
		uint64_t ib = 0;

		uint64_t oArange = 0;

		unsigned int const bucketshift = 6;
		uint64_t const bucketsize = 1ull << bucketshift;

		assert ( ra.size() );
		assert ( rb.size() );

		int64_t const mindiag = - (rb.size() - 1);
		int64_t const maxdiag = ra.size() - 1;
		uint64_t const numdiag = static_cast<uint64_t>(maxdiag - mindiag + 1);
		uint64_t const abuckets = (numdiag + bucketsize - 1)/bucketsize;
		Ascore.ensureSize(abuckets);
		std::fill(Ascore.begin(),Ascore.begin()+abuckets,0ull);
		Alast.ensureSize(abuckets);
		std::fill(Alast.begin(),Alast.begin()+abuckets,-static_cast<int64_t>(k));

		// look for matching k-mers
		while ( ia < oa && ib < ob )
		{
			if ( KA[ia].v < KB[ib].v )
				++ia;
			else if ( KB[ib].v < KA[ia].v )
				++ib;
			else
			{
				assert ( KA[ia].v == KB[ib].v );

				uint64_t ja = ia+1;
				while ( ja < oa && KA[ja].v == KA[ia].v )
					++ja;
				uint64_t jb = ib+1;
				while ( jb < ob && KB[jb].v == KB[ib].v )
					++jb;

				for ( uint64_t z = ia; z < ja; ++z )
					Arange.push(oArange,KmerRange(KA[z].v,z,KA[z].p,ib,jb));

				#if 0
				// reverse region in B so positions on B are in decreasing order
				std::reverse(KB.begin()+ib,KB.begin()+jb);

				for ( uint64_t z = ia; z < ja; ++z )
					Q.push(HeapTodo(KA.begin()+z,KB.begin()+ib,KB.begin()+jb));
				#endif

				ia = ja;
				ib = jb;
			}
		}

		std::sort(Arange.begin(),Arange.begin() + oArange,KmerRangePosComparator());

		for ( uint64_t i = 0; i < oArange; ++i )
		{
			KmerRange const & K = Arange[i];
			int64_t const ap = K.ap;

			for ( uint64_t z = K.bfrom; z < K.bto; ++z )
			{
				int64_t const bp = KB[z].p;
				int64_t const d = ap-bp;
				assert ( d - mindiag >= 0 );
				uint64_t const bucket = (d - mindiag) >> bucketshift;
				// std::cerr << Arange[i] << std::endl;

				int64_t const preva = Alast [ bucket ];
				assert ( preva <= ap );
				int64_t const score = std::min(ap-preva,static_cast<int64_t>(k));
				assert ( score >= 0 );

				// std::cerr << ap << " " << bp << " " << d << " " << bucket << " " << score << std::endl;

				Alast[bucket] = ap;
				Ascore[bucket] += score;
			}
		}

		std::fill(Alast.begin(),Alast.begin()+abuckets,0);

		int64_t maxlen = 0;
		for ( uint64_t i = 0; i < oArange; ++i )
		{
			KmerRange const & K = Arange[i];
			int64_t const ap = K.ap;

			for ( uint64_t z = K.bfrom; z < K.bto; ++z )
			{
				int64_t const bp = KB[z].p;
				int64_t const d = ap-bp;
				assert ( d - mindiag >= 0 );
				uint64_t const bucket = (d - mindiag) >> bucketshift;

				int64_t score = Ascore[bucket];
				if ( bucket > 0 )
					score += Ascore[bucket-1];
				if ( bucket+1 < abuckets )
					score += Ascore[bucket+1];

				if ( score >= minscore && ap >= Alast[bucket] )
				{
					libmaus2::lcs::NNPAlignResult res;

					if ( aid == bid )
						res = nnp.align(
							ra.begin(),ra.end(),ap,
							ra.begin(),ra.end(),bp,
							nnptracecontainer,
							aid==bid,
							nnp.getDefaultMinDiag(),
							nnp.getDefaultMaxDiag(),
							true /* runsuffixpositive */
						);
					else
						res = nnp.align(
							ra.begin(),ra.end(),ap,
							rb.begin(),rb.end(),bp,
							nnptracecontainer,
							aid==bid,
							nnp.getDefaultMinDiag(),
							nnp.getDefaultMaxDiag(),
							true /* runsuffixpositive */
						);

					int64_t const span = res.aepos - res.abpos;

					if ( span >= minlen && span > maxlen )
					{
						// std::cerr << ap << " " << bp << " " << bucket << " " << res << " " << span << std::endl;

						// compute dense trace
						nnptracecontainer.computeTrace(ATC);

						libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
							0 /* flags */,
							aid,
							bid,
							res.abpos,
							res.aepos,
							res.bbpos,
							res.bepos,
							tspace,
							ATC
						);

						ROVL = OVL;
						maxlen = span;
					}

					Alast[bucket] = res.aepos;
				}
			}
		}

		return maxlen > 0;
	}
};

struct Read
{
	uint64_t id;
	std::string data;

	Read() {}
	Read(uint64_t const rid, std::string const & rdata) : id(rid), data(rdata) {}

	bool operator<(Read const & O) const
	{
		return id < O.id;
	}
};

struct ReadPackage
{
	Read read;
	std::vector<Read> Vdata;
	bool gotread;

	ReadPackage()
	: gotread(false)
	{

	}

	void sort()
	{
		std::sort(Vdata.begin(),Vdata.end());
	}

	void reset()
	{
		gotread = false;
		Vdata.resize(0);
	}
};

struct Handler
{
	uint64_t numthreads;
	std::string const tmpfilenamebase;
	libmaus2::autoarray::AutoArray<SeedAndExtend::unique_ptr_type> ASAE;
	libmaus2::dazzler::db::DatabaseFile const & DB;
	std::string outname;
	libmaus2::dazzler::align::AlignmentWriterArray AWA;

	Handler(int64_t const minlen, int64_t const w_err, int64_t const w_back, int64_t const k, int64_t const tspace, libmaus2::dazzler::db::DatabaseFile const & rDB, std::string const & routname, uint64_t const rnumthreads, std::string const & rtmpfilenamebase)
	: numthreads(rnumthreads),  tmpfilenamebase(rtmpfilenamebase), ASAE(numthreads), DB(rDB), outname(routname), AWA(tmpfilenamebase+"_AWA",numthreads,tspace)
	{
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			SeedAndExtend::unique_ptr_type tptr(new SeedAndExtend(minlen,w_err,w_back,k,tspace));
			ASAE[i] = UNIQUE_PTR_MOVE(tptr);
		}
	}

	void merge()
	{
		AWA.merge(outname,tmpfilenamebase + "_mergetmp");
	}

	void handle(std::vector<ReadPackage> & VRP)
	{
		#if defined(_OPENMP)
		#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
		#endif
		for ( uint64_t j = 0; j < VRP.size(); ++j )
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			SeedAndExtend & SAE = *(ASAE[tid]);
			libmaus2::dazzler::align::AlignmentWriter & AW = AWA[tid];

			ReadPackage & P = VRP[j];

			assert ( P.gotread );
			assert ( P.read.id );

			P.sort();

			uint64_t ida1 = P.read.id-1;
			::std::size_t aoff = std::string::npos;
			std::string sa = DB[ida1];
			bool const aok = (aoff=sa.find(P.read.data)) != std::string::npos;
			assert ( aok );

			SAE.setupA(P.read.data,ida1);

			for ( uint64_t z = 0; z < P.Vdata.size(); ++z )
			{
				assert ( P.Vdata[z].id );
				uint64_t idb1 = P.Vdata[z].id-1;

				::std::size_t boff = std::string::npos;
				bool rok = false;
				bool inv = false;
				std::string sb = DB[idb1];

				rok = rok || ((boff = sb.find(P.Vdata[z].data)) != std::string::npos);

				if ( ! rok )
				{
					sb = libmaus2::fastx::reverseComplementUnmapped(sb);
					inv = true;
				}

				rok = rok || ((boff=sb.find(P.Vdata[z].data)) != std::string::npos);

				assert ( rok );

				libmaus2::dazzler::align::Overlap OVL;
				bool const ok = SAE.seedAndExtend(P.Vdata[z].data,idb1,OVL);

				OVL.path.abpos += aoff;
				OVL.path.aepos += aoff;
				OVL.path.bbpos += boff;
				OVL.path.bepos += boff;
				if ( inv )
					OVL.flags |= libmaus2::dazzler::align::Overlap::getInverseFlag();

				if ( ok )
				{
					// std::cerr << OVL << std::endl;
					// std::cerr << read.id << " " << Vdata[z].id << std::endl;
					AW.put(OVL);
				}
				else
				{
					std::cerr << "[W] no significant overlap found for " << P.read.id << " -> " << P.Vdata[z].id << std::endl;
				}
			}

			// std::cerr << P.read.id << std::endl;

			P.reset();
		}

		VRP.resize(0);
	}
};

int canutolas(libmaus2::util::ArgParser const & arg)
{
	libmaus2::timing::RealTimeClock rtc;
	rtc.start();

	libmaus2::util::LineBuffer LB(std::cin);

	std::string const outname = arg[1];
	std::string const dbname = arg[0];
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	std::cerr << "[V] copying " << dbname << " to memory...";
	// libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbname,"mem:db1prefix"));
	libmaus2::dazzler::db::DatabaseFile::DBArrayFileSet::unique_ptr_type dbptr(
		libmaus2::dazzler::db::DatabaseFile::copyToArrays(dbname)
	);
	std::cerr << "done." << std::endl;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->getDBURL()));
	// PDB->computeTrimVector();

	char const * a;
	char const * e;

	char const * tokenread = "read";
	::std::size_t l_tokenread = strlen(tokenread);
	char const * tokendata = "data";
	::std::size_t l_tokendata = strlen(tokendata);
	libmaus2::util::DigitTable const DT;
	char const * tokenrecend = "+ +";
	::std::size_t l_tokenrecend = strlen(tokenrecend);

	ReadPackage P;

	int const w_err = arg.uniqueArgPresent("w") ? arg.getUnsignedNumericArg<uint64_t>("w") : libmaus2::lcs::NNP::getDefaultMaxWindowError();
	int const w_back = arg.uniqueArgPresent("b") ? arg.getUnsignedNumericArg<uint64_t>("b") : libmaus2::lcs::NNP::getDefaultMaxBack();
	int64_t const minlen = 100;
	unsigned int const k = 14;
	int64_t const tspace = 100;
	uint64_t const batchsize = 1024;

	Handler handler(minlen, w_err, w_back, k, tspace, *PDB, outname, numthreads, tmpfilebase);
	uint64_t processed = 0;

	std::vector<ReadPackage> VRP;
	while ( LB.getline(&a,&e) )
	{
		if ( strncmp(a,tokenread,l_tokenread) == 0 )
		{
			assert ( ! P.gotread );
			assert ( P.Vdata.size() == 0 );

			uint64_t id = 0;
			char const * p = a + l_tokenread;
			while ( p != e && DT[*p] )
			{
				id *= 10;
				id += (*(p++)) - '0';
			}

			assert ( p != e && *p == ' ');
			p++;

			P.read = Read(id,std::string(p,e));
			P.gotread = true;
		}
		else if ( strncmp(a,tokendata,l_tokendata) == 0 )
		{
			assert ( P.gotread );

			uint64_t id = 0;
			char const * p = a + l_tokendata;
			while ( p != e && DT[*p] )
			{
				id *= 10;
				id += (*(p++)) - '0';
			}

			assert ( p != e && *p == ' ');
			p++;

			P.Vdata.push_back(Read(id,std::string(p,e)));
		}
		else if ( strncmp(a,tokenrecend,l_tokenrecend) == 0 )
		{
			VRP.push_back(P);
			P.reset();

			if ( VRP.size() >= batchsize )
			{
				processed += VRP.size();
				handler.handle(VRP);

				std::cerr << "[V] processed " << processed << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " " << processed / rtc.getElapsedSeconds() << std::endl;
			}
		}
	}

	if ( VRP.size() )
	{
		processed += VRP.size();
		handler.handle(VRP);
		std::cerr << "[V] processed " << processed << " " << rtc.formatTime(rtc.getElapsedSeconds()) << " " << processed / rtc.getElapsedSeconds() << std::endl;
	}

	handler.merge();

	std::cerr << "[V] total time " << rtc.formatTime(rtc.getElapsedSeconds()) << std::endl;

	return EXIT_SUCCESS;
}

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

static std::string helpMessage(libmaus2::util::ArgParser const & arg)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));

	uint64_t maxlhs = 0;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		assert ( ita->first.size() );

		if ( ita->first.size() == 1 )
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+1));
		else
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+2));
	}

	std::ostringstream messtr;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		std::string const key = ita->first;

		messtr << "\t";
		messtr << std::setw(maxlhs) << std::setfill(' ');
		if ( key.size() == 1 )
			messtr << (std::string("-")+key);
		else
			messtr << (std::string("--")+key);

		messtr << std::setw(0);

		messtr << ": ";

		messtr << ita->second;
		messtr << "\n";
	}

	return messtr.str();
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( (arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help")) || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] in.db out.las\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return canutolas(arg);
		}


	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
