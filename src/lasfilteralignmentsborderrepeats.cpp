/*
    daccord
    Copyright (C) 2016-2017 German Tischler

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
#include <RepeatIdComparator.hpp>
#include <config.h>

#include <libmaus2/dazzler/align/OverlapProperCheck.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/dazzler/db/InqualContainer.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>
#include <libmaus2/aio/SerialisedPeeker.hpp>
#include <libmaus2/dazzler/align/OverlapInfoIndexer.hpp>
#include <libmaus2/dazzler/align/LasIntervals.hpp>

struct IsProper
{
	typedef IsProper this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPVA;
	libmaus2::autoarray::AutoArray<libmaus2::dazzler::align::TracePoint> TPVB;
	libmaus2::dazzler::db::DatabaseFile const & DB;
	int64_t tspace;
	libmaus2::lcs::NNP nnp;
	libmaus2::lcs::NNPTraceContainer nnptracecontainer;
	libmaus2::lcs::AlignmentTraceContainer ATC;
	libmaus2::dazzler::align::OverlapProperCheck const & OPC;

	IsProper(
		libmaus2::dazzler::db::DatabaseFile const & rDB,
		libmaus2::dazzler::align::OverlapProperCheck const & rOPC,
		int64_t const rtspace
	) : DB(rDB), tspace(rtspace), nnp(30,125), OPC(rOPC) {}

	/*
	 * check whether there is a proper overlap by the B reads given in overlaps A and B
	 */
	bool isProper(
		libmaus2::dazzler::align::Overlap const & A,
		libmaus2::dazzler::align::Overlap const & B,
		libmaus2::lcs::NNPAlignResult * pres = 0
	)
	{
		// get trace points
		uint64_t oa = A.getTracePoints(tspace,A.bread,TPVA,0);
		uint64_t ob = B.getTracePoints(tspace,B.bread,TPVB,0);

		// sync up on A
		uint64_t ia = 0, ib = 0;
		while (
			ia < oa && ib < ob && (TPVA[ia].apos != TPVB[ib].apos || TPVA[ia].apos % tspace != 0)
		)
		{
			if ( TPVA[ia].apos < TPVB[ib].apos )
				++ia;
			else
				++ib;
		}

		// no sync point, return false
		if ( ia == oa || ib == ob )
			return false;

		assert ( ia != oa );
		assert ( ib != ob );

		// get position on C and D reads
		uint64_t const cpos = TPVA[ia].bpos;
		uint64_t const dpos = TPVB[ib].bpos;

		// get read data
		bool ainv = A.isInverse();
		std::string cread = DB.decodeRead( A.bread, ainv );
		bool binv = B.isInverse();
		std::string dread = DB.decodeRead( B.bread, binv );

		// compute alignment
		libmaus2::lcs::NNPAlignResult res = nnp.align(
			cread.begin(),cread.end(),cpos,
			dread.begin(),dread.end(),dpos,
			nnptracecontainer,
			false /* self */
		);

		// compute dense trace
		nnptracecontainer.computeTrace(ATC);

		// apply reverse complement if necessary so C read is not revcomp
		if ( ainv )
		{
			ATC.reverse();
			std::swap(res.abpos,res.aepos);
			res.abpos = cread.size() - res.abpos;
			res.aepos = cread.size() - res.aepos;
			std::swap(res.bbpos,res.bepos);
			res.bbpos = dread.size() - res.bbpos;
			res.bepos = dread.size() - res.bepos;
			ainv ^= 1;
			binv ^= 1;

			cread = libmaus2::fastx::reverseComplementUnmapped(cread);
			dread = libmaus2::fastx::reverseComplementUnmapped(dread);
		}

		// checks
		assert ( ! ainv );
		assert ( cread == DB.decodeRead(A.bread,ainv) );
		assert ( dread == DB.decodeRead(B.bread,binv) );
		assert (
			libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
				ATC.ta,
				ATC.te,
				cread.begin()+res.abpos,dread.begin()+res.bbpos
			)
		);

		#if 0
		if ( VOVL[a].isInverse() )
			libmaus2::lcs::AlignmentPrint::printAlignmentLines(
				std::cerr,
				cread.begin()+res.abpos,
				res.aepos-res.abpos,
				dread.begin()+res.bbpos,
				res.bepos-res.bbpos,
				80,
				ATC.ta,
				ATC.te
			);
		#endif

		// compute dazzler style overlap data structure
		libmaus2::dazzler::align::Overlap const NOVL = libmaus2::dazzler::align::Overlap::computeOverlap(
			binv ? libmaus2::dazzler::align::Overlap::getInverseFlag() : 0,
			A.bread,
			B.bread,
			res.abpos,
			res.aepos,
			res.bbpos,
			res.bepos,
			tspace,
			ATC
		);

		//bool const proper = OPC(NOVL,tspace,1,&std::cerr);

		// check whether alignment is proper
		bool const proper = OPC(NOVL,tspace).proper;

		if ( pres )
			*pres = res;

		return proper;
	}
};

void handle(
	std::ostream & err,
	libmaus2::autoarray::AutoArray<Repeat> const & R,
	libmaus2::dazzler::db::DatabaseFile const & DB,
	libmaus2::dazzler::align::OverlapProperCheck const & OPC,
	IsProper & isproper,
	std::vector < libmaus2::math::IntegerInterval<int64_t> > & RIV,
	std::vector < libmaus2::dazzler::align::Overlap > & VOVL,
	int64_t const prevaread,
	int64_t const lthres,
	libmaus2::dazzler::align::AlignmentWriter & AW,
	std::ostream & symkill,
	int64_t const tspace,
	double const termval
)
{
	RepeatIdComparator comp;

	// merge repeats on prevaread
	std::pair<Repeat const *, Repeat const *> const P = ::std::equal_range(R.begin(),R.end(),Repeat(prevaread),comp);
	for ( Repeat const * p = P.first; p != P.second; ++p )
		RIV.push_back(libmaus2::math::IntegerInterval<int64_t>(p->abpos,p->aepos-1));

	// merge intervals
	RIV = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(RIV);

	// print repeat intervals
	for ( uint64_t i = 0; i < RIV.size(); ++i )
		err << "[VR] " << prevaread << " " << RIV[i] << " " << DB[prevaread].size() << std::endl;

	// boolean vector of alignments we keep
	std::vector<bool> Vkeep(VOVL.size(),true);

	// if there are any repeats on A read
	if ( RIV.size() )
	{
		// vector marking alignments which are totally contained in repeat regions
		std::vector<bool> replonly(VOVL.size(),false);

		// number of alignments which have repeat free segments on A
		uint64_t numreplfree = 0;
		// number of alignments covered by repeats on A
		uint64_t numreplonly = 0;

		// iterate over alignments and mark such which are covered/contained by repeats on A
		for ( uint64_t i = 0; i < VOVL.size(); ++i )
		{
			// get A interval
			libmaus2::math::IntegerInterval<int64_t> IA(VOVL[i].path.abpos,VOVL[i].path.aepos-1);
			// intersect with repeat intervals
			std::vector < libmaus2::math::IntegerInterval<int64_t> > QIV;
			for ( uint64_t j = 0; j < RIV.size(); ++j )
				QIV.push_back(IA.intersection(RIV[j]));
			// merge intervals
			QIV = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(QIV);
			// sum up diameters
			int64_t diam = 0;
			for ( uint64_t j = 0; j < QIV.size(); ++j )
				diam += QIV[j].diameter();
			// sanity check
			assert ( diam <= IA.diameter() );

			// if most is covered by repeat areas
			if ( diam+lthres >= static_cast<int64_t>(IA.diameter()) )
			{
				// err << "repoverlapdiam=" << diam << " Q=" << IA.diameter() << " " << VOVL[i] << std::endl;
				replonly[i] = true;
				numreplonly++;
			}
			else
			{
				replonly[i] = false;
				numreplfree++;
			}
		}

		std::vector<bool> replcheck = replonly;
		std::vector<uint64_t> Vnonrepl;
		uint64_t const foundextensionthres = 1;

		for ( uint64_t i = 0; i < VOVL.size(); ++i )
			if ( ! replonly[i] )
				Vnonrepl.push_back(i);

		// iterate over alignments
		for ( uint64_t i = 0; i < VOVL.size(); ++i )
			// if alignment is covered by repeats on A
			if ( replonly[i] )
			{
				// get interval
				libmaus2::math::IntegerInterval<int64_t> IC(VOVL[i].path.abpos,VOVL[i].path.aepos-1);

				uint64_t foundextension = 0;

				// check for an extension into a non repeat area
				for ( uint64_t jj = 0; foundextension < foundextensionthres && jj < Vnonrepl.size(); ++jj )
				{
					// non repeat alignment
					uint64_t const j = Vnonrepl[jj];

					assert ( ! replonly[j] );

					// get interval
					libmaus2::math::IntegerInterval<int64_t> ID(VOVL[j].path.abpos,VOVL[j].path.aepos-1);
					// compute intersection
					libmaus2::math::IntegerInterval<int64_t> IE = IC.intersection(ID);

					// if alignment covers most of the repeat covered alignment
					if ( !IE.isEmpty() && IE.diameter() >= 0.95 * IC.diameter() )
					{
						// read length i
						uint64_t const isize = (*(OPC.RL))[VOVL[i].bread];
						// remaining on the right of i
						uint64_t const iright = isize - VOVL[i].path.bepos;
						// read length j
						uint64_t const jsize = (*(OPC.RL))[VOVL[j].bread];
						// remaining on the right of j
						uint64_t const jright = jsize - VOVL[j].path.bepos;
						// remaining on the right of both
						uint64_t const eright = std::min(iright,jright);

						// remaining on the left of i
						uint64_t const ileft = VOVL[i].path.bbpos;
						// remaining on the left of j
						uint64_t const jleft = VOVL[j].path.bbpos;
						// remaining on both
						uint64_t const eleft = std::min(ileft,jleft);

						// if both have >= lthres remaining on the left
						if ( static_cast<int64_t>(eleft) >= lthres )
						{
							// get error for i on the left
							double const errleft_i = OPC.getMaxErrorForRange(
								VOVL[i].bread,
								VOVL[i].path.bbpos-lthres,
								VOVL[i].path.bbpos,
								tspace,
								VOVL[i].isInverse()
							);
							//get error for j on the left
							double const errleft_j = OPC.getMaxErrorForRange(
								VOVL[j].bread,
								VOVL[j].path.bbpos-lthres,
								VOVL[j].path.bbpos,
								tspace,
								VOVL[j].isInverse()
							);
							// get combined error
							double const errleft = std::max(errleft_i,errleft_j);

							// if both are below termval, then reads should align in this region
							if ( errleft <= termval)
							{
								// err << "errleft=" << errleft << std::endl;
								foundextension++;
							}
						}
						// if both have >= lthres remaining on the right
						if ( static_cast<int64_t>(eright) >= lthres )
						{
							// get error for i on the right
							double const errright_i = OPC.getMaxErrorForRange(
								VOVL[i].bread,
								VOVL[i].path.bepos,
								VOVL[i].path.bepos+lthres,
								tspace,
								VOVL[i].isInverse()
							);
							// get error for j on the right
							double const errright_j = OPC.getMaxErrorForRange(
								VOVL[j].bread,
								VOVL[j].path.bepos,
								VOVL[j].path.bepos+lthres,
								tspace,
								VOVL[j].isInverse()
							);
							// get combined error
							double const errright = std::max(errright_i,errright_j);

							// if both are below termval, then reads should align in this region
							if ( errright <= termval )
							{
								// err << "errright=" << errright << std::endl;
								foundextension++;
							}
						}
					}
				}

				// err << "foundextension=" << foundextension << std::endl;

				// if we found insufficient extensions, then remove read from check list
				if ( foundextension < foundextensionthres ) // could have an evidence threshold here
					replcheck[i] = false;
			}

		err << "[V] " << prevaread << " numreplfree=" << numreplfree << " numreplonly=" << numreplonly << std::endl;

		// iterate over reads we should check
		for ( uint64_t i = 0; i < VOVL.size(); ++i )
			if ( replcheck[i] )
			{
				// found extension?
				bool found = false;

				// see if we find a proper extension
				for ( uint64_t jj = 0; !found && jj < Vnonrepl.size(); ++jj )
				{
					uint64_t const j = Vnonrepl[jj];

					assert ( ! replonly[j] );

					libmaus2::lcs::NNPAlignResult res;
					if ( isproper.isProper(VOVL[i],VOVL[j],&res) )
					{
						libmaus2::math::IntegerInterval<int64_t> IA(VOVL[i].path.bbpos,VOVL[i].path.bepos-1);
						libmaus2::math::IntegerInterval<int64_t> IQ(res.bbpos,res.bepos-1);

						if ( IQ.diameter() >= IA.diameter() + lthres )
						{
							found = true;

							#if 0
							std::cerr << "x found\n";
							std::cerr << VOVL[j] << "\n";
							std::cerr << VOVL[i] << "\n";
							std::cerr << res << std::endl;
							std::cerr << IA << std::endl;
							std::cerr << IQ << std::endl;
							#endif
						}
					}
				}

				// if no extension is found, then do not keep alignment
				if ( !found )
					Vkeep[i] = false;
			}

		err << "[V] " << prevaread << " replcheck done" << std::endl;
	}

	uint64_t written = 0;
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
		if ( Vkeep[i] )
		{
			AW.put(VOVL[i]);
			++written;
		}
		else
		{
			libmaus2::dazzler::align::OverlapInfo info = VOVL[i].getHeader().getInfo().swapped();

			if ( VOVL[i].isInverse() )
			{
				uint64_t const alen = (*(OPC.RL))[info.aread >> 1];
				uint64_t const blen = (*(OPC.RL))[info.bread >> 1];
				info = info.inverse(alen,blen);
			}

			assert ( (info.aread & 1) == 0 );

			info.serialise(symkill);
		}

	err << "[V] processing " << prevaread << " complete, kept " << written << " / " << VOVL.size() << std::endl;

	RIV.resize(0);
	VOVL.resize(0);
}

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

static double getDefaultTermVal()
{
	return 0.35;
}

int lasfilteralignmentsborderrepeats(libmaus2::util::ArgParser const & arg)
{
	double const termval = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : getDefaultTermVal();
	int64_t const lthres = 1000;
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);

	std::string const outfn = arg[0];
	std::string const dbfn = arg[1];
	std::string const repfn = arg[2];

	std::vector < std::string > Vinfn;
	for ( uint64_t i = 3; i < arg.size(); ++i )
		Vinfn.push_back(arg[i]);
	// std::string const infn = arg[3];

	for ( uint64_t i = 0; i < Vinfn.size(); ++i )
	{
		std::string const infn = Vinfn[i];
		std::string const lasindexname = libmaus2::dazzler::align::OverlapIndexer::getIndexName(infn);
		std::string const dalignerlasindexname = libmaus2::dazzler::align::DalignerIndexDecoder::getDalignerIndexName(infn);

		if (
			! libmaus2::util::GetFileSize::fileExists(lasindexname)
			||
			libmaus2::util::GetFileSize::isOlder(lasindexname,infn)
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(infn,&std::cerr);
		}

		if (
			! libmaus2::util::GetFileSize::fileExists(dalignerlasindexname)
			||
			libmaus2::util::GetFileSize::isOlder(dalignerlasindexname,infn)
		)
		{
			libmaus2::dazzler::align::OverlapIndexer::constructIndex(infn,&std::cerr);
		}
	}


	libmaus2::autoarray::AutoArray<Repeat> const R = Repeat::loadArray(repfn);

	std::cerr << "[V] copying " << dbfn << " to memory...";
	std::vector<std::string> Vtracks;
	Vtracks.push_back("inqual");
	libmaus2::dazzler::db::DatabaseFile::DBFileSet::unique_ptr_type dbptr(libmaus2::dazzler::db::DatabaseFile::copyToPrefix(dbfn,"mem:db1prefix",true /*regtmp */,&Vtracks));
	std::cerr << "done." << std::endl;
	libmaus2::dazzler::db::DatabaseFile::unique_ptr_type PDB(new libmaus2::dazzler::db::DatabaseFile(dbptr->fn));
	libmaus2::dazzler::db::DatabaseFile & DB = *PDB;
	if ( DB.part != 0 )
	{
		std::cerr << "Partial databases are not supported." << std::endl;
		return EXIT_FAILURE;
	}

	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getTSpace(Vinfn);
	DB.computeTrimVector();

	std::vector<uint64_t> RL;
	DB.getAllReadLengths(RL);

	libmaus2::dazzler::align::LasIntervals LAI(Vinfn,DB.size(),std::cerr);
	std::pair<int64_t,int64_t> const LAIint = LAI.getInterval();

	libmaus2::dazzler::db::Track::unique_ptr_type const Ptrack(DB.readTrack("inqual",0));
	libmaus2::dazzler::align::OverlapProperCheck const OPC(RL,*Ptrack,termval);

	libmaus2::autoarray::AutoArray < IsProper::unique_ptr_type > Aisproper(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		IsProper::unique_ptr_type tptr(new IsProper(DB,OPC,tspace));
		Aisproper[i] = UNIQUE_PTR_MOVE(tptr);
	}

	libmaus2::dazzler::align::AlignmentWriterArray::unique_ptr_type AWA(new libmaus2::dazzler::align::AlignmentWriterArray(tmpfilebase + "_simrep_tmp", numthreads, tspace));

	std::vector < std::string > Vsymkillfn(numthreads);
	libmaus2::autoarray::AutoArray < libmaus2::aio::OutputStreamInstance::unique_ptr_type > Asymkill(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		std::ostringstream ostr;
		ostr << tmpfilebase << "_symkill_" << i << ".tmp";
		std::string const fn = ostr.str();
		Vsymkillfn[i] = fn;
		libmaus2::util::TempFileRemovalContainer::addTempFile(fn);
		libmaus2::aio::OutputStreamInstance::unique_ptr_type tptr(
			new libmaus2::aio::OutputStreamInstance(fn)
		);
		Asymkill[i] = UNIQUE_PTR_MOVE(tptr);
	}

	int64_t minaread = LAIint.first;
	int64_t maxaread = LAIint.second;

	if ( maxaread < minaread )
	{
		AWA->merge(outfn,tmpfilebase + "_simrep_mergetmp");
		return EXIT_SUCCESS;
	}

	if ( arg.uniqueArgPresent("J") )
	{

		std::string const Js = arg["J"];
		std::istringstream istr(Js);
		int64_t Icnt;
		istr >> Icnt;

		if ( ! istr )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Js << std::endl;
			lme.finish();
			throw lme;
		}

		int const c = istr.get();

		if ( ! istr || c == std::istream::traits_type::eof() || c != ',' )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Js << std::endl;
			lme.finish();
			throw lme;
		}

		int64_t Idiv;
		istr >> Idiv;

		if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Js << std::endl;
			lme.finish();
			throw lme;
		}

		int64_t const toparead = maxaread + 1;
		int64_t const readspan = (toparead > minaread) ? (toparead-minaread) : 0;

		if ( readspan && ! Idiv )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] denominator of J argument cannot be zero" << std::endl;
			lme.finish();
			throw lme;
		}

		if ( toparead > minaread )
		{
			int64_t const partsize = Idiv ? (readspan + Idiv - 1)/Idiv : 0;

			int64_t const ilow = std::min(minaread + Icnt * partsize,toparead);
			int64_t const ihigh = std::min(ilow+partsize,toparead);

			if ( ihigh > ilow )
			{
				minaread = ilow;
				maxaread = ihigh-1;
			}
			else
			{
				minaread = 0;
				maxaread = -1;
			}
		}
	}
	else if ( arg.uniqueArgPresent("I") )
	{
		std::string const Is = arg["I"];
		std::istringstream istr(Is);
		int64_t Iminaread;
		istr >> Iminaread;

		if ( ! istr )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		int const c = istr.get();

		if ( ! istr || c == std::istream::traits_type::eof() || c != ',' )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		int64_t Imaxaread;
		istr >> Imaxaread;

		if ( ! istr || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "[E] unable to parse " << Is << std::endl;
			lme.finish();
			throw lme;
		}

		minaread = std::max(Iminaread,minaread);
		maxaread = std::min(Imaxaread,maxaread);
	}

	int64_t const toparead = maxaread+1;

	std::cerr << "[V] handling [" << minaread << "," << toparead << ")" << std::endl;

	int volatile gfailed = 0;
	libmaus2::parallel::PosixSpinLock gfailedlock;

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( int64_t r = minaread; r < toparead; ++r )
	{
		try
		{
			#if defined(_OPENMP)
			uint64_t const tid = omp_get_thread_num();
			#else
			uint64_t const tid = 0;
			#endif

			IsProper & isproper = *(Aisproper[tid]);
			libmaus2::dazzler::align::AlignmentWriter & AW = (*AWA)[tid];
			std::ostream & symkill = *(Asymkill[tid]);

			// libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Plas(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileRegion(infn,r,r+1));
			libmaus2::dazzler::align::AlignmentFileCat::unique_ptr_type Plas(LAI.openRange(r,r+1));
			libmaus2::dazzler::align::Overlap OVL;

			std::vector < libmaus2::dazzler::align::Overlap > VOVL;

			while ( Plas->getNextOverlap(OVL) )
				VOVL.push_back(OVL);

			if ( VOVL.size() )
			{
				// check whether A read has any marked repeats
				RepeatIdComparator const rcomp;
				std::pair<Repeat const *, Repeat const *> const P = ::std::equal_range(R.begin(),R.end(),Repeat(VOVL.front().aread),rcomp);
				std::vector < libmaus2::math::IntegerInterval<int64_t> > RIV;
				for ( Repeat const * p = P.first; p != P.second; ++p )
				{
					assert ( static_cast<int64_t>(p->id) == VOVL.front().aread );
					RIV.push_back(libmaus2::math::IntegerInterval<int64_t>(p->abpos,p->aepos-1));
				}

				std::ostringstream err;
				handle(err,R,DB,OPC,isproper,RIV,VOVL,VOVL.back().aread,lthres,AW,symkill,tspace,termval);

				{
					libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
					std::cerr << err.str();
				}
			}
		}
		catch(std::exception const & ex)
		{
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << ex.what() << std::endl;
			}

			gfailedlock.lock();
			gfailed = 1;
			gfailedlock.unlock();
		}
	}

	if ( gfailed )
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << "[E] parallel processing loop failed" << std::endl;
		lme.finish();
		throw lme;
	}

	// AWA->merge(outfn,tmpfilebase + "_simrep_mergetmp",libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::getDefaultMergeFanIn(),1);
	AWA->merge(outfn,tmpfilebase + "_simrep_mergetmp",libmaus2::dazzler::align::SortingOverlapOutputBuffer<>::getDefaultMergeFanIn(),numthreads);

	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		Asymkill[i]->flush();
		Asymkill[i].reset();
	}

	// merge sym files
	std::string const symkillmerge = tmpfilebase + "_merge_symkill.tmp";
	libmaus2::util::TempFileRemovalContainer::addTempFile(symkillmerge);

	libmaus2::sorting::SerialisingSortingBufferedOutputFile<libmaus2::dazzler::align::OverlapInfo>::reduce(Vsymkillfn,symkillmerge);

	for ( uint64_t i = 0; i < numthreads; ++i )
		libmaus2::aio::FileRemoval::removeFile(Vsymkillfn[i]);

	libmaus2::aio::OutputStreamFactoryContainer::rename(symkillmerge,outfn + ".sym");
	libmaus2::dazzler::align::OverlapInfoIndexer::createInfoIndex(outfn + ".sym",DB.size());
	libmaus2::dazzler::align::OverlapIndexer::constructIndex(outfn);

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
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("I", formatRHS("reads interval",std::string("0," + libmaus2::util::NumberSerialisation::formatNumber(std::numeric_limits<uint64_t>::max(),0)))));
	optionMap . push_back ( std::pair < std::string, std::string >("J", formatRHS("reads part",std::string("0,1"))));
	optionMap . push_back ( std::pair < std::string, std::string >("e", formatRHS("proper termination error threshold",getDefaultTermVal())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("prefix for temporary files",libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname))));
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
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 3 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las in.db in.rep in.las\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return lasfilteralignmentsborderrepeats(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
