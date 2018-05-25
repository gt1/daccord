/*
    daccord
    Copyright (C) 2018 German Tischler-Höhle

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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/dazzler/db/DatabaseFile.hpp>
#include <libmaus2/fastx/FastaPeeker.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/lcs/NPLinMem.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/parallel/NumCpus.hpp>
#include <libmaus2/dazzler/align/SimpleOverlapParser.hpp>
#include <libmaus2/dazzler/align/OverlapDataInterface.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/dazzler/align/DalignerIndexDecoder.hpp>
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/ProgramHeaderLineSet.hpp>
#include <libmaus2/bambam/BamAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamMultiAlignmentDecoderFactory.hpp>
#include <libmaus2/bambam/BamPeeker.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>

#include <config.h>

static uint64_t getDefaultNumThreads()
{
	return libmaus2::parallel::NumCpus::getNumLogicalProcessors();
}

std::string getTmpFileBase(libmaus2::util::ArgParser const & arg)
{
	std::string const tmpfilebase = arg.uniqueArgPresent("T") ? arg["T"] : libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	return tmpfilebase;
}

::libmaus2::bambam::BamHeader::unique_ptr_type updateHeader(
	::libmaus2::util::ArgInfo const & arginfo,
	::libmaus2::bambam::BamHeader const & header
)
{
	std::string const headertext(header.text);

	// add PG line to header
	std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
		headertext,
		"bammerge", // ID
		"bammerge", // PN
		arginfo.commandline, // CL
		::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
		std::string(PACKAGE_VERSION) // VN
	);
	// construct new header
	::libmaus2::bambam::BamHeader::unique_ptr_type uphead(new ::libmaus2::bambam::BamHeader(upheadtext));

	return UNIQUE_PTR_MOVE(uphead);
}

struct SimpleThreadPoolTerminate
{
	libmaus2::parallel::SimpleThreadPool * STP;

	SimpleThreadPoolTerminate(libmaus2::parallel::SimpleThreadPool * rSTP)
	: STP(rSTP) {}

	~SimpleThreadPoolTerminate()
	{
		if ( STP )
		{
			STP->terminate();
			STP->join();
		}
	}
};

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);

		// number of threads
		uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();

		std::string const tmpfilebase = getTmpFileBase(arg);

		libmaus2::util::ArgInfo arginfo(argc,argv);
		{
			std::ostringstream ostr;
			ostr << numthreads;
			arginfo.replaceKey("inputthreads",ostr.str());
			arginfo.replaceKey("outputthreads",ostr.str());
		}

		libmaus2::parallel::SimpleThreadPool::unique_ptr_type PSTP;
		if ( numthreads > 1 )
		{
			libmaus2::parallel::SimpleThreadPool::unique_ptr_type TSTP(
				new libmaus2::parallel::SimpleThreadPool(numthreads)
			);
			PSTP = UNIQUE_PTR_MOVE(TSTP);

			libmaus2::bambam::BamAlignmentDecoderFactory::setThreadPool(PSTP.get());
		}

		SimpleThreadPoolTerminate STPT(PSTP ? PSTP.get() : 0);

		libmaus2::bambam::BamAlignmentDecoderWrapper::unique_ptr_type decwrapper(
			libmaus2::bambam::BamMultiAlignmentDecoderFactory::construct(
				arginfo,false, // do not put rank
				0, /* copy stream */
				std::cin, /* standard input */
				true, /* concatenate instead of merging */
				false /* streaming */
			)
		);
		::libmaus2::bambam::BamAlignmentDecoder * ppdec = &(decwrapper->getDecoder());
		::libmaus2::bambam::BamAlignmentDecoder & dec = *ppdec;
		::libmaus2::bambam::BamHeader const & header = dec.getHeader();

		std::string const headertext(header.text);

		// add PG line to header
		std::string const upheadtext = ::libmaus2::bambam::ProgramHeaderLineSet::addProgramLine(
			headertext,
			"bamfilterlongest", // ID
			"bamfilterlongest", // PN
			arg.commandline, // CL
			::libmaus2::bambam::ProgramHeaderLineSet(headertext).getLastIdInChain(), // PP
			std::string(PACKAGE_VERSION) // VN
		);
		// construct new header
		::libmaus2::bambam::BamHeader uphead(upheadtext);

		libmaus2::bambam::BamPeeker BP(dec);

		libmaus2::bambam::BamAlignment algn;

		std::vector< ::libmaus2::lz::BgzfDeflateOutputCallback * > * Pcbs = 0;
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(
			libmaus2::bambam::BamBlockWriterBaseFactory::construct(uphead,arginfo,Pcbs));

		while ( BP.peekNext(algn) )
		{
			int64_t const refid = algn.getRefID();

			std::vector < libmaus2::bambam::BamAlignment::shared_ptr_type > V;

			std::map < std::string, uint64_t> M;

			while ( BP.peekNext(algn) && algn.getRefID() == refid )
			{
				BP.getNext(algn);
				libmaus2::bambam::BamAlignment::shared_ptr_type sptr(algn.sclone());
				V.push_back(sptr);

				std::string const name = algn.getName();

				if ( M.find(name) == M.end() )
				{
					M[name] = V.size()-1;
				}
				else
				{
					uint64_t const previd = M.find(name)->second;

					if ( algn.getReferenceLength() > V[previd]->getReferenceLength() )
						M[name] = V.size()-1;
				}
			}

			for ( std::map < std::string, uint64_t>::const_iterator it = M.begin();
				it != M.end(); ++it )
			{
				Pwriter->writeAlignment(*V[it->second]);
			}
		}

		Pwriter.reset();
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
