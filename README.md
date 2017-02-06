d'accord
========

d'accord is a non hybrid long read consensus program based on local de Bruijn graph assembly. The package contains the following programs

 - daccord: main consensus computation program
 - computeintrinsicqv: compute intrinsic quality value track
 - lasdetectsimplerepeats: simple repeat detection for repeats completed contained inside a long read
 - lasfilteralignments: remove improper alignments (likely repeat induced  alignments)
 - lasfilteralignmentsborderrepeats: remove repeat induced alignments at read borders 
 - maftobam: convert MAF files produced by PBsim to BAM
 - bamidrename: replace read names in BAM file by unique identifiers
 - generateperfectpile: generate perfect alignment pile based on read to reference alignments
 - checklas: compare aligner generated data with perfect alignment piles
 - checkconsensus: check consensus accurary based on ground truth data
 - sortfasta: sort a FastA file by read name
 - mapconstoraw: insert corrected fragments into uncorrected reads
 - fillfasta: insert missing reads into FastA file

A short list of options is available for each program by calling it
with the -h parameter, e.g.

	daccord -h

Source
------

The daccord source code is hosted on github:

	git@github.com:gt1/daccord.git

Release packages can be found at

	https://github.com/gt1/daccord/releases

Please make sure to choose a package containing the word "release" in it's name if you
intend to compile daccord for production (i.e. non development) use.

Compilation of daccord
----------------------

daccord needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then daccord can be compiled and
installed in ${HOME}/daccord using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/daccord
	- make install

The release packages come with a configure script included (making the autoreconf call unnecessary for source obtained via one of those).

daccord
-------

daccord is the main consensus computation program. It requires the name of an LAS file and a Dazzler database as input arguments:

	daccord reads.las reads.db

Please see DAZZ_DB (https://github.com/thegenemyers/DAZZ_DB) for creating Dazzler databases and
DALIGNER (https://github.com/thegenemyers/DALIGNER) for creating LAS alignment files.

daccord has a set of optional parameters which can be set on the command line:

 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -w: window size for building De Bruijn graphs
 * -a: advance size for De Bruijn graph windows. Window i is at positions [i*a,i*a+w)
 * -d: maximum depth used. If set to some v then daccord will use only the v longest alignments for each read
 * -f: produce full alignments. Set by calling with -f1. If set then daccord will output one record per corrected read.
       Corrected stretches are printed as upper case letters, uncorrected stretches as lower case letters.
 * -V: verbosity. Smaller values make the program more quiet.
 * -I: read interval. Given as -Ii,j . This means the program will process the read id interval [i,j).
 * -J: batch interval. Given as -Ji,j for 0 <= i < j. The program will process batch id i out of j packages.
       More precisely, if the input LAS file starts at A-read id a and ends at A-read z, then each batch package
       has size s=(z-a+j-1)/j and batch package i contains reads [a+i*s,a+(i+1)*s)
 * -E: sets the name of the error profile to be used. By default this is the name of the input LAS file with .eprof appended.
       If this file does not exist then the program will estime the error profile based on the alignments in the LAS file.
 * -m: minimum window coverage. This denotes the minimum number of reads which need to cover a window so consensus construction
       will be attempted.
 * -e: maximum window error. If this is set then the windows which have an average error between the reads and the consensus of
       more than this will not be used for the consensus.
 * -l: minimum length of output. Corrected stretches shorter than this length will not be produced. This is ignored if the -f parameter is set.

Please note that the parameters need to be provided before any other arguments (i.e. before the LAS file name and the name of the Dazzler DB).

daccord produces a FastA file on standard output. The read names follow the
scheme

	read-id/unique-id/from_to

An example for this is

	2/7/500_23141

This designates that the respective sequence is the corrected version of
read id 2 (in the trimmed Dazzler database) from base 500 to 23141. The
unique-id field provides a unique identifier for each corrected fragment and
has otherwise no specific meaning.

computeintrinsicqv
------------------

computeintrinsicqv computes intrinsic quality values and requires a numerical depth 
parameter -d, a Dazzler database and an LAS file to compute intrinsic quality values. An example call is:

	computeintrinsicqv -d40 reads.db reads.las

The argument for the depth parameter is the sequencing depth, which is given by the number of sequenced bases
divided by the size of the genome sequenced. The reads file can refer to a single block of the database, e.g.

	computeintrinsicqv -d40 reads.db reads.1.las

is a valid call for the first block of reads.db . The database however needs to be complete, i.e.

	computeintrinsicqv -d40 reads.1 reads.1.las

is not a valid call. When intrinsic quality values have been computed for all blocks of a database, then the
quality tracks for the single blocks can be concatenated to the quality track of the complete database using
the Catrack tool of DAZZ_DB.

lasdetectsimplerepeats
----------------------

lasdetectsimplerepeats is a program which detects simple repeats (such which are completely contained in a read)
by virtue of checking intrinsic quality values. The program takes a Dazzler database and an LAS file as input.
An example call is

	lasdetectsimplerepeats reads.db reads.las >reads.rep

The database given has to be complete, the LAS file can refer to a single block of the database.

lasdetectsimplerepeats has a set of optional parameters which can be set on the command line:

 * -e: error threshold for proper alignment ends. The end of an alignment is considered proper if it either extends to a read
       end or to a region with an intrinsic quality corresponding to an error of at least this value. By default
       this value is set to 0.35, i.e 35%.
 * -d: depth threshold for repeat detection. A region of a read will be considered as a repeat if at least this many
       other reads give evidence that the region is a repeat. By default this is set to 10, i.e. 10.

lasdetectsimplerepeats produces it's output on the standard output channel. The output for single blocks can be concatenated
to the output for the whole database using the Unix cat program. This requires the outputs of lasdetectsimplerepeats
to be concatenated in increasing order of the block id, i.e. run

	cat reads.1.rep reads.2.rep ... reads.n.rep >reads.rep

if n is the number of blocks in the database.

lasfilteralignments
-------------------

lasfilteralignments is a program which removes improper alignments from an LAS file. It requires intrinsic qualities to do so.
The program takes a Dazzler database and an LAS file as input.
An example call is

	lasfilteralignments reads.db reads.las

This will create the output file reads_filtered.las .

lasfilteralignments has a set of optional parameters which can be set on the command line:

 * -e: error threshold for proper alignment ends. The end of an alignment is considered proper if it either extends to a read
       end or to a region with an intrinsic quality corresponding to an error of at least this value. By default
       this value is set to 0.35, i.e 35%.

lasfilteralignmentsborderrepeats
--------------------------------

lasfilteralignmentsborderrepeats is a program which removes (probably)
repeat induced alignments involving a prefix or suffix of a read. The
repeats are not detected by the program but need to be provided in the form
of an output file as produced by lasdetectsimplerepeats. It expect four
input arguments:

 * an output file name (output in LAS format)
 * a Dazzler database
 * a repeats file (as produced by lasdetectsimplerepeats)
 * an input LAS file name

An example call is

	lasfilteralignmentsborderrepeats out.las in.db in.rep in.las

lasfilteralignmentsborderrepeats has a set of optional parameters which can be set on the command line:

 * -e: error threshold for proper alignment ends. The end of an alignment is considered proper if it either extends to a read
       end or to a region with an intrinsic quality corresponding to an error of at least this value. By default
       this value is set to 0.35, i.e 35%.
 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -I: read interval. Given as -Ii,j . This means the program will process the read id interval [i,j).
 * -J: batch interval. Given as -Ji,j for 0 <= i < j. The program will process batch id i out of j packages.
       More precisely, if the input LAS file starts at A-read id a and ends at A-read z, then each batch package
       has size s=(z-a+j-1)/j and batch package i contains reads [a+i*s,a+(i+1)*s)
 * -T: prefix for temporary files. By default temporary files will be
       created in the current directory. The program tries to create file
       names avoiding collisions with other program runs based on the name
       of the program, the host name it runs on and the process id.

maftobam
--------

maftobam converts alignments produced in the MAF format by PBsim to the BAM
format. It requires one argument designating a FastA file containing the
underlying reference. An example call is

	maftobam ref.fasta <in.maf >out.bam

An optional second file name can be given to perform online reference id
replacements. This is useful as PBsim does not output the name of the
original reference sequence.

bamidrename
-----------

bamidrename processes a BAM file an replaces all read names after the scheme

	L0/id/0_len

where id is the index of the alignment record (i.e. 0 for the first one, 1
for the second, etc.) and len is the length of the query sequence. This
allows to identify the index of a record in the output file
after the file has later been sorted to a different order.

generateperfectpiles
--------------------

generateperfectpiles generates alignments contained in perfect alignment piles
from a set of read to reference alignments given in a BAM file. The alignments
produced are output in DALIGNER's LAS format. The input BAM file is assumed
to have geen generated by first running bamidrename on it and subsequently
having been sorted by coordinate order (use for instance biobambam2's
bamsort utility to obtain this order).

A sample call is

	generateperfectpiles out.las in.bam

generateperfectpiles has a set of optional parameters which can be set on the command line:

 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -T: prefix for temporary files. By default temporary files will be
       created in the current directory. The program tries to create file
       names avoiding collisions with other program runs based on the name
       of the program, the host name it runs on and the process id.

checklas
--------

checklas compares an aligner produced LAS file with an LAS file containing
perfect piles only and prints a statstic comparing both for each read. A
sample call is

	checklas in.db perfectpiles.las in.bam in.las

It requires four arguments

 * in.db: a Dazzler database on which perfectpiles.las and in.las are based
 * perfectpiles.las: ground truth perfect alignment piles in LAS format
 * in.bam: ground truth reads to reference alignments
 * in.las: aligner generated data in LAS format

The file in.bam is expected to contain read names as produced by
bamidrename. It has to be sorted by read name using biobambam2's bamsort
utility.

Optional parameters are:

 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -T: prefix for temporary files. By default temporary files will be
   created in the current directory. The program tries to create file
   names avoiding collisions with other program runs based on the name
   of the program, the host name it runs on and the process id.
 * --verbose: a higher value leads to more verbosity (default: 0)
 * --minsiglen: minimum length of ground truth alignemtns considered (default 1000)
   The ground truth alignments file may contain arbitrarily
   short alignments because two reads overlap by a very small
   amount of bases on the reference. As such short alignments
   do not designate significant alignment events no real aligner will produce them.
 * --tlow: minimum read id considered (by default all reads in the input file are considered)
 * --high: maximum read id considered (by default all reads in the input file are considered)
 * --mark: produce file in.las.mark.las in which ground truth alignments are
   marked using flag 2^30, any other alignments in in.las are copied
   as they are (i.e. no flag is added)
 * --mis: produce file in.las.mis.las containing all ground truth read against read alignments
   which are not contained in in.las
 * --keep: produce file in.las.keep.las containing all ground truth read against read alignment which are contained in in.las
 * --seqdepth: average sequencing coverage

Statistic lines are printed on the standard output channel. A sample output
line is

	[R] 1941 missing 0 got 62 extra 0 erate 0.118001 U00096.2:409510,426938 0 1 1 1 180 regular

The columns are

 * line id ([R])
 * read id (1941)
 * string missing
 * number of missing ground truth alignments in in.las (0)
 * string got
 * number of ground truth alignments present in in.las (62)
 * string extra
 * number of extraneous alignments in in.las which are not true overlaps (0)
 * string erate
 * error rate of ground truth alignment (0.118001)
 * alignment coordinates of ground truth mapping in format refid:from,to (U00096.2:409510,426938)
 * minimum true rate of any trace block on read
 * average true rate of any trace block on read
 * maximum true rate of any trace block on read
 * number of blocks average is based on
 * line category "regular" or complete miss (regular). A read is marked as
   regular if in.las contains any alignments for it, no matter whether
   they are contained in the ground truth set or not. Otherwise it is a complete miss.

The true rate values are computed as follows. For each trace block (usually
of size 100) the (at most) seqdepth best (i.e. lowest number of errors) aligning reads
are collected. Then the number n_t of true (contained in the ground truth) and
n_f of false alignments in this list of seqdepth are counted. Based on this
a fraction f = n_t / (n_t+n_f) is computed which designates how many of the
(at most) seqdepth minimum error alignments are true. The program then
outputs the minimum, average and maximum over all trace blocks for the read.
In essence the minimum column should be 1 for any reads not stemming from
repetetive regions.

checkconsensus
--------------

checkconsensus checks consensus accuracy by comparing consensus sequences to
ground truth data. The program expects four arguments

 * a sub command argument. The available sub commands are index, check,
   batchlist. These sub commands are explained below.
 * ref.fasta: a FastA file containing the reference sequences
 * reads.bam: a BAM file containing ground truth read to reference alignments
 * reads_cons.fasta: a FastA file containing the consensus sequences

The names of the sequences in ref.fasta must match the sequence lines in the
header of reads.bam and the sequences need to appear in the same order in
ref.fasta and the header of reads.bam. The read names in reads.bam need to
follow the scheme described for bamidrename. reads.bam needs to have been sorted to
query name order using biobambam2's bamsort. The read names in
reads_cons.fasta need to follow the scheme used in the about of daccord as
described above.

Before any other sub commands can be run the program needs to be run with
teh index sub command. This computes indices for the various input files
which will be used by the other sub command.

If the program is to be run on a single node, then the check sub command can
be used subsequently to perform the accuracy comparison. Optional arguments
are

 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -T: prefix for temporary files. By default temporary files will be created in the current directory. The program tries to create file
   names avoiding collisions with other program runs based on the name
   of the program, the host name it runs on and the process id.

For larger read sets the checking can be distributed over the nodes of a
compute cluster. To this end the program can be run with the batchlist
sub command. This will produce a shell script containing program calls which
need to be called to obtain the same result as if the check sub command
would have been called. An example call is

	checkconsensus -p4 batchlist ref.fasta reads.bam reads_cons.fasta

This will produce a shell script similar to

	checkconsensus <opts> batchprocess ref.fasta reads.bam reads_cons.fasta reads_cons.fasta_check_000000
	checkconsensus <opts> batchprocess ref.fasta reads.bam reads_cons.fasta reads_cons.fasta_check_000001
	checkconsensus <opts> batchprocess ref.fasta reads.bam reads_cons.fasta reads_cons.fasta_check_000002
	checkconsensus <opts> batchprocess ref.fasta reads.bam reads_cons.fasta reads_cons.fasta_check_000003
	checkconsensus <opts> batchmerge ref.fasta reads.bam reads_cons.fasta
	checkconsensus <opts> cleanup ref.fasta reads.bam reads_cons.fasta

using the sub commands batchprocess, batchmerge and cleanup (which should
not be called directly). All the lines for batchprocess can be processed in
parallel. After the batchprocess lines have been processed the batchmerge
sub command will merge the resulting partial data and produce the output of
the accuracy check. The cleanup sub command removes partial files produced
by the batchprocess commands.

The batchlist sub command has the following optional parameters:

 * -t: sets the number of threads used by each sub command in the shell
   script produced. By default this is the number of logical CPU cores detected
 * -T: prefix for temporary files. By default temporary files will be created in the current directory. The program tries to create file
   names avoiding collisions with other program runs based on the name
   of the program, the host name it runs on and the process id.
 * -p: number of work packages produced, i.e. number of batchprocess calls
   in the resulting script.

The output of checkconsensus contains several types of lines:

 * single read lines
 * coverage lines
 * global statistics line(s)
 * EP lines
 * EM lines

An example for a single read line is

	6       AlignmentStatistics(matches=14732,mismatches=0,insertions=2,deletions=0,editdistance=2,erate=0.0001357405)      AlignmentStatistics(matches=112436,mismatches=0,insertions=10,deletions=0,editdistance=10,erate=0.0000889316)   [0,14731]       0:97757,112489

The columns are

 * read id (6)
 * accuracy of read fragment in comparison with ground truth (AlignmentStatistics(matches=14732,mismatches=0,insertions=2,deletions=0,editdistance=2,erate=0.0001357405))
 * accumulated average accuracy for this read and reads with smaller read id
 * interval on read corrected ([0,14731])
 * reference region corrected fragment was aligned to for comparison in the
   format refid:from-to (0:97757,112489)

An example for a coverage line is

	[C]     16      13645   13642   0.99978

The columns are

 * line identifier ([C])
 * read id (16)
 * length of ground truth in reference bases (13645)
 * number of bases covered by reconstruction fragments
 * column four / column three (i.e. fraction of ground truth bases covered
   by reconstruction fragments)

An example of a global statistics line is

	[G]     29460533        29447890        0.999571        AlignmentStatistics(matches=29447461,mismatches=195,insertions=2199,deletions=234,editdistance=2628,erate=0.0000892357)

The columns are

 * line identifier ([G])
 * reference base sum (29460533, sum over third columns of all coverage lines)
 * covered base sum (29447890, sum over fourth columns of all coverage lines)
 * column three / column two (0.999571, i.e. fraction of ground truth bases covered by reconstruction fragments)
 * accumulated error statistics over all single read lines  (AlignmentStatistics(matches=29447461,mismatches=195,insertions=2199,deletions=234,editdistance=2628,erate=0.0000892357))

An example of an EM line is

	[EM]    0.0004741958    18      0.0091047

The columns are

 * line identifier ([EM])
 * error rate (0.0004741958)
 * number of reads having this error rate compared to ground truth or higher (18)
 * fraction of reads having this error rate compare to ground truth or  higher (0.0091047)

An example of an EP line is

	[EP]    0.0000485201    669     0.338392

The columns are

 * line identifier ([EP])
 * error rate (0.0000485201)
 * number of reads having this error rate compared to ground truth or lower (669)
 * fraction of reads having this error rate compare to ground truth or lower (0.338392)

sortfasta
---------

sortfasta reads a FastA file from standard input, sorts it by read name and
outputs the sorted data on standard output. Read names are compared in the
same way as biobambam2's bamsort does for query name sorting.

mapconstoraw
------------

mapconstoraw inserts corrected fragments into uncorrected reads by mapping the fragments onto the reads.
It requires two arguments, both in FastA format:

 * the uncorrected reads
 * the corected read fragments

Both of these input files are required to have been sorted using the sortfasta program.

The read names in the uncorrected reads file need to follow the scheme

	read_{id}

where {id} is a numerical id. A valid example would be

	read_5

The read names in the corrected read fragments file need to follow the scheme

	read_{id}_{subid}

where {id} and {subid} are numerical values. A valid example would be

	read_5_0

which would designate the first fragment (sub id 0) of read 5. The sub ids for each read need to be
consecutive and start from 0.

mapconstoraw outputs a FastA file with the mapped corrected fragments inserted into the uncorrected reads. Corrected
fragments are in upper case, uncorrected regions in lower case.

mapconstoraw has the following optional parameters:

 * -t: sets the number of threads used. By default this is the number of logical CPU cores detected
 * -T: prefix for temporary files. By default temporary files will be created in the current directory. The program tries to create file
   names avoiding collisions with other program runs based on the name
   of the program, the host name it runs on and the process id.

fillfasta
---------

fillfasta processes two sorted FastA files and computes their union such that the second file takes precedence (i.e. if a read name
is present in both files, then the data from the second file will be kept and the one from the first one discarded).
The read names contained in the second file need to be a subset of the read names found in the first file.
The read names in both files need to follow the schem

	read_{id}

where {id} is a numerical id. A valid example would be

	read_5
