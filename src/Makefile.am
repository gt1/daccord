ACLOCAL_AMFLAGS=-I m4

noinst_HEADERS = \
	ActiveElement.hpp \
	ChainSet.hpp \
	ComputeOffsetLikely.hpp \
	DebruijnGraphBase.hpp \
	DebruijnGraphContainer.hpp \
	DebruijnGraph.hpp \
	DebruijnGraphInterface.hpp \
	DecodedReadContainer.hpp \
	DotProduct.hpp \
	HandleContext.hpp \
	HetThreshold.hpp \
	Links.hpp \
	Node.hpp \
	OffsetLikely.hpp \
	Repeat.hpp \
	RepeatIdComparator.hpp \
	TraceType.hpp

MANPAGES = 

man_MANS = ${MANPAGES}

EXTRA_DIST = ${MANPAGES}

bin_PROGRAMS = daccord computeintrinsicqv lasdetectsimplerepeats lasfilteralignments lasfilteralignmentsborderrepeats maftobam generateperfectpiles bamidrename \
	checklas checkconsensus sortfasta mapconstoraw fillfasta canutolas spikenoise computeextrinsicqv marktrue encodegraph \
	split_agr split_dis filterchains wgsimtobam rlastobam bamfilterlongest computeintrinsicqv2 filtersym wmap


daccord_SOURCES = daccord.cpp
daccord_LDADD = ${LIBMAUS2LIBS}
daccord_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
daccord_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

computeintrinsicqv_SOURCES = computeintrinsicqv.cpp
computeintrinsicqv_LDADD = ${LIBMAUS2LIBS}
computeintrinsicqv_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
computeintrinsicqv_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasdetectsimplerepeats_SOURCES = lasdetectsimplerepeats.cpp
lasdetectsimplerepeats_LDADD = ${LIBMAUS2LIBS}
lasdetectsimplerepeats_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
lasdetectsimplerepeats_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasfilteralignments_SOURCES = lasfilteralignments.cpp
lasfilteralignments_LDADD = ${LIBMAUS2LIBS}
lasfilteralignments_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
lasfilteralignments_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

lasfilteralignmentsborderrepeats_SOURCES = lasfilteralignmentsborderrepeats.cpp
lasfilteralignmentsborderrepeats_LDADD = ${LIBMAUS2LIBS}
lasfilteralignmentsborderrepeats_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
lasfilteralignmentsborderrepeats_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

maftobam_SOURCES = maftobam.cpp
maftobam_LDADD = ${LIBMAUS2LIBS}
maftobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
maftobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

generateperfectpiles_SOURCES = generateperfectpiles.cpp
generateperfectpiles_LDADD = ${LIBMAUS2LIBS}
generateperfectpiles_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
generateperfectpiles_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamidrename_SOURCES = bamidrename.cpp
bamidrename_LDADD = ${LIBMAUS2LIBS}
bamidrename_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
bamidrename_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

checklas_SOURCES = checklas.cpp
checklas_LDADD = ${LIBMAUS2LIBS}
checklas_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
checklas_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

checkconsensus_SOURCES = checkconsensus.cpp
checkconsensus_LDADD = ${LIBMAUS2LIBS}
checkconsensus_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
checkconsensus_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

sortfasta_SOURCES = sortfasta.cpp
sortfasta_LDADD = ${LIBMAUS2LIBS}
sortfasta_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
sortfasta_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

mapconstoraw_SOURCES = mapconstoraw.cpp
mapconstoraw_LDADD = ${LIBMAUS2LIBS}
mapconstoraw_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
mapconstoraw_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

fillfasta_SOURCES = fillfasta.cpp
fillfasta_LDADD = ${LIBMAUS2LIBS}
fillfasta_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
fillfasta_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

canutolas_SOURCES = canutolas.cpp
canutolas_LDADD = ${LIBMAUS2LIBS}
canutolas_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
canutolas_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

spikenoise_SOURCES = spikenoise.cpp
spikenoise_LDADD = ${LIBMAUS2LIBS}
spikenoise_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
spikenoise_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

computeextrinsicqv_SOURCES = computeextrinsicqv.cpp
computeextrinsicqv_LDADD = ${LIBMAUS2LIBS}
computeextrinsicqv_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
computeextrinsicqv_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

marktrue_SOURCES = marktrue.cpp
marktrue_LDADD = ${LIBMAUS2LIBS}
marktrue_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
marktrue_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

encodegraph_SOURCES = encodegraph.cpp
encodegraph_LDADD = ${LIBMAUS2LIBS}
encodegraph_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
encodegraph_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

split_agr_SOURCES = split_agr.cpp
split_agr_LDADD = ${LIBMAUS2LIBS}
split_agr_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
split_agr_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

split_dis_SOURCES = split_dis.cpp
split_dis_LDADD = ${LIBMAUS2LIBS}
split_dis_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
split_dis_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

filterchains_SOURCES = filterchains.cpp
filterchains_LDADD = ${LIBMAUS2LIBS}
filterchains_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
filterchains_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

wgsimtobam_SOURCES = wgsimtobam.cpp
wgsimtobam_LDADD = ${LIBMAUS2LIBS}
wgsimtobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
wgsimtobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

rlastobam_SOURCES = rlastobam.cpp
rlastobam_LDADD = ${LIBMAUS2LIBS}
rlastobam_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
rlastobam_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

bamfilterlongest_SOURCES = bamfilterlongest.cpp
bamfilterlongest_LDADD = ${LIBMAUS2LIBS}
bamfilterlongest_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
bamfilterlongest_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

computeintrinsicqv2_SOURCES = computeintrinsicqv2.cpp
computeintrinsicqv2_LDADD = ${LIBMAUS2LIBS}
computeintrinsicqv2_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
computeintrinsicqv2_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

filtersym_SOURCES = filtersym.cpp
filtersym_LDADD = ${LIBMAUS2LIBS}
filtersym_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
filtersym_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}

wmap_SOURCES = wmap.cpp
wmap_LDADD = ${LIBMAUS2LIBS}
wmap_LDFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS} ${LIBMAUS2LDFLAGS}
wmap_CPPFLAGS = ${AM_CPPFLAGS} ${LIBMAUS2CPPFLAGS}
