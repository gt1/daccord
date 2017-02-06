#! /bin/bash
VERSION=`grep AC_INIT < configure.ac | awk -F',' '{print $2}'`
FIRST=`echo $VERSION | awk -F'.' '{print $1}'`
SECOND=`echo $VERSION | awk -F'.' '{print $2}'`
THIRD=`echo $VERSION | awk -F'.' '{print $3}'`
NEXTTHIRD=`expr ${THIRD} + 1`
export DEBEMAIL=tischler@mpi-cbg.de
export DEBFULLNAME="German Tischler"

function cleanup
{
	if [ ! -z "${COMMITFILE}" ] ; then
		if [ -f "${COMMITFILE}" ] ; then
			rm -f "${COMMITFILE}"
		fi
	fi
}

COMMITFILE=commit_msg_$$.txt

trap cleanup EXIT SIGINT SIGTERM

# make sure we have the latest version
git pull

# create commit log message
joe "${COMMITFILE}"

if [ ! -s "${COMMITFILE}" ] ; then
	echo "Empty commit log, aborting"
	exit 1
fi

# update to next minor version
awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^AC_INIT/ {gsub(first"."second"."third,first"."second"."third+1);print} ; !/^AC_INIT/{print}' < configure.ac > configure.ac.tmp
mv configure.ac.tmp configure.ac

# update change log
CHANGELOG=ChangeLog dch --distribution unstable -v ${FIRST}.${SECOND}.${NEXTTHIRD}-1

# commit files
git add configure.ac ChangeLog

git commit -F "${COMMITFILE}"
git push

TAG=daccord_${FIRST}_${SECOND}_${NEXTTHIRD}
git tag -a ${TAG} -m "daccord version ${FIRST}_${SECOND}_${NEXTTHIRD}"
git push origin ${TAG}

git checkout master
VERSION=`grep <configure.ac "AC_INIT" | perl -p -e "s/.*AC_INIT\(//" | awk -F ',' '{print $2}'`
DATE=`date +"%Y%m%d%H%M%S"`
RELEASE=${VERSION}-release-${DATE}
git checkout -b ${RELEASE}-branch master
autoreconf -i -f
ADDFILES="INSTALL Makefile.in aclocal.m4 autom4te.cache compile config.guess config.h.in config.sub configure depcomp install-sh ltmain.sh m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4 missing src/Makefile.in"
mv .gitignore .gitignore_
git add ${ADDFILES}
git commit -m "Release ${RELEASE}"
mv .gitignore_ .gitignore
git tag ${RELEASE}
git push origin ${RELEASE}
git checkout master
git branch -D ${RELEASE}-branch

exit 0
