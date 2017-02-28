AC_INIT(daccord,0.0.3,[tischler@mpi-cbg.de],[daccord],[http://www.mpi-cbg.de])
AC_CANONICAL_SYSTEM
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE([std-options subdir-objects parallel-tests])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

AC_CANONICAL_SYSTEM

AC_ARG_ENABLE(optimization,
        AS_HELP_STRING([--enable-optimization],[use compiler optimization (default yes)]),
        [optimization=${enableval}],[optimization=yes])
AC_ARG_ENABLE(debug,
        AS_HELP_STRING([--enable-debug],[use compiler debug flags (default no)]),
        [debug=${enableval}],[debug=no])
AC_ARG_ENABLE(profile,
        AS_HELP_STRING([--enable-profile],[use compiler profiling flags (default no)]),
        [profile=${enableval}],[profile=no])
AC_ARG_ENABLE(fast,
        AS_HELP_STRING([--enable-fast],[disable evaluation of assertions (default no)]),
        [fast=${enableval}],[fast=no])

# check for C compiler
CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
AC_PROG_CC
CFLAGS="${CFLAGS_SAVE}"

CFLAGS_SAVE="${CFLAGS}"
CFLAGS=""
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AM_PROG_CC_C_O
CFLAGS="${CFLAGS_SAVE}"
CXXFLAGS="${CXXFLAGS_SAVE}"

# check for C++ compiler
CXXFLAGS_SAVE="${CXXFLAGS}"
CXXFLAGS=""
AC_PROG_CXX
CXXFLAGS="${CXXFLAGS_SAVE}"

if test "${debug}" = "yes" ; then
	CXXDEBUGFLAGS="-g -O0"
	CDEBUGFLAGS="-g -O0"

	case ${CC} in
		gcc)
			CXXDEBUGFLAGS="${CXXDEBUGFLAGS} -rdynamic"
			CDEBUGFLAGS="${CDEBUGFLAGS} -rdynamic"
			;;
		*-mingw32msvc-gcc)
			CXXDEBUGFLAGS="${CXXDEBUGFLAGS} -rdynamic"
			CDEBUGFLAGS="${CDEBUGFLAGS} -rdynamic"
			;;
	esac

        CXXFLAGS="${CXXFLAGS} ${CXXDEBUGFLAGS}"
        CFLAGS="${CFLAGS} ${CDEBUGFLAGS}"
	openmp="no"

	AC_LANG_PUSH([C++])
	CXXFLAGS_SAVE="${CXXFLAGS}"
	CXXFLAGS="-fsanitize=address"
	AC_MSG_CHECKING([whether the C++ compiler supports the -fsanitize=address flag])
	AC_TRY_COMPILE([],[],[sanitizeaddress=yes],[sanitizeaddress=no])
	AC_MSG_RESULT([${sanitizeaddress}])
	CXXFLAGS="${CXXFLAGS_SAVE}"
	AC_LANG_POP

	if test "${sanitizeaddress}" = "yes" ; then
		CXXFLAGS="${CXXFLAGS} -fsanitize=address"
		CFLAGS="${CFLAGS} -fsanitize=address"
	fi

	AC_LANG_PUSH([C++])
	CXXFLAGS_SAVE="${CXXFLAGS}"
	CXXFLAGS="-fstack-protector"
	AC_MSG_CHECKING([whether the C++ compiler supports the -fstack-protector flag])
	AC_TRY_COMPILE([],[],[stackprotect=yes],[stackprotect=no])
	AC_MSG_RESULT([${stackprotect}])
	CXXFLAGS="${CXXFLAGS_SAVE}"
	AC_LANG_POP

	if test "${stackprotect}" = "yes" ; then
		STACKPROTECT="-fstack-protector"
		CXXFLAGS="${CXXFLAGS} ${STACKPROTECT}"
		CFLAGS="${CFLAGS} ${STACKPROTECT}"
	fi


	AC_MSG_NOTICE([Adding debug options ${CDEBUGFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
	AC_MSG_NOTICE([Adding debug options ${CXXDEBUGFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
else
	if test "${profile}" = "yes" ; then
	        CXXPROFFLAGS="-g -pg"
	        CPROFFLAGS="-g -pg"

	        CXXFLAGS="${CXXPROFFLAGS} -O3 -rdynamic"
	        CFLAGS="${CPROFFLAGS} -O3 -rdynamic"

		AC_MSG_NOTICE([Adding profiling options ${CPROFFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])
		AC_MSG_NOTICE([Adding profiling options ${CXXPROFFLAGS} to C compiler ${CXX} flags, full list now ${CXXFLAGS}])
	else
		if test "${optimization}" = "yes" ; then
			COPTFLAGS=
			case ${CC} in
				gcc)
					COPTFLAGS="-O3 -g -rdynamic"
					;;
				*-mingw32msvc-gcc)
					COPTFLAGS="-O3 -g -rdynamic"
					;;
				cl.exe)
					COPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CFLAGS="${CFLAGS} ${COPTFLAGS}"
			AC_MSG_NOTICE([Adding optimisation options ${COPTFLAGS} to C compiler ${CC} flags, full list now ${CFLAGS}])

			CXXOPTFLAGS=
			case ${CXX} in
				g++)
					CXXOPTFLAGS="-O3 -g -rdynamic"
					;;
				*-mingw32msvc-g++)
					CXXOPTFLAGS="-O3 -g -rdynamic"
					;;
				cl.exe)
					CXXOPTFLAGS="-O2 -Ob2 -Ot -Oy"
					;;
			esac
			CXXFLAGS="${CXXFLAGS} ${CXXOPTFLAGS}"

			AC_MSG_NOTICE([Adding optimisation options ${CXXOPTFLAGS} to C++ compiler ${CXX} flags, full list now ${CXXFLAGS}])
		fi
	fi
fi


# check for pkg-config tool
PKG_PROG_PKG_CONFIG

if test "${fast}" = "yes" ; then
	CPPFLAGS="${CPPFLAGS} -DNDEBUG"
fi

# check for standard headers
AC_LANG_PUSH([C++])
AC_HEADER_STDC
AC_LANG_POP

if test "${CXX}" = "g++" ; then
	CPPFLAGS="${CPPFLAGS} -W -Wall"
fi

AC_ARG_WITH([libmaus2],
            [AS_HELP_STRING([--with-libmaus2@<:@=PATH@:>@], [path to installed libmaus2 library @<:@default=@:>@])],
            [with_libmaus2=${withval}],
            [with_libmaus2=])

if test ! -z "${with_libmaus2}" ; then
	PKGCONFIGPATHSAVE="${PKG_CONFIG_PATH}"
	if test -z "${PKG_CONFIG_PATH}" ; then
		export PKG_CONFIG_PATH="${with_libmaus2}/lib/pkgconfig"
	else
		export PKG_CONFIG_PATH="${with_libmaus2}/lib/pkgconfig:${PKG_CONFIG_PATH}"
	fi
fi

PKG_CHECK_MODULES([libmaus2],[libmaus2 >= 2.0.303])

if test ! -z "${with_libmaus2}" ; then
	if test ! -z "${PKGCONFIGPATHSAVE}" ; then
		export PKG_CONFIG_PATH="${PKGCONFIGPATHSAVE}"
	fi
fi

LIBMAUS2CPPFLAGS="${libmaus2_CFLAGS}"
LIBMAUS2LIBS="${libmaus2_LIBS}"

CPPFLAGS_SAVE="${CPPFLAGS}"
LDFLAGS_SAVE="${LDFLAGS}"
LIBS_SAVE="${LIBS}"

CPPFLAGS="${CPPFLAGS} ${libmaus2_CFLAGS}"
LIBS="${LIBS} ${libmaus2_LIBS}"

AC_LANG_PUSH([C++])
AC_MSG_CHECKING([whether we can compile a program using libmaus2])
	AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <libmaus2/util/ArgInfo.hpp>

				int main(int argc, char * argv[[]]) {
					::libmaus2::util::ArgInfo const arginfo(argc,argv);
					return 0;
			}])],
			have_libmaus2=yes,
			have_libmaus2=no
		)
AC_LANG_POP
AC_MSG_RESULT($have_libmaus2)

if test "${have_libmaus2}" = "no" ; then
	AC_MSG_ERROR([Required libmaus2 is not available.])
fi

AC_MSG_NOTICE([Using flags ${CFLAGS} for C compiler ${CC}])
AC_MSG_NOTICE([Using flags ${CXXFLAGS} for C++ compiler ${CXX}])

PACKAGE_NAME=${PACKAGE}
PACKAGE_VERSION=${VERSION}

AC_SUBST([PACKAGE_NAME])
AC_SUBST([PACKAGE_VERSION])
# libmaus2
AC_SUBST([LIBMAUS2CPPFLAGS])
AC_SUBST([LIBMAUS2LIBS])
# 
AC_OUTPUT(Makefile src/Makefile)
