dnl Process this file with autoconf to produce a configure script.

AC_INIT([bioacoustics], 0.2.5) dnl package name, version

# require at least autoconf 2.62
AC_PREREQ(2.62)

dnl Now find the compiler and compiler flags to use
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

RBIN="${R_HOME}/bin/R"
RVER=$("${R_HOME}/bin/Rscript" -e 'cat((R.version)$major, (R.version)$minor, sep=".")')

: ${CC=$("${RBIN}" CMD config CC)}
: ${CFLAGS=$("${RBIN}" CMD config CFLAGS)}
: ${LDFLAGS=$("${RBIN}" CMD config LDFLAGS)}

if test [ -z "$CXX" ] ; then
	CXX11=CXX11
	AS_VERSION_COMPARE([${RVER}], [3.4.0], [CXX11=CXX1X], [], [])
	CXX11=$("${RBIN}" CMD config "${CXX11}")
fi

if test [ -z "$CXX11STD" ] ; then
	if test [ -z "$CXX1XSTD"]
	then
		CXX11STD=CXX11STD
		AS_VERSION_COMPARE([${RVER}], [3.4.0], [CXX11STD=CXX1XSTD], [], [])
		CXX11STD=$("${RBIN}" CMD config "${CXX11STD}")
	else
		CXX11STD=CXX1XSTD
	fi
fi

: ${CXXPICFLAGS=$("${RBIN}" CMD config CXXPICFLAGS)}
: ${CXX="${CXX11} ${CXX11STD} ${CXXPICFLAGS}"}


AC_PROG_CC
AC_PROG_CXX

#
# Check if compiler supports C++11 features
#
AX_CXX_COMPILE_STDCXX_11()

AC_MSG_NOTICE([CC: ${CC}])
AC_MSG_NOTICE([CXX: ${CXX}])

#
# Check for cmake
#
CMAKE_SEARCH_PATH="${PATH}:/Applications/CMake.app/Contents/bin"
AC_CHECK_PROG([CMAKE_FOUND], [cmake], [cmake], [not found], ${CMAKE_SEARCH_PATH})

if test x"${CMAKE_FOUND}" != x"cmake" ; then
  AC_MSG_ERROR(Could not find cmake)
fi

#
# Check for OpenMP
#
AC_LANG([C++])
AC_OPENMP
AC_LANG([C])
AC_OPENMP

if test x"${OPENMP_CFLAGS}" != "x"; then
	FFTW_OPENMP="--enable-openmp"
fi

#
# Check if MAKE is set, if not set to make
#
if test -z "${MAKE}"; then
  MAKE = make
fi


#
# fftw
#

dnl Select an optional path to the header file of the FFTW library
AC_ARG_WITH([fftw-include],
            AC_HELP_STRING([--with-fftw-include=FFTW_INCL],
                           [the location of fftw3 header file]),
            [fftw_incl_path=$withval])
if test [ -n "$fftw_incl_path" ] ; then
   FFTW_INCL="-I${fftw_incl_path}"
else
  if test [ -n "${LIB_FFTW}" ] ; then
     FFTW_INCL="-I${LIB_FFTW}/include"
  fi
fi

AC_CHECK_HEADERS(fftw3.h, [fftw_ok="yes"], [fftw_ok="no"])


# dnl ditto for a library path
AC_ARG_WITH([fftw-lib],
            AC_HELP_STRING([--with-fftw-lib=FFTW_LIB],
                           [the location of fftw3 library]),
            [fftw_lib_path=$withval])
if test [ -n "$fftw_lib_path" ] ; then
   FFTW_LIB="-L$fftw_lib_path"
else
  if test [ -n "${LIB_FFTW}" ] ; then
     FFTW_LIB="-L${LIB_FFTW}/lib"
  else
    if test -n "${FFTW_CONFIG}"; then
      fftw_lib_path=$(fftw_config --libs)
      FFTW_LIB="${fftw_lib_path}"
    fi
  fi
fi

AC_SEARCH_LIBS(fftw_plan_r2r_1d, fftw3, [fftw_ok="yes"], [fftw_ok="no"])

## If fftw3 is not found, download fftw3 and build it locally
if test x"${fftw_ok}" = x"no"; then
  AC_MSG_NOTICE([Need to download and build fftw3])

  ## define FFTW version
  FFTW_VERSION=3.3.8

  ## define FFTW file and download URL
  FFTW_TGZ="fftw-${FFTW_VERSION}.tar.gz"
  FFTW_URL="http://fftw.org/${FFTW_TGZ}"

  ## C Compiler options
  FFTW_CFLAGS=

  ## additional C Compiler options for linking
  FFTW_CLINKFLAGS=

  ## Libraries necessary to link with the package
  FFTW_LIB="-L\"$(pwd)/fftw-${FFTW_VERSION}/lib\""

  ## Necessary Include dirs
  FFTW_INCL="-I$(pwd)/fftw-${FFTW_VERSION}/include"

  ## Download fftw3 source code
  $("${R_HOME}/bin/Rscript" --vanilla -e "download.file(url='${FFTW_URL}', destfile='${FFTW_TGZ}')")

  ## Extract fftw3 source code and remove .tar.gz
  $("${R_HOME}/bin/Rscript" --vanilla -e "untar(tarfile='${FFTW_TGZ}', exdir='$(pwd)')")
  $(rm -rf ${FFTW_TGZ})

  ## Compile fftw3 source code and clean up
  AC_MSG_NOTICE([Starting to install fftw3 library to $(pwd)/fftw-${FFTW_VERSION}])

  $(set -e; \
    cd fftw-${FFTW_VERSION}; \
    ./configure --prefix=$(pwd) \
                --enable-static \
                --disable-shared \
                ${FFTW_OPENMP} \
                CFLAGS=-fPIC > /dev/null; \
    ${MAKE} > /dev/null; \
    ${MAKE} install > /dev/null;
  )

  FFTW3_INSTALL=$?

  if test x"${FFTW3_INSTALL}" != x"0" ; then
    AC_MSG_ERROR([Could not install fftw3 library, try to install it manually?])
  fi

  if [[ "$LIBS" != *"lfftw3"* ]]; then
      LIBS="$LIBS -lfftw3"
  fi

  AC_MSG_NOTICE([Done installing fftw3 library])

else
  AC_MSG_NOTICE([Suitable fftw3 library found.])
fi

LIBS="${FFTW_LIB} ${LIBS}"
CFLAGS="${CFLAGS} ${FFTW_INCL} -O3 -march=native -mtune=native"
CXXFLAGS="${CXXFLAGS} ${FFTW_INCL} ${CXX11STD}"

dnl substitute C/CXXFLAGS and LIBS
AC_SUBST(CFLAGS)
AC_SUBST(LIBS)
AC_SUBST(CXXFLAGS)
AC_SUBST(OPENMP_CXXFLAGS)
dnl and do substitution in the src/Makevars.in
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

