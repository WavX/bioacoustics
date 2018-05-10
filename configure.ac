dnl Process this file with autoconf to produce a configure script.

AC_INIT([bioacoustics], 0.1.0) dnl package name, version

# require at least autoconf 2.61
AC_PREREQ(2.61)

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

RBIN="${R_HOME}/bin/R"

: ${CC=`"${RBIN}" CMD config CC`}
: ${CXX=`"${RBIN}" CMD config CXX`}
: ${CPP=`"${RBIN}" CMD config CPP`}
: ${CFLAGS=`"${RBIN}" CMD config CFLAGS`}
: ${CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`}
: ${LDFLAGS=`"${RBIN}" CMD config LDFLAGS`}
AC_PROG_CC
AC_PROG_CPP

AC_MSG_NOTICE([CC: ${CC}])
AC_MSG_NOTICE([CXX: ${CXX}])


#
# fftw
#

dnl Select an optional include path, from a configure option
dnl or from an environment variable.
AC_ARG_WITH([fftw-include],
            AC_HELP_STRING([--with-fftw-include=INCLUDE_PATH],
                           [the location of fftw3 header file]),
            [fftw_include_path=$withval])
FFTW_CPPFLAGS="-I."
if test [ -n "$fftw_include_path" ] ; then
   FFTW_CPPFLAGS="-I. -I${fftw_include_path}"
else
  if test [ -n "${FFTW_INCLUDE}" ] ; then
     FFTW_CPPFLAGS="-I. -I${FFTW_INCLUDE}"
  fi
fi

AC_CHECK_HEADERS(fftw3.h,,
    AC_MSG_ERROR([fftw3.h not found in standard or given locations]))


# dnl ditto for a library path
AC_ARG_WITH([fftw-lib],
            AC_HELP_STRING([--with-fftw-lib=LIB_PATH],
                           [the location of fftw library]),
            [fftw_lib_path=$withval])
if test [ -n "$fftw_lib_path" ] ; then
   FFTW_LIBS="-L$fftw_lib_path"
else
  if test [ -n "${FFTW_LIBS}" ] ; then
     FFTW_LIBS="-L${FFTW_LIBS}"
  else
    if test -n "${FFTW_CONFIG}"; then
      fftw_lib_path=`fftw_config --libs | sed s/-lfftw//`
      FFTW_LIBS="${fftw_lib_path}"
    fi
  fi
fi

AC_SEARCH_LIBS(fftw_plan_r2r_1d, fftw3, ,
    AC_MSG_ERROR("Library fftw3 not found in standard or given locations"))


#
# libsoxr-lsr
#

dnl Select an optional include path, from a configure option
dnl or from an environment variable.
AC_ARG_WITH([libsoxr-lsr-include],
            AC_HELP_STRING([--with-libsoxr-lsr-include=INCLUDE_PATH],
                           [the location of libsoxr-lsr header file]),
            [lsr_include_path=$withval])
LSR_CPPFLAGS="-I."
if test [ -n "$lsr_include_path" ] ; then
   LSR_CPPFLAGS="-I. -I${lsr_include_path}"
else
  if test [ -n "${LSR_INCLUDE}" ] ; then
     LSR_CPPFLAGS="-I. -I${LSR_INCLUDE}"
  fi
fi

AC_CHECK_HEADERS(soxr-lsr.h, [lsr_ok="yes"], [lsr_ok="no"])

# dnl ditto for a library path
AC_ARG_WITH([libsoxr-lsr-lib],
            AC_HELP_STRING([--with-libsoxr-lsr-lib=LIB_PATH],
                           [the location of libsoxr-lsr library]),
            [lsr_lib_path=$withval])
if test [ -n "$lsr_lib_path" ] ; then
   LSR_LIBS="-L$lsr_lib_path"
else
  if test [ -n "${LSR_LIBS}" ] ; then
     LSR_LIBS="-L${LSR_LIBS}"
  else
    if test -n "${LSR_CONFIG}"; then
      lsr_lib_path=`lsr_config --libs | sed s/-lsoxr-lsr//`
      LSR_LIBS="${lsr_lib_path}"
    fi
  fi
fi

AC_SEARCH_LIBS(src_simple, soxr-lsr, [lsr_ok="yes"], [lsr_ok="no"])

## If libsoxr-lsr is not found, download libsoxr and build it locally
if test x"${lsr_ok}" = x"no"; then
   AC_MSG_NOTICE([Need to download and build libsoxr])

   ## define LSR version
   LSR_VERSION=0.1.3

   ## define LSR file and download URL
   LSR_TGZ="${LSR_VERSION}.tar.gz"
   LSR_URL="https://github.com/WavX/soxr/archive/${LSR_TGZ}"

   ## C Compiler options
   LSR_CFLAGS=

   ## additional C Compiler options for linking
   LSR_CLINKFLAGS=

   ## Libraries necessary to link with the package
   LSR_LIBS="-L\"$(pwd)/soxr-${LSR_VERSION}/lib\" -lsoxr -lsoxr-lsr"
   LIBS="${LSR_LIBS} ${LIBS}"

   ## Necessary Include dirs
   LSR_INCL="-I$(pwd)/soxr-${LSR_VERSION}/include"

   ## Download libsoxr source code
   $("${R_HOME}/bin/Rscript" --vanilla -e "download.file(url='${LSR_URL}', destfile='${LSR_TGZ}')")

   ## Extract libsoxr source code and remove .tar.gz
   $("${R_HOME}/bin/Rscript" --vanilla -e "untar(tarfile='${LSR_TGZ}', exdir='$(pwd)')")
   $(rm -rf ${LSR_TGZ})

   ## Compile libsoxr source code and clean up
   AC_MSG_NOTICE([Starting to install soxr-lsr library to $(pwd)/soxr-${LSR_VERSION}])
   $(cd soxr-${LSR_VERSION}; \
     cmake $(pwd) -Wno-dev -DCMAKE_BUILD_TYPE=Release \
                  -DCMAKE_INSTALL_PREFIX=$(pwd) \
                  -DBUILD_SHARED_LIBS=OFF \
                  -DCMAKE_POSITION_INDEPENDENT_CODE=ON > /dev/null; \
     make > /dev/null; \
     make install > /dev/null; )

   AC_MSG_NOTICE([Done installing soxr-lsr library])

   LSR_CPPFLAGS="${LSR_INCL}"
else
   AC_MSG_NOTICE([Suitable libsoxr-lsr library found.])
fi

LIBS="${LIBS} ${FFTW_LIBS} ${LSR_LIBS}"
CPPFLAGS="${CPPFLAGS} ${FFTW_CPPFLAGS} ${LSR_CPPFLAGS} -O3 -march=native -mtune=native"

dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
dnl and do substitution in the src/Makevars.in
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

