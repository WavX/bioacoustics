dnl Process this file with autoconf to produce a configure script.

AC_INIT([bioacoustics], 0.1.2) dnl package name, version

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
      fftw_lib_path=`fftw_config --libs`
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
   FFTW_URL="ftp://ftp.fftw.org/pub/fftw/${FFTW_TGZ}"

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
   $(cd fftw-${FFTW_VERSION}; \
     ./configure --prefix=$(pwd) --enable-static --disable-shared CFLAGS=-fPIC > /dev/null; \
     make > /dev/null; \
     make install > /dev/null; )

   AC_MSG_NOTICE([Done installing fftw3 library])

else
   AC_MSG_NOTICE([Suitable fftw3 library found.])
fi


#
# libsoxr-lsr
#

dnl Select an optional path to the header file of the LSR library
AC_ARG_WITH([lsr-include],
            AC_HELP_STRING([--with-lsr-include=LSR_INCL],
                           [the location of libsoxr-lsr header file]),
            [lsr_incl_path=$withval])
if test [ -n "$lsr_incl_path" ] ; then
   LSR_INCL="-I${lsr_incl_path}"
else
  if test [ -n "${LIB_LSR}" ] ; then
     LSR_INCL="-I${LIB_LSR}/include"
  fi
fi

AC_CHECK_HEADERS(soxr-lsr.h, [lsr_ok="yes"], [lsr_ok="no"])

# dnl ditto for a library path
AC_ARG_WITH([lsr-lib],
            AC_HELP_STRING([--with-lsr-lib=LSR_LIB],
                           [the location of libsoxr-lsr library]),
            [lsr_lib_path=$withval])
if test [ -n "$lsr_lib_path" ] ; then
   LSR_LIB="-L$lsr_lib_path"
else
  if test [ -n "${LIB_LSR}" ] ; then
     LSR_LIB="-L${LIB_LSR}/lib"
  else
    if test -n "${LSR_CONFIG}"; then
      lsr_lib_path=`lsr_config --libs`
      LSR_LIB="${lsr_lib_path}"
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
   LSR_LIB="-L\"$(pwd)/soxr-${LSR_VERSION}/lib\""

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

LIBS="${FFTW_LIB} ${LSR_LIB} -lsoxr-lsr -lsoxr -lfftw3 ${LIBS}"
CPPFLAGS="${CPPFLAGS} ${FFTW_INCL} ${LSR_INCL} -O3 -march=native -mtune=native"

dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
dnl and do substitution in the src/Makevars.in
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

