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
   LIBS="-L$fftw_lib_path ${LIBS}"
else
  if test [ -n "${FFTW_LIBS}" ] ; then
     LIBS="-L${FFTW_LIBS} ${LIBS}"
  else
    if test -n "${FFTW_CONFIG}"; then
      fftw_lib_path=`fftw_config --libs | sed s/-lfftw//`
      LIBS="${fftw_lib_path} ${LIBS}"
    fi
  fi
fi

AC_SEARCH_LIBS(fftw_plan_r2r_1d, fftw3, ,
    AC_MSG_ERROR("Library fftw3 not found in standard or given locations"))


if test -n "${FFTW_CONFIG}"; then
  FFTW_CPPFLAGS=`fftw_config --cflags`
fi

#
# libsamplerate
#

dnl Select an optional include path, from a configure option
dnl or from an environment variable.
AC_ARG_WITH([libsamplerate-include],
            AC_HELP_STRING([--with-libsamplerate-include=INCLUDE_PATH],
                           [the location of libsamplerate header file]),
            [lsr_include_path=$withval])
LSR_CPPFLAGS="-I."
if test [ -n "$lsr_include_path" ] ; then
   LSR_CPPFLAGS="-I. -I${lsr_include_path}"
else
  if test [ -n "${LSR_INCLUDE}" ] ; then
     LSR_CPPFLAGS="-I. -I${LSR_INCLUDE}"
  fi
fi

AC_CHECK_HEADERS(samplerate.h,,
   AC_MSG_ERROR([samplerate.h not found in standard or given locations]))

# dnl ditto for a library path
AC_ARG_WITH([libsamplerate-lib],
            AC_HELP_STRING([--with-libsamplerate-lib=LIB_PATH],
                           [the location of libsamplerate library]),
            [lsr_lib_path=$withval])
if test [ -n "$lsr_lib_path" ] ; then
   LIBS="-L$lsr_lib_path ${LIBS}"
else
  if test [ -n "${LSR_LIBS}" ] ; then
     LIBS="-L${LSR_LIBS} ${LIBS}"
  else
    if test -n "${LSR_CONFIG}"; then
      lsr_lib_path=`lsr_config --libs | sed s/-lsamplerate//`
      LIBS="${lsr_lib_path} ${LIBS}"
    fi
  fi
fi

AC_SEARCH_LIBS(src_simple, samplerate, ,
    AC_MSG_ERROR("Library libsamplerate not found in standard or given locations"))

if test -n "${LSR_CONFIG}"; then
  LSR_CPPFLAGS=`lsr_config --cflags`
fi
CPPFLAGS="${CPPFLAGS} ${FFTW_CPPFLAGS} ${LSR_CPPFLAGS} -O3 -march=native -mtune=native"


dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
dnl and do substitution in the src/Makevars.in
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
