AC_PREREQ([2.69])
AC_INIT([bioacoustics], [0.2.10])
AC_CONFIG_SRCDIR([src])   dnl adjust if you have a different top-level file

dnl -- basic checks
AC_PROG_CC
AC_PROG_CXX

dnl -- Find R_HOME and R executables
: ${R_HOME:=$(R RHOME)}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([could not determine R_HOME])
fi

RBIN="${R_HOME}/bin/R"
RVER="$("${R_HOME}/bin/Rscript" -e 'cat(R.version$major, ".", R.version$minor)')"

: ${CC:=$("${RBIN}" CMD config CC)}
: ${CFLAGS:=$("${RBIN}" CMD config CFLAGS)}
: ${LDFLAGS:=$("${RBIN}" CMD config LDFLAGS)}
: ${CXX:=$("${RBIN}" CMD config CXX)}
: ${CXXFLAGS:=$("${RBIN}" CMD config CXXFLAGS)}
: ${CXXPICFLAGS:=$("${RBIN}" CMD config CXXPICFLAGS)}

AC_MSG_NOTICE([CC: ${CC}])
AC_MSG_NOTICE([CXX: ${CXX}])

dnl ---------- CMake detection ----------
AC_ARG_WITH([cmake],
  AS_HELP_STRING([--with-cmake=PATH], [path to cmake executable]),
  [with_cmake_path=$withval],
  [with_cmake_path=""])

dnl If user supplied an explicit path, use it (and verify)
if test -n "${with_cmake_path}"; then
  if test -x "${with_cmake_path}"; then
    CMAKE_PROG="${with_cmake_path}"
    CMAKE_FOUND=yes
  else
    AC_MSG_WARN([--with-cmake specified but '${with_cmake_path}' is not executable. Ignoring.])
    CMAKE_PROG=""
    CMAKE_FOUND=no
  fi
else
  CMAKE_PROG=""
  CMAKE_FOUND=no

  dnl Try common names/locations. Prefer cmake then cmake3.
  AC_PATH_PROG([CMAKE_PROG], [cmake], [none])
  if test "x$CMAKE_PROG" = "xnone"; then
    AC_PATH_PROG([CMAKE_PROG], [cmake3], [none])
  fi

  if test "x$CMAKE_PROG" != "xnone"; then
    CMAKE_FOUND=yes
  else
    CMAKE_PROG=""
    CMAKE_FOUND=no
  fi
fi

if test x"${CMAKE_FOUND}" = x"yes"; then
  AC_MSG_NOTICE([Found cmake: ${CMAKE_PROG}])
else
  AC_MSG_WARN([cmake not found. cmake-dependent components will be skipped. If you need cmake, install it or run configure with --with-cmake=/path/to/cmake])
fi

AC_SUBST([CMAKE_PROG])
AC_SUBST([CMAKE_FOUND])
dnl ---------- end cmake detection ----------


dnl -- simple OpenMP feature test (portable)
AC_MSG_CHECKING([for OpenMP support])
save_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -fopenmp"
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM([], [])],
  [openmp_supported=yes],
  [openmp_supported=no]
)
CFLAGS="$save_CFLAGS"
if test "$openmp_supported" = "yes"; then
  OPENMP_CFLAGS="-fopenmp"
  OPENMP_LIBS="-fopenmp"
  AC_MSG_RESULT([yes])
else
  OPENMP_CFLAGS=""
  OPENMP_LIBS=""
  AC_MSG_RESULT([no])
fi

if test x"${OPENMP_CFLAGS}" != "x"; then
  FFTW_OPENMP="--enable-openmp"
else
  FFTW_OPENMP=""
fi

dnl -- Ensure MAKE has a sensible default
if test -z "${MAKE}"; then
  MAKE=make
fi

#
# fftw selection
#
AC_ARG_WITH([fftw-include],
            AS_HELP_STRING([--with-fftw-include=DIR], [location of fftw3 header file]),
            [fftw_incl_path=$withval],
            [fftw_incl_path=""])

if test -n "${fftw_incl_path}"; then
  FFTW_INCL="-I${fftw_incl_path}"
elif test -n "${LIB_FFTW}"; then
  FFTW_INCL="-I${LIB_FFTW}/include"
else
  FFTW_INCL=""
fi

AC_CHECK_HEADERS([fftw3.h],
                 [fftw_ok="yes"],
                 [fftw_ok="no"])

AC_ARG_WITH([fftw-lib],
            AS_HELP_STRING([--with-fftw-lib=DIR], [location of fftw3 library]),
            [fftw_lib_path=$withval],
            [fftw_lib_path=""])

if test -n "${fftw_lib_path}"; then
  FFTW_LIB="-L${fftw_lib_path}"
elif test -n "${LIB_FFTW}"; then
  FFTW_LIB="-L${LIB_FFTW}/lib"
else
  FFTW_LIB=""
fi

AC_SEARCH_LIBS([fftw_plan_r2r_1d], [fftw3],
               [fftw_ok="yes"],
               [fftw_ok="no"])

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

  AC_MSG_NOTICE([Downloading ${FFTW_URL}])
  "${R_HOME}/bin/Rscript" --vanilla -e "download.file(url='${FFTW_URL}', destfile='${FFTW_TGZ}')"

  AC_MSG_NOTICE([Extracting ${FFTW_TGZ}])
  tar xzf "${FFTW_TGZ}"

  rm -f "${FFTW_TGZ}"

  AC_MSG_NOTICE([Starting to install fftw3 library to $(pwd)/fftw-${FFTW_VERSION}])
  ( set -e
    cd "fftw-${FFTW_VERSION}"
    ./configure --prefix="$(pwd)" --enable-static --disable-shared ${FFTW_OPENMP} CFLAGS=-fPIC > /dev/null
    ${MAKE} > /dev/null
    ${MAKE} install > /dev/null
  )
  FFTW3_INSTALL=$?

  if test x"${FFTW3_INSTALL}" != x"0" ; then
    AC_MSG_ERROR([Could not install fftw3 library, try to install it manually?])
  fi

  case " $LIBS " in
    *" lfftw3"*)
      ;;
    *)
      LIBS="$LIBS -lfftw3"
      ;;
  esac

  AC_MSG_NOTICE([Done installing fftw3 library])
else
  AC_MSG_NOTICE([Suitable fftw3 library found.])
fi

dnl -- Compose final flags and libs
LIBS="${FFTW_LIB} ${LIBS} ${OPENMP_LIBS}"
CFLAGS="${CFLAGS} ${FFTW_INCL} -O3 -march=native -mtune=native ${OPENMP_CFLAGS}"
CXXFLAGS="${CXXFLAGS} ${FFTW_INCL} ${CXXPICFLAGS}"

dnl substitute variables for Makevars.in
AC_SUBST([CFLAGS])
AC_SUBST([LIBS])
AC_SUBST([CXXFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
