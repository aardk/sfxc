dnl Process this file with autoconf to produce a configure script.
AC_INIT(SFXC, 3.5.0)
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR(m4)
AM_INIT_AUTOMAKE()

dnl backup of generic flags
CXXFLAGS_IN=$CXXFLAGS
CFLAGS_IN=$CFLAGS

dnl Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AX_PROG_CXX_MPI
AC_PROG_F77
AC_PROG_INSTALL
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

AC_F77_LIBRARY_LDFLAGS

dnl Checking generic c libraries
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(pthread, sem_init)
dnl Needed for the optimized channel extractor
AC_CHECK_LIB(dl, dlopen)

dnl setting flags for sfxc
SFXC_CXXFLAGS='-I${top_srcdir}/include -std=gnu++03'
SFXC_LDADD=' '
UTILS_CXXFLAGS='-std=gnu++03'

AC_CHECK_FUNCS(sincos)

dnl backup of generic flags
LDFLAGS_IN=$LDFLAGS
LIBS_IN=$LIBS

dnl ---------------------------------------------------------------------------
dnl Checking for ./configure options
dnl ---------------------------------------------------------------------------
dnl debug
AC_ARG_ENABLE(debug,
              [  --enable-debug           Enable debugging information [default=no]],
              USE_DEBUG="$enableval", USE_DEBUG="no")

AC_ARG_ENABLE(optimization,
              [  --disable-optimization   Disable optimization [default=yes]],
              USE_OPTIMIZE="$enableval", USE_OPTIMIZE="yes")

AC_ARG_ENABLE(mpi,
              [  --enable-mpi             Compile with dependency to MPI [default=yes]],
              USE_MPI="$enableval", USE_MPI="yes")

AC_ARG_ENABLE(progress,
              [  --enable-progress        Enable progress messages [default=no]],
              PRINT_PROGRESS="$enableval", PRINT_PROGRESS="no")

AC_ARG_ENABLE(ipp,
              [  --enable-ipp             Use the Intel Performance Primitives [default=no]],
              USE_IPP="$enableval", USE_IPP="no")

AC_ARG_ENABLE(armpl,
              [  --enable-armpl           Use the Arm Performance Libraries [default=no]],
              USE_ARMPL="$enableval", USE_ARMPL="no")

AC_ARG_ENABLE(gsl,
              [  --enable-gsl             Compile with dependency to GSL [default=yes]],
              USE_GSL="$enableval", USE_GSL="yes")

AC_ARG_ENABLE(double,
              [  --enable-double          Use double precision [default=no]],
              USE_DOUBLE="$enableval", USE_DOUBLE="no")

AC_ARG_ENABLE(gendelay,
              [  --enable-gendelay        Compile with the generate_delay_model (require Fortan) [default=yes]],
              WITH_GENDELAY="$enableval", WITH_GENDELAY="yes")

AC_ARG_ENABLE(calc10,
              [  --enable-calc10          Compile with the cacl10 (require Fortan) [default=yes]],
              WITH_CALC="$enableval", WITH_CALC="yes")

AC_ARG_WITH(ipp-path, 
              [  --with-ipp-path=PATH    Location of the ipp libraries [default=$IPPROOT]], 
              IPP_PATH="$withval", IPP_PATH=$IPPROOT)

AC_ARG_WITH(armpl-path, 
              [  --with-armpl-path=PATH  Location of the armpl libraries [default=$ARMPL_DIR]], 
              ARMPL_PATH="$withval", ARMPL_PATH=$ARMPL_DIR)

dnl ---------------------------------------------------------------------------
dnl Checking external libraries and dependency
dnl ---------------------------------------------------------------------------
if test $USE_DEBUG = "yes"; then
  SFXC_CFLAGS="$SFXC_CFLAGS -g -DSFXC_PRINT_DEBUG -rdynamic"
  SFXC_CXXFLAGS="$SFXC_CXXFLAGS -g -DSFXC_PRINT_DEBUG -rdynamic"
  UTILS_CXXFLAGS="$UTILS_CXXFLAGS -g"
else
  SFXC_CFLAGS="$SFXC_CFLAGS -rdynamic"
  SFXC_CXXFLAGS="$SFXC_CXXFLAGS -rdynamic"
fi

if test $USE_OPTIMIZE = "yes"; then
  OPTIMIZE_CFLAGS="-O2 -ffast-math"
  OPTIMIZE_CXXFLAGS="-O2 -ffast-math"
fi

if test $PRINT_PROGRESS = "yes"; then
  SFXC_CFLAGS="$SFXC_CFLAGS -DPRINT_PROGRESS -DSFXC_PRINT_DEBUG"
  SFXC_CXXFLAGS="$SFXC_CXXFLAGS -DPRINT_PROGRESS -DSFXC_PRINT_DEBUG"
fi

if test $USE_MPI = "yes"; then
  SFXC_CFLAGS="$SFXC_CFLAGS -DUSE_MPI"
  SFXC_CXXFLAGS="$SFXC_CXXFLAGS -DUSE_MPI"
fi

if test $USE_DOUBLE = "yes" ; then
  AC_DEFINE([USE_DOUBLE], [1], [Use double precision])
fi

if test $USE_IPP = "yes"; then
  USE_FFTW3="no"
  AC_DEFINE([USE_IPP], [1], [Use the Intel Performance Primitives])
  AC_CHECK_LIB(ippsem64t, ippsGetLibVersion, [OLD_IPP=1], [OLD_IPP=0], [-L$IPP_PATH/sharedlib -Wl,-no-as-needed])
  if test $OLD_IPP == 1; then
    if test "$IPP_PATH" != ""; then
       FFT_LDADD="-L$IPP_PATH/sharedlib"
       FFT_CPPFLAGS="-I$IPP_PATH/include"
       LDFLAGS="$LDFLAGS -L$IPP_PATH/sharedlib"
       CPPFLAGS="$CPPFLAGS -I$IPP_PATH/include"
    fi
    FFT_LDADD="$FFT_LDADD -lippsem64t -lippvmem64t -lippcoreem64t"
    FFT_LDADD_FLOAT="$FFT_LDADD"
  else
    AC_CHECK_LIB(ipps, ippsGetLibVersion, [], [
        echo "Could not find Intel IPP library."
        exit -1
        ], -L$IPP_PATH/lib/intel64 -L$IPP_PATH/../compiler/lib/intel64 -liomp5)
    if test "$IPP_PATH" != ""; then
       FFT_LDADD="-L$IPP_PATH/lib/intel64 -L$IPP_PATH/../compiler/lib/intel64/"
       FFT_CPPFLAGS="-I$IPP_PATH/include"
       LDFLAGS="$LDFLAGS -L$IPP_PATH/lib/intel64 -L$IPP_PATH/../compiler/lib/intel64/"
       CPPFLAGS="$CPPFLAGS -I$IPP_PATH/include"
    fi
    FFT_LDADD="$FFT_LDADD -lipps -lippvm -lippcore -liomp5"
    FFT_LDADD_FLOAT="$FFT_LDADD"
  fi
elif test $USE_ARMPL = "yes"; then
  USE_FFTW3="no"
  if test "$ARMPL_PATH" != ""; then
    FFT_LDADD="-L$ARMPL_PATH/lib -Wl,-rpath-link,$ARMPL_PATH/lib"
    FFT_CPPFLAGS="-I$ARMPL_PATH/include"
  fi
  FFT_LDADD_FLOAT="$FFT_LDADD -larmpl -lamath"
  FFT_LDADD="$FFT_LDADD -larmpl -lamath"
else
  USE_FFTW3="yes"
  AC_CHECK_LIB(fftw3f, fftwf_execute, [], [
        echo "Couldn't find fftw3f library."
        exit -1
        ])
  FFT_LDADD_FLOAT="-lfftw3f"
  if test $USE_DOUBLE = "yes"; then
    AC_CHECK_LIB(fftw3, fftw_execute, [], [
        echo "Couldn't find fftw3 library (double precision). "
        exit -1
        ])
    FFT_LDADD="-lfftw3"
  else
    FFT_LDADD="-lfftw3f"
  fi
fi
SFXC_CPPFLAGS="$SFXC_CPPGLAGS $FFT_CPPFLAGS"
SFXC_LDADD="$SFXC_LDADD $FFT_LDADD"

dnl ---------------------------------------------------------------------------
dnl setting shipped libraries
dnl ---------------------------------------------------------------------------
dnl setting Damien's libraries
DAMIEN_CXXFLAGS='-I${top_srcdir}/lib/common/src -I${top_srcdir}/lib/containers/src -I${top_srcdir}/lib/netlib/src -I${top_srcdir}/lib/testunit/src '
DAMIEN_LDADD='-L${top_srcdir}/lib/testunit -ltestunit -L${top_srcdir}/lib/netlib -lnetlib -L${top_srcdir}/lib/containers -lcontainers -L${top_srcdir}/lib/common -lcommon -ltestunit'

SFXC_CXXFLAGS="$SFXC_CXXFLAGS $DAMIEN_CXXFLAGS"
SFXC_LDADD="$SFXC_LDADD $DAMIEN_LDADD"

dnl setting json-cpp
JSONCPP_CXXFLAGS=' -I${top_srcdir}/lib/jsoncpp/include '
JSONCPP_LDADD=' -L${top_srcdir}/lib/jsoncpp -ljsoncpp '

SFXC_CXXFLAGS="$SFXC_CXXFLAGS $JSONCPP_CXXFLAGS"
SFXC_LDADD="$SFXC_LDADD $JSONCPP_LDADD"

dnl if test $USE_MPI = "yes"; then
dnl   AC_CHECK_LIB(mpi, MPI_Init, [], [
dnl                echo "Couldn't find MPI library. "
dnl                exit -1
dnl               ])
dnl fi

dnl setting gsl
dnl first check gslcblas
if test $USE_GSL = "yes"; then
  AC_CHECK_LIB([gslcblas],[cblas_dgemm])
  AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
  dnl found gsl, use the global version
  GSL_CXXFLAGS=''
  GSL_LDADD='-lgsl -lgslcblas'
fi

WITH_SFXC="no"
WITH_SFXC_UTILS="no"
if test $USE_GSL = "yes"; then
  WITH_SFXC_UTILS="yes"
  if test $USE_MPI = "yes"; then
    WITH_SFXC="yes"
  fi
fi

if test $USE_GSL = "no"; then
      WITH_GENDELAY="no"
fi


SFXC_CXXFLAGS="$SFXC_CXXFLAGS $GSL_CXXFLAGS"
SFXC_LDADD="$SFXC_LDADD $GSL_LDADD"

dnl setting vex_parser
VEX_PARSER_CXXFLAGS='-I${top_srcdir}/lib/vex_parser/install/include '
VEX_PARSER_LDADD='-L${top_srcdir}/lib/vex_parser/install/lib -lvex_parser '
AC_SUBST(VEX_PARSER_CXXFLAGS)
AC_SUBST(VEX_PARSER_LDADD)

SFXC_CXXFLAGS="$SFXC_CXXFLAGS $VEX_PARSER_CXXFLAGS $OPTIMIZE_CXXFLAGS"
SFXC_LDADD="$SFXC_LDADD $VEX_PARSER_LDADD $OPTIMIZE_LDADD"

BASE_CXXFLAGS='-I${top_srcdir}/include '
BASE_LDADD=' '

UTILS_CPPFLAGS="$UTILS_CPPFLAGS $FFT_CPPFLAGS"
UTILS_CXXFLAGS="$UTILS_CXXFLAGS $BASE_CXXFLAGS $VEX_PARSER_CXXFLAGS $JSONCPP_CXXFLAGS $DAMIEN_CXXFLAGS $GSL_CXXFLAGS $OPTIMIZE_CXXFLAGS"
UTILS_LDADD="$BASE_LDADD $VEX_PARSER_LDADD $JSONCPP_LDADD $DAMIEN_LDADD $FFT_LDADD_FLOAT $GSL_LDADD $OPTIMIZE_LDADD"

dnl resetting generic flags
CXXFLAGS=$CXXFLAGS_IN
CFLAGS=$CFLAGS_IN
LDFLAGS=$LDFLAGS_IN
LIBS=$LIBS_IN

dnl export these variable (so Makefile substitutions
dnl can be made).
dnl AC_SUBST(FFT_CPPFLAGS)
dnl AC_SUBST(FFT_LDADD)
dnl AC_SUBST(FFT_LDADD_FLOAT)
dnl AC_SUBST(SFXCFFT_CPPFLAGS)
dnl AC_SUBST(SFXCFFT_LDADD)
dnl AC_SUBST(SFXCFFT_CPPFLAGS_FLOAT)
dnl AC_SUBST(SFXCFFT_LDADD_FLOAT)
AC_SUBST(DAMIEN_CXXFLAGS)
AC_SUBST(DAMIEN_LDADD)
AC_SUBST(JSONCPP_CXXFLAGS)
AC_SUBST(JSONCPP_LDADD)

AC_SUBST(UTILS_CPPFLAGS)
AC_SUBST(UTILS_CXXFLAGS)
AC_SUBST(UTILS_LDADD)
AC_SUBST(BASE_CXXFLAGS)
AC_SUBST(BASE_LDADD)
AC_SUBST(GSL_CXXFLAGS)
AC_SUBST(GSL_LDADD)

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LDADD)

AC_SUBST(SFXC_CFLAGS)
AC_SUBST(SFXC_LDFLAGS)
AC_SUBST(SFXC_CPPFLAGS)
AC_SUBST(SFXC_CXXFLAGS)
AC_SUBST(SFXC_LDADD)

AC_SUBST(WITH_CALC)
AC_SUBST(WITH_GENDELAY)
AC_SUBST(WITH_SFXC)
AC_SUBST(WITH_SFXC_UTILS)

dnl echo "SFXC_CFLAGS: $SFXC_CFLAGS"
dnl echo "SFXC_CXXFLAGS: $SFXC_CXXFLAGS"
dnl echo "SFXC_LDFLAGS: $SFXC_LDFLAGS"
dnl echo "SFXC_LDADD: $SFXC_LDADD"


AM_CONDITIONAL(DOUBLE_PRECISION, [test "$USE_DOUBLE" = "yes"])
AM_CONDITIONAL(GENDELAY, [test "$WITH_GENDELAY" = "yes"])
AM_CONDITIONAL(CALCTEN, [test "$WITH_CALC" = "yes"])
AM_CONDITIONAL(SFXC_UTILS, [test "$WITH_SFXC_UTILS" = "yes"])
AM_CONDITIONAL(SFXC, [test "$WITH_SFXC" = "yes"])

# Fix compilation of calc10 with gfortran 10+
if test "$WITH_CALC" == "yes"; then
  AC_LANG_PUSH([Fortran 77])
  AX_CHECK_COMPILE_FLAG([-fallow-argument-mismatch], [gfortran_version_10plus=yes], [gfortran_version_10plus=no])
  if test "$gfortran_version_10plus" == "yes"; then
    FFLAGS="$FFLAGS -fallow-argument-mismatch"
  fi
  AC_LANG_POP([Fortran 77])
fi

AC_OUTPUT(Makefile
          lib/Makefile
          lib/calc10/Makefile
          lib/calc10/include/Makefile
          lib/calc10/lib/Makefile
          lib/calc10/src/Makefile
          lib/common/Makefile
          lib/common/test/Makefile
          lib/containers/Makefile
          lib/containers/test/Makefile
          lib/netlib/Makefile
          lib/testunit/Makefile
          lib/jsoncpp/Makefile
          lib/vex_parser/Makefile
          src/Makefile
          include/Makefile
          utils/Makefile
          utils/delay/Makefile
          doc/Makefile
          doc/html/Makefile
         )

echo "=================== compilation summary =================="
echo "Use mpi: $USE_MPI"
echo "Use fftw3: $USE_FFTW3"
echo "Use ipp: $USE_IPP"
echo "Use armpl: $USE_ARMPL"
echo "Use gsl: $USE_GSL"
echo "Compile sfxc with: $SFXC_CXXFLAGS $CXXFLAGS"
echo "Compile utils with: $UTILS_CXXFLAGS $CXXFLAGS"
echo "Build sfxc: $WITH_SFXC"
echo "Build utils/*: $WITH_SFXC_UTILS"
echo "Build generate_delay_model: $WITH_GENDELAY"
echo "Build calc10: $WITH_CALC"
echo "=========================================================="

