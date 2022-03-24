AC_INIT([CaTSM], [1.0], [jmbrown2@nps.edu])
AM_INIT_AUTOMAKE
AC_LANG([Fortran])
AC_FC_SRCEXT([F90])

AC_ARG_WITH([mpi],
AS_HELP_STRING([--with-mpi],[build with mpi (default)])
AS_HELP_STRING([--without-mpi],[build without mpi])
,[
    with_mpi=$withval
    with_pnetcdf=$withval
],[
    with_mpi="yes"
    with_pnetcdf="yes"
])

my_save_fcflags=$FCFLAGS
FCFLAGS="$FCFLAGS -cpp"
AS_IF([test "x$FC" != "x"],[
    find_compiler=no
])

AS_IF([test "x$find_compiler" = "xno"], [
    AC_MSG_CHECKING([for mpif.h])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[#include "mpif.h"])],
    [AC_MSG_RESULT(yes)], 
    [AC_MSG_RESULT(no)
    with_mpi=no])
],[
    AS_IF([test "x$with_mpi" != "xno"], [
        AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
        AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
        acx_mpi_save_FC="$FC"
        FC="$MPIFC"
        AC_SUBST(MPIFC)
    ],[
        AC_PROG_FC
    ])
])

AS_IF([test "x$FFTW_HOME" != "x"], [
    AC_MSG_NOTICE([FFTW_HOME set, updating paths])
    FCFLAGS="-I$FFTW_HOME/include $FCFLAGS"
    LDFLAGS="-L$FFTW_HOME/lib $LDFLAGS"
])

AS_IF([test "x$PNETCDF_HOME" != "x"], [
    AC_MSG_NOTICE([PNETCDF_HOME set, updating paths])
    FCFLAGS="-I$PNETCDF_HOME/include -cpp -C $FCFLAGS"
    LDFLAGS="-L$PNETCDF_HOME/lib $LDFLAGS"
    with_pnetcdf="yes"
])

AS_IF([test "x$NETCDF_HOME" != "x"], [
    AC_MSG_NOTICE([NETCDF_HOME set, updating paths])
    FCFLAGS="-I$NETCDF_HOME/include $FCFLAGS"
    LDFLAGS="-L$NETCDF_HOME/lib $LDFLAGS"
])

AS_IF([test "x$with_mpi" = "xyes"], [
    with_pnetcdf="yes"
],[
    AS_IF([test "x$with_pnetcdf" = "xyes"], [
    AC_MSG_ERROR(
        [Cannot run pnetcdf without MPI]
    )])
])

AS_IF([test "x$with_mpi" != "xyes"], [
    AS_IF([test "x$with_pnetcdf" = "xyes"], [
    AC_MSG_ERROR(
        [Cannot run netcdf with MPI]
    )])
])

AS_IF([test "x$with_pnetcdf" = "xyes"], [
    AC_CHECK_LIB([pnetcdf], [nfmpi_close],
    [AC_SUBST([LIBNC], ["-lpnetcdf"])
    AC_DEFINE([PNETCDF],[])],
    [AC_MSG_ERROR(
        [library pnetcdf test failed])
    with_mpi=no],
    [])
],[
    AC_CHECK_LIB([netcdff], [nf_close],[
        AC_SUBST([LIBNC], ["-lnetcdff"])
    ],[
        AC_MSG_ERROR([netcdf test failed])
    ])
    AC_DEFINE([NETCDF],[])
])

AS_IF([test "x$with_mpi" = "xyes"], [
    AC_DEFINE([MPI],[])
])

my_save_fcflags="$FCFLAGS"
FCFLAGS="$FCFLAGS -fallow-argument-mismatch"
AC_MSG_CHECKING([whether FC supports -fallow-argument-mismatch])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="-fallow-argument-mismatch"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS="$my_save_fcflags $AM_FCFLAGS -O3"
AC_MSG_CHECKING([whether FC supports -O3])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="$AM_FCFLAGS -O3"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS="$my_save_fcflags $AM_FCFLAGS -g"
AC_MSG_CHECKING([whether FC supports -g])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="$AM_FCFLAGS -g"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS=$my_save_fcflags
AC_SUBST([AM_FCFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CHECK_LIB([fftw3], [dfftw_plan_many_dft],
    [AC_SUBST([LIBFFTW3], ["-lfftw3"])],
    [AC_MSG_ERROR(
        [fftw3 test failed (try specifying "LDFLAGS='-Lfftw3/path/lib'")])],
    [])

AM_CONDITIONAL([MPI], [test x$with_mpi = "xyes"])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT