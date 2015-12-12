#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(mfix, 2015.2, develop@mfix.netl.doe.gov)

m4_include([build-aux/ax_check_compile_flag.m4])
m4_include([build-aux/ax_prog_fc_mpi.m4])

AC_MSG_NOTICE([==============================================================])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([       MFIX: Multiphase Flow with Interphase eXchanges])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                      mfix.netl.doe.gov])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                        Version ${MFIX_VERSION}])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([==============================================================])

AC_LANG(Fortran)

AC_ARG_ENABLE([dmp],
[  --enable-dmp    Turn on MPI support],
[case "${enableval}" in
  yes) dmp=true ;;
  no)  dmp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-dmp]) ;;
esac],[dmp=false])

AX_PROG_FC_MPI([test x"$dmp" = xtrue],[use_mpi=yes],[
       use_mpi=no
       if test x"$dmp" = xtrue; then
         AC_MSG_FAILURE([DMP support requested, but no MPI compiler found.  Specify an MPI compiler as a configure_mfix option, such as FC=mpifort or FC=mpif90])
       fi
     ])

# Checks for programs.
AC_PROG_GREP
AC_PROG_LN_S
AC_PROG_MKDIR_P
AC_PROG_RANLIB
AC_PROG_FC
AC_PROG_F77([$FC])
AC_PROG_CC
AC_PROG_SED
AC_OPENMP()
AX_CHECK_COMPILE_FLAG([-assume byterecl], [
                        FCFLAGS="$FCFLAGS -assume byterecl"])

AC_CONFIG_SRCDIR([model/mfix.f])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Wno-portability no-dependencies])

AC_ARG_ENABLE([mkl],
[  --enable-mkl    Turn on Math Kernel Library support],
[case "${enableval}" in
  yes) mkl=true ;;
  no)  mkl=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-mkl]) ;;
esac],[mkl=false])
AM_CONDITIONAL([MKL], [test x$mkl = xtrue])

AC_ARG_VAR([NETCDF_LIB],[ Path to NETCDF module ])
AC_ARG_VAR([NETCDF_INCLUDE],[ Path to NETCDF include file ])
AM_CONDITIONAL([USE_NETCDF], [test ! -z "$NETCDF_INCLUDE"])

AC_ARG_VAR([TAU_MAKEFILE],[ Path to TAU Makefile ])
AC_ARG_VAR([TAU_F90],[ Actual Fortran compiler used by TAU compiler wrapper ])

AC_ARG_ENABLE([smp],
[  --enable-smp    Turn on OpenMP support],
[case "${enableval}" in
  yes) smp=true ;;
  no)  smp=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-smp]) ;;
esac],[smp=false])
AM_CONDITIONAL([OPENMP], [test x$smp = xtrue])

AM_CONDITIONAL([MPI], [test x$dmp = xtrue])

AC_FC_PP_DEFINE
AC_FC_PP_SRCEXT(f)
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH([unlimited])
AC_FC_IMPLICIT_NONE

AC_FC_MODULE_FLAG
AC_FC_MODULE_OUTPUT_FLAG

AC_CONFIG_FILES([Makefile build-aux/Makefile.usr])

# Define Autoconf variables so Automake doesn't see GNU Make conditionals
make_ifnotclean='
ifneq ($(MAKECMDGOALS),clean)
'
AC_SUBST([make_ifnotclean])
AM_SUBST_NOTMAKE([make_ifnotclean])

make_endif='
endif
'
AC_SUBST([make_endif])
AM_SUBST_NOTMAKE([make_endif])
AC_CONFIG_COMMANDS([cases],
[cp $ac_srcdir/build-aux/test.sh $ac_builddir/build-aux; for casedir in benchmarks legacy_tests tests tutorials; do [[ "$ac_builddir" != "$ac_srcdir" ]] && cp -r $ac_srcdir/$casedir $ac_builddir; chmod --recursive +w $ac_builddir/$casedir; done ],
)

AC_OUTPUT

AC_MSG_NOTICE([==============================================================])
AC_MSG_NOTICE([Now run "make" to build mfix with:])
AC_MSG_NOTICE([FC=$FC])
AC_MSG_NOTICE([FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([F77=$F77])
AC_MSG_NOTICE([FFLAGS=$FFLAGS])
AC_MSG_NOTICE([==============================================================])
