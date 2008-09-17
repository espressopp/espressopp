#
# SYNOPSIS
#
#   AXES_BOOST_MPI
#
# DESCRIPTION
#
#   Test for the mpi library from the Boost C++ libraries.
#   The macro requires a preceding call to AXES_BOOST_BASE,
#   AXES_BOOST_SERIALIZATION and AXES_MPI.
#   After the macro has been executed, axes_cv_boost_mpi denotes,
#   whether boost mpi has been found. If it was found,
#   axes_cv_boost_lib_suffix contains the suffix of the boost
#   libraries, and axes_cv_boost_lib_path contains the path where the
#   libraries can be found.
#   axes_cv_boost_mpi_links denotes, whether the specified parameters
#   can be used to link a boost mpi program. If is "no", this usually
#   denotes that boost mpi was compiled with an mpi compiler that is
#   not compatible with the compiler specified here.
#   
#
#   This macro calls:
#
#     AC_SUBST(BOOST_MPI_LIBS)
#
#   And defines:
#
#     HAVE_BOOST_MPI
#
#   And sets the ac variables:
#
#     axes_cv_boost_mpi
#     axes_cv_boost_mpi_links
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
# LAST MODIFICATION
#
#   2008-09-17
#
# COPYLEFT
#
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#
#   based on AX_BOOST_SERIALIZATION, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_BOOST_MPI],
[
AC_REQUIRE([AXES_MPI])
AC_REQUIRE([AXES_BOOST_SERIALIZATION])

dnl save current flags
axes_boost_mpi_saved_LIBS="$LIBS"
axes_boost_mpi_saved_LDFLAGS="$LDFLAGS"
axes_boost_mpi_saved_CPPFLAGS="$CPPFLAGS"

dnl and preload serialization and MPI libs
LIBS="$LIBS $BOOST_SERIALIZATION_LIBS $MPI_LIBS"
LDFLAGS="$BOOST_LDFLAGS $LDFLAGS"
CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS"

AXES_BOOST_LIB(mpi, MPI, AC_LANG_PROGRAM([[
        @%:@include <boost/mpi.hpp>
]], [[
	boost::mpi::communicator world;
]]))

LIBS="$LIBS $BOOST_MPI_LIBS"
BOOST_MPI_LIBS="$BOOST_MPI_LIBS $BOOST_SERIALIZATION_LIBS $MPI_LIBS"
AC_SUBST(BOOST_MPI_LIBS)

dnl now check whether a boost mpi program can be linked
dnl this test should detect when another compiler was used to
dnl compile boost mpi

if test "x$axes_cv_boost_mpi" != "xno"; then
  AC_CACHE_CHECK([whether a boost mpi program can be linked],
  axes_cv_boost_mpi_links,
    [
    AC_LINK_IFELSE(AC_LANG_PROGRAM([[
          @%:@include <boost/mpi.hpp>
    ]], [[
          boost::mpi::communicator world;
          int i;
          world.recv(0,0,i);
    ]]), 
      axes_cv_boost_mpi_links=yes,
      axes_cv_boost_mpi=no
      axes_cv_boost_mpi_links=no
    )])
fi

dnl restore flags
LIBS="$axes_boost_mpi_saved_LIBS"
LDFLAGS="$axes_boost_mpi_saved_LDFLAGS"
CPPFLAGS="$axes_boost_mpi_saved_CPPFLAGS"

])
