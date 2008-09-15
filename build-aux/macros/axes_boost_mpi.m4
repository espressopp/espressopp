#
# SYNOPSIS
#
#   AXES_BOOST_MPI
#
# DESCRIPTION
#
#   Test for the mpi library from the Boost C++ libraries.
#   The macro requires a preceding call to AXES_BOOST_BASE,
#   ACX_BOOST_SERIALIZATION and ACX_MPI.
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
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
# LAST MODIFICATION
#
#   2008-08-15
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
dnl and preload serialization and MPI libs
LIBS="$LIBS $BOOST_SERIALIZATION_LIBS $MPI_LIBS"

AXES_BOOST_LIB(mpi, MPI, AC_LANG_PROGRAM([[
        @%:@include <boost/mpi.hpp>
]], [[
        boost::mpi::communicator comm;
        return 0;
]]))

BOOST_MPI_LIBS="$BOOST_MPI_LIBS $BOOST_SERIALIZATION_LIBS $MPI_LIBS"
AC_SUBST(BOOST_MPI_LIBS)

dnl restore flags
LIBS="$axes_boost_mpi_saved_LIBS"

])
