#
# SYNOPSIS
#
#   AXES_MPI
#
# DESCRIPTION
#
#   This macro checks whether MPI (Message Passing Interface) works, a
#   standard API for parallel process communication
#   (see http://www-unix.mcs.anl.gov/mpi/). This macro does _not_ test
#   for a MPI compiler, since there is a too big variety of them; it rather
#   just tests whether the current compiler is able to compile and link
#   MPI programs. In addition, this macro checks whether linking against
#   libmpi is necessary.
#   
#   This macro calls:
#
#     AC_SUBST(MPI_LIBS)
#
#   And defines:
#
#     HAVE_MPI
#
#   And sets the cache value:
#
#     axes_cv_mpi
#
# LAST MODIFICATION
#
#   2008-08-27
#
# COPYLEFT
#
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#   based on ACX_MPI, copyrighted as follows:
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Julian C. Cummings <cummings@cacr.caltech.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AXES_MPI], [

AC_CACHE_CHECK(whether MPI is available,
                axes_cv_mpi,[
    dnl save current flags
    axes_mpi_libs_saved="$LIBS"

    dnl check whether we need -lmpi or no
    for lib in yes -lmpi; do
        if test "$lib" != "yes"; then
            LIBS="$axes_mpi_libs_save $lib"
        fi
        AC_LINK_IFELSE(AC_LANG_PROGRAM([[
            @%:@include <mpi.h>
        ]], [[
            MPI_Init((int *)0, (char ***)(0));
            return 0;
        ]]), axes_cv_mpi=$lib, axes_cv_mpi=no)
        if test $axes_cv_mpi != no; then break; fi
    done

    dnl restore flags
    LIBS="$axes_mpi_libs_saved"

    dnl fixup output if -lmpi is necessary
    if test $axes_cv_mpi != no; then
        if test $axes_cv_mpi != yes; then
            axes_cv_mpi="requires $axes_cv_mpi"
        fi
    fi
])

if test "x$axes_cv_mpi" != "xno"; then
    AC_DEFINE(HAVE_MPI,1, [define if MPI is available])
    if test "x$axes_cv_mpi" != "xyes"; then
        MPI_LIBS=`echo $axes_cv_mpi | sed -e 's/^requires //'`
    fi
    AC_SUBST(MPI_LIBS)
fi
])
