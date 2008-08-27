#
# SYNOPSIS
#
#   AXES_IMPORT_PREFIXES
#
# DESCRIPTION
#
#   A macro to set CPPFLAGS and LDFLAGS from a prefix. E.g., if you
#   configure libxx with --prefix=/home/axel/software, then EXTRA_PREFIXES
#   will add includes for /home/axel/software/lib and .../include, enabling
#   use of libxx. Note that lib can also be lib32, lib64 or shlib, depending
#   what is available on your system.
#
#   This macro sets:
#
#     CPPFLAGS
#     LDFLAGS
#
#
# LAST MODIFICATION
#
#   2008-08-27
#
# COPYLEFT
#
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
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

AC_DEFUN([AXES_IMPORT_PREFIXES],[

AC_ARG_VAR(EXTRA_PREFIXES, [list of locations where software is installed. All paths,
extended by /include, are added to the include paths in CPPFLAGS, and to LDFLAGS
extended by /lib, /lib32, /lib64 or /shlib, if existent])

for axes_import_dir in $EXTRA_PREFIXES; do
    dnl look for possible includes
    for axes_import_prefix in $EXTRA_PREFIXES; do
        axes_import_path="$axes_import_prefix/include"
        if test -d "$axes_import_path"; then
            dnl make sure we import each path only once
            case " $CPPFLAGS " in
                *" -I$axes_import_path "*) ;;
                *) CPPFLAGS="$CPPFLAGS -I$axes_import_path" ;; 
            esac
        fi
    done
    dnl and for possible libpaths
    for axes_import_prefix in $EXTRA_PREFIXES; do
        for axes_import_tail in lib lib64 lib32 shlib; do
            axes_import_path="$axes_import_prefix/$axes_import_tail"
            if test -d "$axes_import_path"; then
                dnl make sure we import each path only once
                case " $LDFLAGS " in
                    *" -L$axes_import_path "*) ;;
                    *) LDFLAGS="$LDFLAGS -L$axes_import_path" ;; 
                esac
            fi
        done
    done
done

])

