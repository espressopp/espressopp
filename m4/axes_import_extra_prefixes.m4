#
# SYNOPSIS
#
#   AXES_IMPORT_EXTRA_PREFIXES
#
# DESCRIPTION
#
#   This macro is a convenience macro that helps the user to add
#   nonstandard locations of libraries and programs, so that autoconf
#   can find them.
#   Instead of having to specify the CPPFLAGS and LDFLAGS of such
#   libraries/programs, this macro will read directory names from the
#   environment variable EXTRA_PREFIXES. CPPFLAGS will be
#   extended by -I<extra prefix>/include (if the directory exists), 
#   LDFLAGS by -L<extra prefix>/lib, -L<extra prefix>/lib32, 
#   -L<extra prefix>/lib64 or -L<extra prefix>/shlib, if any of them
#   exists.
#
#   This macro modifies:
#
#     CPPFLAGS
#     LDFLAGS
#
#
# LAST MODIFICATION
#
#   2008-08-29
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

AC_DEFUN([AXES_IMPORT_EXTRA_PREFIXES],[

AC_ARG_VAR(EXTRA_PREFIXES, [list of directories where software is installed;  All paths,
extended by /include, are added to the include paths in CPPFLAGS, and to LDFLAGS
extended by /lib, /lib32, /lib64 or /shlib, if existent])

for axes_extra_dir in $EXTRA_PREFIXES; do
    dnl look for possible includes
    for axes_extra_prefix in $EXTRA_PREFIXES; do
        axes_extra_path="$axes_extra_prefix/include"
        if test -d "$axes_extra_path"; then
            dnl make sure we import each path only once
            case " $CPPFLAGS " in
                *" -I$axes_extra_path "*) ;;
                *) CPPFLAGS="$CPPFLAGS -I$axes_extra_path" ;; 
            esac
        fi
    done
    dnl and for possible libpaths
    for axes_extra_prefix in $EXTRA_PREFIXES; do
        for axes_extra_tail in lib lib64 lib32 shlib; do
            axes_extra_path="$axes_extra_prefix/$axes_extra_tail"
            if test -d "$axes_extra_path"; then
                dnl make sure we extra each path only once
                case " $LDFLAGS " in
                    *" -L$axes_extra_path "*) ;;
                    *) LDFLAGS="$LDFLAGS -L$axes_extra_path" ;; 
                esac
            fi
        done
    done
done

])

