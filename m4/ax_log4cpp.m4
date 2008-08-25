===========================================================================
#
# SYNOPSIS
#
#   AX_LOG4CPP([MINIMUM-VERSION])
#
# DESCRIPTION
#
#   Test for the log4cpp library of a particular version (or newer)
#
#   If no path to the installed log4cpp library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $LOG4CPP_ROOT environment variable. 
#
#   This macro calls:
#
#     AC_SUBST(LOG4CPP_CPPFLAGS) / AC_SUBST(LOG4CPP_LDFLAGS)
#
#   And sets:
#
#     HAVE_LOG4CPP
#
# LAST MODIFICATION
#
#   2008-08-08
#
# COPYLEFT
#
#   Copyright (c) 2008 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_LOG4CPP],
[
AC_ARG_WITH([log4cpp],
	AS_HELP_STRING([--with-log4cpp@<:@=DIR@:>@], [use log4cpp (default is yes) - it is possible to specify the root directory for log4cpp (optional)]),
	[
    if test "$withval" = "no"; then
		want_log4cpp="no"
    elif test "$withval" = "yes"; then
        want_log4cpp="yes"
        ac_log4cpp_path=""
    else
	    want_log4cpp="yes"
        ac_log4cpp_path="$withval"
	fi
    ],
    [want_log4cpp="yes"])


AC_ARG_WITH([log4cpp-libdir],
        AS_HELP_STRING([--with-log4cpp-libdir=LIB_DIR],
        [Force given directory for log4cpp libraries. Note that this will overwrite library path detection, so use this parameter only if default library detection fails and you know exactly where your log4cpp libraries are located.]),
        [
        if test -d $withval
        then
                ac_log4cpp_lib_path="$withval"
        else
                AC_MSG_ERROR(--with-log4cpp-libdir expected directory name)
        fi
        ],
        [ac_log4cpp_lib_path=""]
)

if test "x$want_log4cpp" = "xyes"; then
	log4cpp_lib_version_req=ifelse([$1], ,1.20.0,$1)
	log4cpp_lib_version_req_shorten=`expr $log4cpp_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
	log4cpp_lib_version_req_major=`expr $log4cpp_lib_version_req : '\([[0-9]]*\)'`
	log4cpp_lib_version_req_minor=`expr $log4cpp_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
	log4cpp_lib_version_req_sub_minor=`expr $log4cpp_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
	if test "x$log4cpp_lib_version_req_sub_minor" = "x" ; then
		log4cpp_lib_version_req_sub_minor="0"
    	fi
	WANT_LOG4CPP_VERSION=`expr $log4cpp_lib_version_req_major \* 100000 \+  $log4cpp_lib_version_req_minor \* 100 \+ $log4cpp_lib_version_req_sub_minor`
	AC_MSG_CHECKING(for log4cpplib >= $log4cpp_lib_version_req)
	succeeded=no

	dnl first we check the system location for log4cpp libraries
	dnl this location ist chosen if log4cpp libraries are installed with the --layout=system option
	dnl or if you install log4cpp with RPM
	if test "$ac_log4cpp_path" != ""; then
		LOG4CPP_LDFLAGS="-L$ac_log4cpp_path/lib"
		LOG4CPP_CPPFLAGS="-I$ac_log4cpp_path/include"
	else
		for ac_log4cpp_path_tmp in /usr /usr/local /opt /opt/local ; do
			if test -d "$ac_log4cpp_path_tmp/include/log4cpp" && test -r "$ac_log4cpp_path_tmp/include/log4cpp"; then
				LOG4CPP_LDFLAGS="-L$ac_log4cpp_path_tmp/lib"
				LOG4CPP_CPPFLAGS="-I$ac_log4cpp_path_tmp/include"
				break;
			fi
		done
	fi

    dnl overwrite ld flags if we have required special directory with
    dnl --with-log4cpp-libdir parameter
    if test "$ac_log4cpp_lib_path" != ""; then
       LOG4CPP_LDFLAGS="-L$ac_log4cpp_lib_path"
    fi

	CPPFLAGS_SAVED="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $LOG4CPP_CPPFLAGS"
	export CPPFLAGS

	LDFLAGS_SAVED="$LDFLAGS"
	LDFLAGS="$LDFLAGS $LOG4CPP_LDFLAGS"
	export LDFLAGS

	AC_LANG_PUSH(C++)
     	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
	@%:@include <log4cpp/Category.hh>
	]], [[
        log4cpp::Category& myLogger = log4cpp::Category::getRoot();
        LOG4CPP_DEBUG(myLogger, "debug message");
	]])],[
        AC_MSG_RESULT(yes)
	succeeded=yes
	found_system=yes
       	],[
       	])
	AC_LANG_POP([C++])


	dnl if we found no log4cpp with system layout we search for log4cpp libraries
	dnl built and installed without the --layout=system option or for a staged(not installed) version
	if test "x$succeeded" != "xyes"; then
		_version=0
		if test "$ac_log4cpp_path" != ""; then
			if test -d "$ac_log4cpp_path" && test -r "$ac_log4cpp_path"; then
				for i in `ls -d $ac_log4cpp_path/include/log4cpp-* 2>/dev/null`; do
					_version_tmp=`echo $i | sed "s#$ac_log4cpp_path##" | sed 's/\/include\/log4cpp-//' | sed 's/_/./'`
					V_CHECK=`expr $_version_tmp \> $_version`
					if test "$V_CHECK" = "1" ; then
						_version=$_version_tmp
					fi
					VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
					LOG4CPP_CPPFLAGS="-I$ac_log4cpp_path/include/log4cpp-$VERSION_UNDERSCORE"
				done
			fi
		else
			for ac_log4cpp_path in /usr /usr/local /opt /opt/local ; do
				if test -d "$ac_log4cpp_path" && test -r "$ac_log4cpp_path"; then
					for i in `ls -d $ac_log4cpp_path/include/log4cpp-* 2>/dev/null`; do
						_version_tmp=`echo $i | sed "s#$ac_log4cpp_path##" | sed 's/\/include\/log4cpp-//' | sed 's/_/./'`
						V_CHECK=`expr $_version_tmp \> $_version`
						if test "$V_CHECK" = "1" ; then
							_version=$_version_tmp
	               					best_path=$ac_log4cpp_path
						fi
					done
				fi
			done

			VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
			LOG4CPP_CPPFLAGS="-I$best_path/include/log4cpp-$VERSION_UNDERSCORE"
            if test "$ac_log4cpp_lib_path" = ""
            then
               LOG4CPP_LDFLAGS="-L$best_path/lib"
            fi

	    		if test "x$LOG4CPP_ROOT" != "x"; then
				if test -d "$LOG4CPP_ROOT" && test -r "$LOG4CPP_ROOT" && test -d "$LOG4CPP_ROOT/lib" && test -r "$LOG4CPP_ROOT/lib"; then
					AC_MSG_NOTICE(We will use a log4cpp library from $LOG4CPP_ROOT)
					LOG4CPP_CPPFLAGS="-I$LOG4CPP_ROOT/include"
					LOG4CPP_LDFLAGS="-L$LOG4CPP_ROOT/lib"
				fi
	    		fi
		fi

		CPPFLAGS="$CPPFLAGS $LOG4CPP_CPPFLAGS"
		export CPPFLAGS
		LDFLAGS="$LDFLAGS $LOG4CPP_LDFLAGS"
		export LDFLAGS

		AC_LANG_PUSH(C++)
	     	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
	        @%:@include <log4cpp/logger.h>
                using namespace log4cpp;
	        ]], [[
                LoggerPtr rootLogger = Logger::getRootLogger();
                LOG4CPP_DEBUG(rootLogger, "debug message");
		]])],[
        	AC_MSG_RESULT(yes)
		succeeded=yes
		found_system=yes
       		],[
	       	])
		AC_LANG_POP([C++])
	fi

	if test "$succeeded" != "yes" ; then
		if test "$_version" = "0" ; then
			AC_MSG_WARN([[We could not detect the log4cpp libraries (version $log4cpp_lib_version_req_shorten or higher). Please specify \$LOG4CPP_ROOT in your environment or give a PATH to --with-log4cpp option.  If you are sure you have log4cpp installed, then check your version number looking in <log4cpp/version.hpp>. See http://randspringer.de/log4cpp for more documentation.]])
		else
			AC_MSG_NOTICE([Your log4cpp libraries seems to old (version $_version).])
		fi
                HAVE_LOG4CPP=0
	else
                HAVE_LOG4CPP=1
		AC_SUBST(LOG4CPP_CPPFLAGS)
		AC_SUBST(LOG4CPP_LDFLAGS)
		AC_DEFINE(HAVE_LOG4CPP,1,[define if the log4cpp library is available])
	fi

        CPPFLAGS="$CPPFLAGS_SAVED"
       	LDFLAGS="$LDFLAGS_SAVED"
fi

])
