#!/bin/sh

if test "$1" = "-wipe"; then
    # instead of setting up, clean the boost imports from the directory
    rm -f LICENSE_1_0.txt
    rm -rf boost
    find src ! \( -name ".svn" -prune \) -a ! -name "Makefile.am" -a ! -type d -a -exec rm {} \;

    exit 1
fi

BOOSTPATH=$1

if test ! -f "LICENSE_1_0.txt" -a \( "x$BOOSTPATH" = "x" -o ! -f "$BOOSTPATH/LICENSE_1_0.txt" \) ; then
    echo "usage: ./boostrap.sh <path to boost source directory>"
    echo "this will extract the boost sources for Boost:MPI"
    exit -1
fi

if test "x$BOOSTPATH" != "x"; then
    # copy the files from boost
    mkdir -p ./boost
    for f in LICENSE_1_0.txt boost/mpi.hpp boost/mpi; do
        cp -r $BOOSTPATH/$f ./$f;
    done
    cp -r $BOOSTPATH/libs/mpi/src .
fi

# and generate the autoconf/automake files

rm -f config.cache acconfig.h

# make directory for autoconf-specific helper scripts
if test ! -d build-aux; then
        mkdir build-aux
fi

if [ -n "`which libtoolize`" ]; then
    libtoolize --copy
fi

aclocal -I ../../../build-aux/macros	&& \
autoheader				&& \
automake --copy --foreign --add-missing	&& \
autoconf				&& \
exit 0

exit 1


