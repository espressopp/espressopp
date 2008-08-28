#!/bin/sh

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

aclocal -I ../../../config		&& \
autoheader				&& \
automake --copy --foreign --add-missing	&& \
autoconf				&& \
exit 0

exit 1


