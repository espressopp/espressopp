#!/bin/sh

BOOSTPATH=$1

if test "x$BOOSTPATH" = "x" -o ! -f "$BOOSTPATH/LICENSE_1_0.txt"; then
    echo "usage: ./boostrap.sh <path to boost source directory>"
    echo "this will extract the boost sources for Boost:MPI"
    exit -1
fi

# copy the files from boost
mkdir -p ./boost
for f in LICENSE_1_0.txt boost/mpi.hpp boost/mpi; do
    cp -r $BOOSTPATH/$f ./$f;
done
cp -r $BOOSTPATH/libs/mpi/src .

# and generate the autoconf/automake files
libtoolize --automake
autoheader
aclocal
automake --foreign --add-missing
autoconf
