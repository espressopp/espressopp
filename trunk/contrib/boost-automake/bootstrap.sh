#!/bin/bash


BOOSTPATH=$1

# find the script directory
scriptdir=${0%/*}
if test "$scriptdir" = "$0"; then
    scriptdir=.
fi
scriptdir=`cd $scriptdir && pwd`

# where to find the autoconf macros
macrodir="$scriptdir/../../build-aux/macros"

###################### configuration

# operate in the boost directory, otherwise
cd "$BOOSTPATH"

# check that this seems to be boost
if test ! -f "LICENSE_1_0.txt"; then
    echo "usage: ./boostrap.sh <path to boost source directory that should get be autoconf'd>"
    echo "this will generate an autoconf interface for boost libraries needs"
    exit -1
fi

# make directory for autoconf-specific helper scripts
if test ! -d "build-aux"; then
    mkdir "build-aux"
fi

# copy Makefile.ams
for f in $scriptdir/Makefile*.am; do
    mfile=${f#$scriptdir/}
    cp $f $mfile
done

# copy configure.ac
sed -e "s,@macrodir@,$macrodir,g" $scriptdir/configure.ac > configure.ac
# copy ltmain.sh
cp $macrodir/../ltmain.sh build-aux
# copy acconfig.hpp.in
cp $scriptdir/acconfig.hpp.in .

# finally, create configure & Co.

rm -f config.cache acconfig.hpp

test -f README    || ln -s index.htm README
test -f AUTHORS   || ln -s people/people.htm AUTHORS
test -f NEWS      || ln -s index.htm NEWS
test -f ChangeLog || ln -s index.htm ChangeLog

aclocal -I $macrodir                    && \
automake --copy --foreign --add-missing	&& \
autoconf                                || \
( echo "auto-generation failed" 1>&2; exit 1 )

exit 0


