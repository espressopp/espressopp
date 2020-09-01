#! /bin/sh
#
# build.sh
# Copyright (C) 2019 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

DIRNAME=$(cd `dirname $0` && pwd)
DOCKERFILE=${DIRNAME}/Dockerfile.local
DIST=ubuntu

if [[ -n $1 ]]; then
    DIST=$1
fi

DIST_DOCKERFILE=${DIRNAME}/Dockerfile.${DIST}

cp -v $DOCKERFILE $DIST_DOCKERFILE
sed -i "s/buildenv:.*/buildenv:$DIST/g" $DIST_DOCKERFILE

docker build -t espressopp:$DIST -f $DIST_DOCKERFILE .

rm -v $DIST_DOCKERFILE
