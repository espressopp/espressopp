#!/bin/bash
#
# build.sh
# Copyright (C) 2019 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#


cp -vr docker ${HOME}
sed -i "1s/ubuntu/${DISTRO}/" ${HOME}/docker/Dockerfile
cd ../../
mv -v ${TRAVIS_REPO_SLUG} $HOME/docker
cp -r $HOME/.ccache ${HOME}/docker/ccache
docker pull $(sed -n '1s/from //p' ${HOME}/docker/Dockerfile)
docker build --build-arg COVERAGE=${COVERAGE} --build-arg EXTERNAL=${EXTERNAL} \
                --build-arg CC=${CC} --build-arg CXX=${CXX} --build-arg BUILD_UG=${BUILD_UG} \
                --build-arg TRAVIS_BRANCH=${TRAVIS_BRANCH} --build-arg TRAVIS_JOB_NUMBER=${TRAVIS_JOB_NUMBER} \
                --build-arg TRAVIS_PULL_REQUEST=${TRAVIS_PULL_REQUEST} --build-arg TRAVIS_JOB_ID=${TRAVIS_JOB_ID} \
                --build-arg TRAVIS_TAG=${TRAVIS_TAG} --build-arg TRAVIS_REPO_SLUG=${TRAVIS_REPO_SLUG} \
                --build-arg TRAVIS_COMMIT=${TRAVIS_COMMIT} \
                -t espressopp/espressopp:${DOCKER_TAG} ${HOME}/docker/ && \
   rm -rf $HOME/.ccache && \
   CON=$(docker run -d espressopp/espressopp:${DOCKER_TAG} /bin/bash) && \
   docker cp ${CON}:/home/espressopp/.ccache ${HOME}/
