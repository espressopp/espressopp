from espressopp/buildenv:ubuntu

ARG COVERAGE
ARG EXTERNAL
ARG CC
ARG CXX
ARG BUILD_UG

#for coverage
ENV CI true
ENV TRAVIS true
ARG TRAVIS_BRANCH
ARG TRAVIS_JOB_NUMBER
ARG TRAVIS_PULL_REQUEST
ARG TRAVIS_JOB_ID
ARG TRAVIS_TAG
ARG TRAVIS_REPO_SLUG
ARG TRAVIS_COMMIT
ARG TRAVIS_COMMIT_RANGE
ARG CODESTYLE

RUN rm -rf /home/espressopp/.ccache
COPY espressopp /home/espressopp/espressopp
COPY ccache/ /home/espressopp/.ccache
USER root
RUN chown -R espressopp:espressopp /home/espressopp/espressopp /home/espressopp/.ccache
USER espressopp

WORKDIR /home/espressopp/espressopp
RUN maintainer/check_codestyle.sh ${CODESTYLE} ${TRAVIS_BRANCH} ${TRAVIS_PULL_REQUEST} ${TRAVIS_COMMIT} ${TRAVIS_COMMIT_RANGE}

RUN if [ "${CODESTYLE}" = ON ]; then \
      exit 0; \
  else \
      mkdir build && \
      cd build && \
      ccache -z && \
      cmake -DCMAKE_INSTALL_PREFIX=/usr -DEXTERNAL_BOOST=$EXTERNAL -DEXTERNAL_MPI4PY=$EXTERNAL -DWITH_XTC=ON ${COVERAGE:+-DUSE_GCOV=ON} .. && \
      make VERBOSE=1 -j2 all ${BUILD_UG:+ug doc ug-pdf} && \
      ccache -s && \
      make test && \
      make install DESTDIR=${PWD} && \
      cd .. && \
      if [ ${COVERAGE} ]; then \
        if [ ${CC} = clang ]; then \
          $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
        else \
          $HOME/.local/bin/codecov; \
        fi; \
      fi; \
  fi
WORKDIR /home/espressopp/espressopp
USER root
RUN [ "${CODESTYLE}" = ON ] && exit 0 || make install
USER espressopp
WORKDIR ../examples
