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

RUN rm -rf /home/espressopp/.ccache
COPY espressopp /home/espressopp/espressopp
COPY ccache/ /home/espressopp/.ccache
USER root
RUN chown -R espressopp:espressopp /home/espressopp/espressopp /home/espressopp/.ccache
USER espressopp

WORKDIR /home/espressopp/espressopp
RUN mkdir build

WORKDIR build
RUN ccache -z
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr -DEXTERNAL_BOOST=$EXTERNAL -DEXTERNAL_MPI4PY=$EXTERNAL -DWITH_XTC=ON ${COVERAGE:+-DUSE_GCOV=ON} ..
RUN make VERBOSE=1 -j2 all ${BUILD_UG:+ug doc ug-pdf}
RUN ccache -s
RUN if [ -z ${COVERAGE} ] && [ -z $BUILD_UG ]; then make test CTEST_OUTPUT_ON_FAILURE=1 TEST_H5MD=ON; fi
RUN make install DESTDIR=${PWD}
RUN cd .. && if [ ${COVERAGE} ]; then \
  if [ ${CC} = clang ]; then \
    $HOME/.local/bin/codecov --gcov-exec "llvm-cov gcov"; \
  else \
    $HOME/.local/bin/codecov; \
  fi; \
fi && cd -
USER root
RUN make install
USER espressopp
WORKDIR ../examples
