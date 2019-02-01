from espressopp/buildenv:ubuntu

ARG COVERAGE
ARG EXTERNAL
ARG CC
ARG CXX

COPY . /home/espressopp/espressopp/

WORKDIR /home/espressopp/espressopp
USER root

RUN mkdir build

WORKDIR build
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr -DEXTERNAL_BOOST=$EXTERNAL -DEXTERNAL_MPI4PY=$EXTERNAL -DWITH_XTC=ON ..
RUN make -j4 all
RUN make install
WORKDIR ../examples
