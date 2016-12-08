from espressopp/buildenv:latest

RUN rm -rf /home/espressopp/.ccache
COPY espressopp /home/espressopp/espressopp
COPY ccache/ /home/espressopp/.ccache
USER root
RUN chown -R espressopp:espressopp /home/espressopp/espressopp /home/espressopp/.ccache
USER espressopp

WORKDIR /home/espressopp/espressopp
RUN mkdir build
WORKDIR build
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr ..
RUN make VERBOSE=1 -j2
RUN make test
RUN make install DESTDIR=${PWD}
USER root
RUN make install
USER espressopp
WORKDIR ../examples
