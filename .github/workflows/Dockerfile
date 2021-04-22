from espressopp/buildenv:ubuntu

COPY espressopp/ /home/espressopp/espressopp
RUN cmake -S espressopp -B espressopp-build -DCMAKE_BUILD_TYPE=Release -DEXTERNAL_MPI4PY=OFF
RUN cmake --build espressopp-build -j2
ENV PYTHONPATH=/home/espressopp/espressopp-build:/home/espressopp/espressopp-build/contrib
