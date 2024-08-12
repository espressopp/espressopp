from espressopp/buildenv:ubuntu

COPY espressopp/ /home/espressopp/espressopp
RUN pip3 install -r espressopp/requirements.txt --break-system-packages
RUN cmake -S espressopp -B espressopp-build -DCMAKE_BUILD_TYPE=Release -DESPP_WERROR=ON
RUN cmake --build espressopp-build -j2
ENV PYTHONPATH=/home/espressopp/espressopp-build:/home/espressopp/espressopp-build/contrib
