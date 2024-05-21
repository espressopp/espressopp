from espressopp/buildenv:ubuntu

COPY espressopp/ /home/espressopp/espressopp
RUN python3 -m venv .venv
RUN source .venv/bin/activate
RUN pip3 install -r espressopp/requirements.txt
RUN cmake -S espressopp -B espressopp-build -DCMAKE_BUILD_TYPE=Release -DESPP_WERROR=ON
RUN cmake --build espressopp-build -j2
ENV PYTHONPATH=/home/espressopp/espressopp-build:/home/espressopp/espressopp-build/contrib
