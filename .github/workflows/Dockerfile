from espressopp/buildenv:ubuntu

COPY espressopp/ /home/espressopp/espressopp
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip3 install -r espressopp/requirements.txt
RUN cmake -S espressopp -B espressopp-build -DCMAKE_BUILD_TYPE=Release -DESPP_WERROR=ON
RUN cmake --build espressopp-build -j2
ENV PYTHONPATH=/home/espressopp/espressopp-build:/home/espressopp/espressopp-build/contrib
