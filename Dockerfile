FROM alexo91/groot-base:alpha-0.2 as builder0
LABEL stage="builder"
SHELL ["/bin/bash", "-c"]

RUN source /opt/geant4/install/bin/geant4.sh \
    && source /opt/root/install/bin/thisroot.sh \
    && geant4-config --cflags \
    && root-config --cflags


RUN mkdir -p /opt/GROOT/source \
    && mkdir /opt/GROOT/build \
    && mkdir /opt/GROOT/install 

COPY ./GROOT.cc /opt/GROOT/source/
COPY ./CMakeLists.txt /opt/GROOT/source/
COPY ./inputFileDetectors.txt /opt/GROOT/source/
COPY ./run.mac /opt/GROOT/source/
COPY ./vis.mac /opt/GROOT/source/
COPY ./src/ /opt/GROOT/source/src/
COPY ./include/ /opt/GROOT/source/include/
COPY ./cmake/ /opt/GROOT/source/cmake/
COPY ./build-GROOT.sh /opt/GROOT/source/

RUN chmod +x /opt/GROOT/source/build-GROOT.sh
RUN "/opt/GROOT/source/build-GROOT.sh"

#RUN cd /opt/GROOT/build \
#    && cmake \
#        -DCMAKE_INSTALL_PREFIX=/opt/GROOT/install \
#        -DBUILD_STATIC_LIBS=ON \
#        -DBUILD_SHARED_LIBS=OFF \
#        ../source \
#    && make -j$(nproc) \
#    && make install

#COPY ./entrypoint.sh /
#RUN chmod +x /entrypoint.sh
#ENTRYPOINT ["/entrypoint.sh"]
