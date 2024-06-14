FROM alexo91/alma9-base:latest as builder0
LABEL stage="builder"

ENV GEANT4_VERSION 11.2.1 
#maybe set this in the docker-compose file

#compile and install geant4
RUN wget https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VERSION}/geant4-v${GEANT4_VERSION}.tar.gz \
    && mkdir -p /opt/geant4/src \
    && mkdir /opt/geant4/build \
    && mkdir /opt/geant4/install \
    && mkdir /opt/geant4/data \
    && tar -xvzf geant4-v${GEANT4_VERSION}.tar.gz -C /opt/geant4/src \
    && cd /opt/geant4/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/geant4/install \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_INSTALL_DATADIR=/opt/geant4/data \
        -DGEANT4_BUILD_MULTITHREADED=OFF \
        -DGEANT4_INSTALL_EXAMPLES=OFF \
        -DGEANT4_USE_SYSTEM_EXPAT=OFF \
        -DBUILD_STATIC_LIBS=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DGEANT4_USE_OPENGL_X11=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_QT=ON \
        ../src/geant4-v${GEANT4_VERSION} \
    && make -j$(nproc) \
    && make install

#compile and install root
RUN mkdir -p /opt/root/ \
    && mkdir /opt/root/build \
    && mkdir /opt/root/install \
    && cd /opt/root/ \
    && git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git src \
    && cd /opt/root/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/root/install \
        -DBUILD_STATIC_LIBS=ON \
        -DBUILD_SHARED_LIBS=OFF \
        ../src \
    && cmake --build . --target install -- -j$(nproc)


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



#--------------------------------------------------

FROM alexo91/alma9-base:latest

COPY --from=builder0 /opt/geant4/install /opt/geant4/
COPY --from=builder0 /opt/root/install /opt/root/
COPY --from=builder0 /opt/GROOT/install /opt/GROOT/

COPY ./entrypoint.sh /
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
