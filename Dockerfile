FROM almalinux:latest

RUN dnf install epel-release -y \
    && dnf config-manager --set-enabled crb \
    && dnf update -y

RUN dnf install gcc \
    make \
    cmake \
    wget \
    xerces-c-devel \
    qt5-devel \
    binutils \
    libX11-devel libXpm-devel libXmu-devel libXft-devel libXext-devel python openssl-devel \
    xrootd-client-devel xrootd-libs-devel \
    mesa-libGL-devel mesa-libGLU-devel glew-devel ftgl-devel mysql-devel \
    fftw-devel cfitsio-devel graphviz-devel libuuid-devel \
    avahi-compat-libdns_sd-devel openldap-devel python-devel python3-numpy \
    libxml2-devel gsl-devel readline-devel qt5-qtwebengine-devel \
    git \
    nano -y


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
        -DGEANT4_BUILD_MULTITHREADED=ON \
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

SHELL ["/bin/bash", "-c"] 

RUN source /opt/geant4/install/bin/geant4.sh
RUN source /opt/root/install/bin/thisroot.sh
RUN geant4-config --cflags
RUN root-config --cflags

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

RUN cd /opt/GROOT/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/GROOT/install \
        -DBUILD_STATIC_LIBS=ON \
        -DBUILD_SHARED_LIBS=OFF \
        ../source \
    && make -j$(nproc) \
    && make install
