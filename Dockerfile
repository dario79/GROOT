FROM almalinux:latest

RUN yum install epel-release -y && yum update -y

RUN yum install -y \
    gcc \
    gcc-c++ \
    make \
    cmake \
    wget \
    xerces-c-devel \
    qt5-devel \
    libXmu-devel \
    git \
    nano


ENV GEANT4_VERSION 11.2.1 
#maybe set this in the docker-compose file

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

RUN mkdir -p /opt/root/ \
    && mkdir /opt/root/build \
    && mkdir /opt/root/install \
    && cd /opt/root/ \
    && git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git src \
    && cd /opt/root/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt/root/install \
        -Dxrootd=OFF \
        -Dbuiltin_xrootd=OFF \
        -Dminimal=ON \
        ../src \
    && cmake --build . --target install -- -j$(nproc)
