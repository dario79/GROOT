#!/bin/bash

source /opt/geant4/install/bin/geant4.sh
source /opt/root/install/bin/thisroot.sh

cd /opt/GROOT/build
cmake -DCMAKE_INSTALL_PREFIX=/opt/GROOT/install -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF ../source
make -j$(nproc)
make install