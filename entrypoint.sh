#!/bin/bash

source /opt/geant4/bin/geant4.sh
source /opt/root/bin/thisroot.sh
source /opt/GROOT/source/build-GROOT.sh

export PATH="/opt/GROOT/install/bin:$PATH"

exec "$@"
