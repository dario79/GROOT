#!/bin/bash

source /opt/geant4/install/bin/geant4.sh
source /opt/root/install/bin/thisroot.sh

export PATH="/opt/GROOT/install/bin:$PATH"

exec "$@"
