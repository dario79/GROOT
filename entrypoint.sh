#!/bin/bash

source /opt/geant4/bin/geant4.sh
source /opt/root/bin/thisroot.sh


export PATH="/opt/GROOT/bin:$PATH"

exec "$@"
