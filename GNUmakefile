# $Id: GNUmakefile 67981 2013-03-13 10:34:11Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GROOT
G4TARGET := $(name)
G4EXLIB := true


ifndef G4INSTALL
 G4INSTALL = "/opt/geant4/geant4.10.06.p01-install-new/usr/local/share/Geant4-10.6.1/geant4make"
# G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

LDFLAGS += $(shell /usr/bin/root-config --glibs)

CPPFLAGS += $(shell /usr/bin/root-config --cflags) -std=gnu++11

# include ../../cadmesh.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
