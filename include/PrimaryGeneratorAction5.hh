//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file eventgenerator/particleGun/include/PrimaryGeneratorAction5.hh
/// \brief Definition of the PrimaryGeneratorAction5 class
//
//
// $Id: PrimaryGeneratorAction5.hh 67332 2013-02-14 17:12:13Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction5_h
#define PrimaryGeneratorAction5_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"

#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TROOT.h"
#include <random>


class G4GeneralParticleSource;
class G4SingleParticleSource;
class PrimaryGeneratorAction;
class DetectorConstruction;


class G4ParticleGun;
class G4Event;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction5
{
  public:
    PrimaryGeneratorAction5(G4ParticleGun*);
   ~PrimaryGeneratorAction5();

  public:
    void GeneratePrimaries(G4Event*, G4double, G4String, G4double);
    G4ParticleDefinition*      out_particle2;
    G4ParticleDefinition*      out_particle1;
    G4double                   targetphi;

  private:

    G4ParticleGun*  fParticleGun;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
