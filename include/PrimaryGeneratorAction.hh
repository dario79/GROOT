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
/// \file electromagnetic/TestEm1/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
// $Id: PrimaryGeneratorAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
//#include "G4GeneralParticleSource.hh"

#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TROOT.h"
#include <random>
class G4Event;
class G4GeneralParticleSource;
class G4SingleParticleSource;
class PrimaryGeneratorMessenger;
class PrimaryGeneratorAction;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

   static const double     pi  = 3.14159265358979323846;
   static const double  twopi  = 2*pi;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);    
   ~PrimaryGeneratorAction();

  public:
    G4int fuffa;
    void SetDefaultKinematic(G4int);
    void SetRndmBeam(G4double val)  {fRndmBeam = val;}   
    void GeneratePrimaries(G4Event*);
    G4double GetPrimaryEnergy();
    static PrimaryGeneratorAction* GetInstance();


//    G4SingleParticleSource* GetParticleGun();
    G4ParticleGun* GetParticleGun() {return fParticleGun;}
    G4double                   energiaIniziale;
    G4double                   targetphi;
    G4double                   g_Energy;
  private:
    G4ParticleGun*             fParticleGun;
//    G4GeneralParticleSource* fParticleGun;
    G4double                   fRndmBeam;
    PrimaryGeneratorMessenger* fGunMessenger;
    static PrimaryGeneratorAction* instance;
    DetectorConstruction*      fDetector;


    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


