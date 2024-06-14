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
/// \file eventgenerator/particleGun/include/PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.hh 83919 2014-09-23 08:40:35Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"

#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TROOT.h"
#include <random>


class G4GeneralParticleSource;
class G4SingleParticleSource;
class PrimaryGeneratorMessenger;
class PrimaryGeneratorAction;
class DetectorConstruction;

class PrimaryGeneratorAction0;
class PrimaryGeneratorAction2;
class PrimaryGeneratorAction3;
class PrimaryGeneratorAction4;
class PrimaryGeneratorAction5;

class G4ParticleGun;
class G4Event;

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
  G4String outputType, outputName;
  G4int Part1A;
  G4int Part1Z;
  G4int Part2A;
  G4int Part2Z;
  G4int Part3A;
  G4int Part3Z;
  G4double Part1Mass;
  G4double Part2Mass;
  G4double Part3Mass;
  G4double ProjectileMass;
  G4double TargetMass;
  G4double rbeam;
  G4double state;
  void SetDefaultKinematic();
  void SetParticle2A(G4int);
  void SetParticle2Z(G4int);
  void SetParticle3Z(G4int);
  void SetParticle3A(G4int);
  void SetParticle1Mass(G4double);
  void SetParticle2Mass(G4double);
  void SetParticle3Mass(G4double);
  void SetProjectileMass(G4double);
  void SetTargetMass(G4double);
  void SetEnergy(G4double);
  void Setrbeam(G4double);
  void SetParticle1A(G4int);
  void SelectOutputName(G4String);
  void SetParticle1Z(G4int);
  void SetRndmBeam(G4double val)  {fRndmBeam = val;}
  void GeneratePrimaries(G4Event*);
  void SelectOutputType(G4String);
  void SetState(G4double);
  G4double GetPrimaryEnergy();
  static PrimaryGeneratorAction* GetInstance();


  G4ParticleGun* GetParticleGun() {return fParticleGun;}
  G4double                   energiaIniziale;
  G4double                   targetphi;
  G4double                   g_Energy;
  G4ParticleDefinition*      out_particle2;
  G4ParticleDefinition*      out_particle1;
  G4ParticleDefinition*      out_particle3;

private:
  G4ParticleGun*             fParticleGun;

  G4double                   fRndmBeam;
  PrimaryGeneratorMessenger* fGunMessenger;
  static PrimaryGeneratorAction* instance;
  DetectorConstruction*      fDetector;



  public:


    PrimaryGeneratorAction0*  GetAction0() { return fAction0; };
    PrimaryGeneratorAction2*  GetAction2() { return fAction2; };
    PrimaryGeneratorAction3*  GetAction3() { return fAction3; };
    PrimaryGeneratorAction4*  GetAction4() { return fAction4; };
    PrimaryGeneratorAction5*  GetAction5() { return fAction5; };

  private:

    PrimaryGeneratorAction0* fAction0; 
    PrimaryGeneratorAction2* fAction2;
    PrimaryGeneratorAction3* fAction3;
    PrimaryGeneratorAction4* fAction4;
    PrimaryGeneratorAction5* fAction5;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
