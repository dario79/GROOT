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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 98304 2016-07-05 15:29:08Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction0.hh"
#include "PrimaryGeneratorAction2.hh"
#include "PrimaryGeneratorAction3.hh"
#include "PrimaryGeneratorAction4.hh"
#include "PrimaryGeneratorAction5.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "SteppingAction.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction* PrimaryGeneratorAction::instance = 0;
DetectorConstruction* det;

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fAction0(0),
   fAction2(0),
   fAction3(0),
   fAction4(0),
   fAction5(0),
   fGunMessenger(0),
   fRndmBeam(0)
{

    det = DC; // get the pointer to the current DC
    // Definition of the General particle Source

    fParticleGun  = new G4ParticleGun(1);
    SetDefaultKinematic();  // 2 = 2body, 3 = 3body  else if standard
    fRndmBeam = 0.;
  
  fAction0 = new PrimaryGeneratorAction0(fParticleGun);
  fAction2 = new PrimaryGeneratorAction2(fParticleGun);
  fAction3 = new PrimaryGeneratorAction3(fParticleGun);
  fAction4 = new PrimaryGeneratorAction4(fParticleGun);
  fAction5 = new PrimaryGeneratorAction5(fParticleGun);

  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fAction0;
  delete fAction2;
  delete fAction3;
  delete fAction4;
  delete fAction5;
  delete fParticleGun;  
  delete fGunMessenger;      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematic()
{

    ProjectileMass=0;           // massa proiettile             //
    TargetMass=111.904826;   // massa target                 //
    Part1Mass=4.002603;   // massa prima part. uscente    //
    Part2Mass=107.90418;    // massa seconda part. uscente  //
    Part3Mass=1.;    // massa terza part. uscente  //
    Part3A = 3;  Part3Z = 3;
    Part2A = 60; Part2Z = 48;
    Part1A = 2;  Part1Z = 2;
    fuffa = 2;
    targetphi = 0.*deg;// -35.*deg;
    g_Energy = 5.; // MeV
    outputType = "both";
    rbeam = 0.75*cm; // 0.5*mm; /* sin(targetrot*deg)* */
    state = 0 *MeV;
    outputName = "GROOT";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetState(G4double states)
{
    state = states;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetEnergy(G4double energy)
{
    g_Energy = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SelectOutputType(G4String t)
{
    outputType = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SelectOutputName(G4String name)
{
    outputName = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle1A(G4int A)
{
    Part1A = A;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle1Z(G4int Z)
{
    Part1Z = Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle2A(G4int A)
{
    Part2A = A;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle2Z(G4int Z)
{
    Part2Z = Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle3A(G4int A)
{
    Part3A = A;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle3Z(G4int Z)
{
    Part3Z = Z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle1Mass(G4double mass)
{
    Part1Mass = mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle2Mass(G4double mass)
{
    Part2Mass = mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetParticle3Mass(G4double mass)
{
    Part3Mass = mass;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetProjectileMass(G4double mass)
{
    ProjectileMass = mass;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetTargetMass(G4double mass)
{
    TargetMass = mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::Setrbeam(G4double beam)
{
   rbeam = beam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    fuffa = det->GetReactionType();

    switch(fuffa)
    {
        case 0:
            fAction0->GeneratePrimaries(anEvent, rbeam, g_Energy, outputType, outputName, det);
        break;
        case 2:
            fAction2->GeneratePrimaries(anEvent, state, g_Energy, outputType, outputName, Part1A, Part1Z, 						Part2A, Part2Z, Part1Mass, Part2Mass, ProjectileMass, TargetMass, 						rbeam, det);
        break;
        case 3:
            fAction3->GeneratePrimaries(anEvent, state, g_Energy, outputType, outputName, Part1A, Part1Z, 						Part2A, Part2Z, Part3A, Part3Z, Part1Mass, Part2Mass, Part3Mass, 						ProjectileMass, TargetMass, rbeam, det);
        break;
        case 4:
            fAction4->GeneratePrimaries(anEvent, rbeam, g_Energy, outputType, outputName, det);
        break;
        case 5:
            fAction5->GeneratePrimaries(anEvent, g_Energy, outputName, rbeam);
        break;
        default:
            G4cerr << "Invalid generator fAction" << G4endl;
    }
  }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
