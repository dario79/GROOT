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
/// \file electromagnetic/Groot1/src/PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
// $Id: PrimaryGeneratorMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"
#include "G4UIdirectory.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                             PrimaryGeneratorAction* Gun)
:G4UImessenger(),Action(Gun),
 fGunDir(0),
 fDefaultCmd(0),
 fGrootDir(0),
 fEnergyCmd(0),
 frbeamCmd(0),
 fP1ACmd(0),
 fP1ZCmd(0),
 fP2ACmd(0),
 fP2ZCmd(0),
 fP3ACmd(0),
 fP3ZCmd(0),
 fmassP1Cmd(0),
 fmassP2Cmd(0),
 fmassP3Cmd(0),
 fmassProjectileCmd(0),
 fmassTargetCmd(0),
 fStateCmd(0),
 fRndmCmd(0)
{
    fGrootDir = new G4UIdirectory("/groot/");
    fGrootDir->SetGuidance("commands specific to this example");

    fGunDir = new G4UIdirectory("/groot/gun/");
    fGunDir->SetGuidance("gun control");

    fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/groot/gun/energy",this);
    fEnergyCmd->SetGuidance("Set energy of the beam");
    fEnergyCmd->SetParameterName("Size",false);
    fEnergyCmd->SetRange("Size>0.");
    fEnergyCmd->SetUnitCategory("Energy");
    fEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    frbeamCmd = new G4UIcmdWithADoubleAndUnit("/groot/gun/rbeam",this);
    frbeamCmd->SetGuidance("Set size of the beam");
    frbeamCmd->SetParameterName("Size",false);
    frbeamCmd->SetRange("Size>0.");
    frbeamCmd->SetUnitCategory("Length");
    frbeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP1ACmd = new G4UIcmdWithAnInteger("/groot/gun/particle1A",this);
    fP1ACmd->SetGuidance("Set A of particle1");
    fP1ACmd->SetParameterName("A",false);
    fP1ACmd->SetRange("A>=0");
    fP1ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP1ZCmd = new G4UIcmdWithAnInteger("/groot/gun/particle1Z",this);
    fP1ZCmd->SetGuidance("Set Z of particle1");
    fP1ZCmd->SetParameterName("Z",false);
    fP1ZCmd->SetRange("Z>=0");
    fP1ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP2ACmd = new G4UIcmdWithAnInteger("/groot/gun/particle2A",this);
    fP2ACmd->SetGuidance("Set number A of particle2");
    fP2ACmd->SetParameterName("A",false);
    fP2ACmd->SetRange("A>=0");
    fP2ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP2ZCmd = new G4UIcmdWithAnInteger("/groot/gun/particle2Z",this);
    fP2ZCmd->SetGuidance("Set Z of particle2");
    fP2ZCmd->SetParameterName("Z",false);
    fP2ZCmd->SetRange("Z>=0");
    fP2ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP3ACmd = new G4UIcmdWithAnInteger("/groot/gun/particle3A",this);
    fP3ACmd->SetGuidance("Set number A of particle3");
    fP3ACmd->SetParameterName("A",false);
    fP3ACmd->SetRange("A>=0");
    fP3ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fP3ZCmd = new G4UIcmdWithAnInteger("/groot/gun/particle3Z",this);
    fP3ZCmd->SetGuidance("Set Z of particle3");
    fP3ZCmd->SetParameterName("Z",false);
    fP3ZCmd->SetRange("Z>=0");
    fP3ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fmassP1Cmd = new G4UIcmdWithADouble("/groot/gun/massPart1",this);
    fmassP1Cmd->SetGuidance("Set mass of Particle1.");
    fmassP1Cmd->SetParameterName("m1",false);
    fmassP1Cmd->SetRange("m1>=0");
    fmassP1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fmassP2Cmd = new G4UIcmdWithADouble("/groot/gun/massPart2",this);
    fmassP2Cmd->SetGuidance("Set mass of Particle2.");
    fmassP2Cmd->SetParameterName("m2",false);
    fmassP2Cmd->SetRange("m2>=0");
    fmassP2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fmassP3Cmd = new G4UIcmdWithADouble("/groot/gun/massPart3",this);
    fmassP3Cmd->SetGuidance("Set mass of Particle3.");
    fmassP3Cmd->SetParameterName("m3",false);
    fmassP3Cmd->SetRange("m3>=0");
    fmassP3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fmassProjectileCmd = new G4UIcmdWithADouble("/groot/gun/massProjectile",this);
    fmassProjectileCmd->SetGuidance("Set mass of Projectile.");
    fmassProjectileCmd->SetParameterName("mp",false);
    fmassProjectileCmd->SetRange("mp>=0");
    fmassProjectileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fmassTargetCmd = new G4UIcmdWithADouble("/groot/gun/massTarget",this);
    fmassTargetCmd->SetGuidance("Set mass of Target.");
    fmassTargetCmd->SetParameterName("mt",false);
    fmassTargetCmd->SetRange("mt>=0");
    fmassTargetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fStateCmd = new G4UIcmdWithADoubleAndUnit("/groot/gun/setState",this);
    fStateCmd->SetGuidance("Set state of the beam.");
    fStateCmd->SetParameterName("nr",false);
    fStateCmd->SetRange("nr>=0");
    fStateCmd->SetUnitCategory("Energy");
    fStateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDefaultCmd = new G4UIcmdWithoutParameter("/groot/gun/setDefault",this);
    fDefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
    fDefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fRndmCmd = new G4UIcmdWithADoubleAndUnit("/groot/gun/rndm",this);
    fRndmCmd->SetGuidance("random lateral extension on the beam");
    fRndmCmd->SetParameterName("rBeam",false);
    fRndmCmd->SetRange("rBeam>=0.");
    fRndmCmd->SetUnitCategory("Length");
    fRndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
    delete fDefaultCmd;
    delete fRndmCmd;
    delete fGunDir;
    delete fEnergyCmd;
    delete frbeamCmd;
    delete fP1ACmd;
    delete fP1ZCmd;
    delete fP2ACmd;
    delete fP2ZCmd;
    delete fP3ACmd;
    delete fP3ZCmd;
    delete fmassP1Cmd;
    delete fmassP2Cmd;
    delete fmassP3Cmd;
    delete fmassProjectileCmd;
    delete fStateCmd;
    delete fmassTargetCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

    if( command == fP1ACmd )
    {
        Action->SetParticle1A(fP1ACmd->GetNewIntValue(newValue));
    }

    if( command == fP1ZCmd )
    {
        Action->SetParticle1Z(fP1ZCmd->GetNewIntValue(newValue));
    }

    if( command == fP2ACmd )
    {
        Action->SetParticle2A(fP2ACmd->GetNewIntValue(newValue));
    }

    if( command == fP2ZCmd )
    {
        Action->SetParticle2Z(fP2ZCmd->GetNewIntValue(newValue));
    }
    if( command == fP3ACmd )
    {
        Action->SetParticle3A(fP3ACmd->GetNewIntValue(newValue));
    }

    if( command == fP3ZCmd )
    {
        Action->SetParticle3Z(fP3ZCmd->GetNewIntValue(newValue));
    }

    if(command == fmassP1Cmd)
    {
        Action->SetParticle1Mass(fmassP1Cmd->GetNewDoubleValue(newValue));
    }

    if(command == fmassP2Cmd)
    {
        Action->SetParticle2Mass(fmassP2Cmd->GetNewDoubleValue(newValue));
    }

    if(command == fmassP3Cmd)
    {
        Action->SetParticle3Mass(fmassP3Cmd->GetNewDoubleValue(newValue));
    }

    if(command == fmassProjectileCmd)
    {
        Action->SetProjectileMass(fmassProjectileCmd->GetNewDoubleValue(newValue));
    }

    if(command == fmassTargetCmd)
    {
        Action->SetTargetMass(fmassTargetCmd->GetNewDoubleValue(newValue));
    }

    if(command == fStateCmd)
    {
        Action->SetState(fStateCmd->GetNewDoubleValue(newValue));
    }

    if(command == fEnergyCmd)
    {
      Action->SetEnergy(fEnergyCmd->GetNewDoubleValue(newValue));
    }

    if(command == frbeamCmd)
    {
      Action->Setrbeam(frbeamCmd->GetNewDoubleValue(newValue));
    }

    if (command == fDefaultCmd)
    {
        Action->SetDefaultKinematic();
    }

    if (command == fRndmCmd)
    {
        Action->SetRndmBeam(fRndmCmd->GetNewDoubleValue(newValue));
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

