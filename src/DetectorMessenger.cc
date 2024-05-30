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
/// \file electromagnetic/Groot1/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),
 fDetector(Det),
 fGrootDir(0),
 fDetDir(0),
 fDetMaterCmd(0),
 fBackingMaterCmd(0),
 fFilesPathCmd(0),
 fDetThicknessCmd(0),
 fBackingThicknessCmd(0),
 fMagFieldCmd(0),
 fUpdateCmd(0),
 fBackingPosCmd(0),
 fTargetPhiCmd(0),
 fReactionTypeCmd(0)
{ 
    fGrootDir = new G4UIdirectory("/groot/");
    fGrootDir->SetGuidance("commands specific to this example");

    fDetDir = new G4UIdirectory("/groot/det/");
    fDetDir->SetGuidance("detector construction commands");
        
    fDetMaterCmd = new G4UIcmdWithAString("/groot/det/setDetMat",this);
    fDetMaterCmd->SetGuidance("Select material of the box.");
    fDetMaterCmd->SetParameterName("choice",false);
    fDetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fBackingMaterCmd = new G4UIcmdWithAString("/groot/det/setBackingMat",this);
    fBackingMaterCmd->SetGuidance("Select material of the Backing.");
    fBackingMaterCmd->SetParameterName("choice",false);
    fBackingMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fDetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/groot/det/setDetThickness",this);
    fDetThicknessCmd->SetGuidance("Set size of the box");
    fDetThicknessCmd->SetParameterName("Size",false);
    fDetThicknessCmd->SetRange("Size>0.");
    fDetThicknessCmd->SetUnitCategory("Length");
    fDetThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fBackingThicknessCmd = new G4UIcmdWithADoubleAndUnit("/groot/det/setBackingThickness",this);
    fBackingThicknessCmd->SetGuidance("Set size of the backing");
    fBackingThicknessCmd->SetParameterName("Size",false);
    fBackingThicknessCmd->SetRange("Size>0.");
    fBackingThicknessCmd->SetUnitCategory("Length");
    fBackingThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/groot/det/setField",this);
    fMagFieldCmd->SetGuidance("Define magnetic field.");
    fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
    fMagFieldCmd->SetParameterName("Bz",false);
    fMagFieldCmd->SetUnitCategory("Magnetic flux density");
    fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fUpdateCmd = new G4UIcmdWithoutParameter("/groot/det/update",this);
    fUpdateCmd->SetGuidance("Update calorimeter geometry.");
    fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
    fUpdateCmd->AvailableForStates(G4State_Idle);

    fFilesPathCmd = new G4UIcmdWithAString("/groot/det/FilesPath",this);
    fFilesPathCmd->SetGuidance("Select the path where the detector is defined");
    fFilesPathCmd->SetParameterName("choice",false);
    fFilesPathCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fBackingPosCmd = new G4UIcmdWithAString("/groot/det/setBackingPos",this);
    fBackingPosCmd->SetGuidance("Select the position where the backing is linked with the detector");
    fBackingPosCmd->SetParameterName("choice",false);
    fBackingPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fTargetPhiCmd = new G4UIcmdWithADouble("/groot/det/setTargetPhi",this);
    fTargetPhiCmd->SetGuidance("Set the angle of the detector");
    fTargetPhiCmd->SetParameterName("phi",false);
    fTargetPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fReactionTypeCmd = new G4UIcmdWithAnInteger("/groot/det/setReactionType",this);
    fReactionTypeCmd->SetGuidance("Set the type of reaction.");
    fReactionTypeCmd->SetParameterName("No",false);
    fReactionTypeCmd->SetRange("No>=0");
    fReactionTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
    delete fDetMaterCmd;
    delete fBackingMaterCmd;
    delete fDetThicknessCmd;
    delete fBackingThicknessCmd;
    delete fMagFieldCmd;
    delete fUpdateCmd;
    delete fDetDir;
    delete fGrootDir;
    delete fFilesPathCmd;
    delete fBackingPosCmd;
    delete fTargetPhiCmd;
    delete fReactionTypeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

    if( command == fDetMaterCmd )
    {
        fDetector->SetDetMaterial(newValue);
    }

    if( command == fBackingMaterCmd )
    {
        fDetector->SetBackingMaterial(newValue);
    }

    if( command == fDetThicknessCmd )
    {
        fDetector->SetDetThickness(fDetThicknessCmd->GetNewDoubleValue(newValue));
    }

    if( command == fBackingThicknessCmd )
    {
        fDetector->SetBackingThickness(fBackingThicknessCmd->GetNewDoubleValue(newValue));
    }

    if( command == fFilesPathCmd )
    {
        fDetector->SetPath(newValue);
    }

    if( command == fBackingPosCmd )
    {
        fDetector->SetBackingPos(newValue);
    }

    if( command == fTargetPhiCmd )
    {
        fDetector->SetTargetPhi(fTargetPhiCmd->GetNewDoubleValue(newValue));
    }

    if( command == fReactionTypeCmd )
    {
        fDetector->SetReactionType(fReactionTypeCmd->GetNewIntValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
