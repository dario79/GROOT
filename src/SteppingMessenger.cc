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
/// \file electromagnetic/Groot1/src/SteppingMessenger.cc
/// \brief Implementation of the SteppingMessenger class
//
// $Id: SteppingMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingMessenger.hh"
#include "G4UIdirectory.hh"
#include "SteppingAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingMessenger::SteppingMessenger( SteppingAction* Step)
:G4UImessenger(),fStepping(Step),
 fStepDir(0),
 fGrootDir(0),
 fOutputTypeCmd(0),
 fOutputNameCmd(0),
 fTrackCmd(0)
{
    fGrootDir = new G4UIdirectory("/groot/");
    fGrootDir->SetGuidance("commands specific to this example");

    fStepDir = new G4UIdirectory("/groot/analysis/");
    fStepDir->SetGuidance("gun control");

    fOutputTypeCmd = new G4UIcmdWithAString("/groot/analysis/outputType",this);
    fOutputTypeCmd->SetGuidance("Select the type of the output");
    fOutputTypeCmd->SetParameterName("choice",false);
    fOutputTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fOutputNameCmd = new G4UIcmdWithAString("/groot/analysis/outputName",this);
    fOutputNameCmd->SetGuidance("Select the name of the output");
    fOutputNameCmd->SetParameterName("choice",false);
    fOutputNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    fTrackCmd = new G4UIcmdWithAString("/groot/track",this);
    fTrackCmd->SetGuidance("Select if you want the track");
    fTrackCmd->SetParameterName("choice",false);
    fTrackCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingMessenger::~SteppingMessenger()
{
    delete fStepDir;
    delete fOutputTypeCmd;
    delete fOutputNameCmd;
    delete fTrackCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{

    if(command == fOutputTypeCmd)
    {
        fStepping->SelectOutputType(newValue);
    }

    if(command == fOutputNameCmd)
    {
        fStepping->SelectOutputName(newValue);
    }

    if(command == fTrackCmd)
    {
        fStepping->SelectTrack(newValue);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

