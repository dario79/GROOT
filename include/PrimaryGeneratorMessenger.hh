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
/// \file electromagnetic/TestEm1/include/PrimaryGeneratorMessenger.hh
/// \brief Definition of the PrimaryGeneratorMessenger class
//
// $Id: PrimaryGeneratorMessenger.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
   ~PrimaryGeneratorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    PrimaryGeneratorAction* Action;

    G4UIdirectory*             fGunDir;
    G4UIcmdWithoutParameter*   fDefaultCmd;
    G4UIdirectory*             fGrootDir;
    G4UIcmdWithADoubleAndUnit* fEnergyCmd;
    G4UIcmdWithADoubleAndUnit* frbeamCmd;
    G4UIcmdWithAnInteger*      fP1ACmd;
    G4UIcmdWithAnInteger*      fP1ZCmd;
    G4UIcmdWithAnInteger*      fP2ACmd;
    G4UIcmdWithAnInteger*      fP2ZCmd;
    G4UIcmdWithAnInteger*      fP3ACmd;
    G4UIcmdWithAnInteger*      fP3ZCmd;
    G4UIcmdWithADouble*        fmassP1Cmd;
    G4UIcmdWithADouble*        fmassP2Cmd;
    G4UIcmdWithADouble*        fmassP3Cmd;
    G4UIcmdWithADouble*        fmassProjectileCmd;
    G4UIcmdWithADouble*        fmassTargetCmd;
    G4UIcmdWithADoubleAndUnit*        fStateCmd;
    G4UIcmdWithADoubleAndUnit* fRndmCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

