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
/// \file hadronic/Hadr03/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
// $Id: SteppingAction.hh 73011 2013-08-15 08:48:30Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <map>

class G4ParticleDefinition;
class RunAction;
class EventAction;
class SteppingMessenger;
class PrimaryGeneratorAction;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(RunAction*, EventAction*, PrimaryGeneratorAction* );
   ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    G4int sw, id;
    G4String PName, outputType, track, outputName, oldPName;
    G4double oldEnergy;
    void SelectOutputType(G4String);
    void SelectOutputName(G4String);
    void SelectTrack(G4String);

    
  private:
    RunAction*              fRunAction;
    EventAction*            fEventAction;
    PrimaryGeneratorAction*    fPrimary;
    SteppingMessenger*      fSteppingMessenger;
    std::map<G4ParticleDefinition*,G4int> fParticleFlag;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
