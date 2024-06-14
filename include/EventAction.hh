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
/// \file hadronic/Hadr03/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 66241 2012-12-13 18:34:42Z gunter $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <vector>

class EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    
    void SetPrintModulo(G4int val) {fPrintModulo = val;};
    void SetDrawFlag(G4String val) {fDrawFlag = val;};
    void AddEdep(G4double, G4ThreeVector, G4String, G4String, G4String, G4int, G4int);

  private:
    G4int                 fPrintModulo;                    
    EventActionMessenger* fEventMessenger;
    G4String               fDrawFlag;
    G4double theta,phi;//,e_rec;
    G4String outputType;
    std::vector <G4int> Type;
    std::vector <G4double> Energy;
    std::vector <G4double> Theta;
    std::vector <G4double> Phi;
    std::vector <G4String> PName;
    std::vector <G4String> VName;
    std::vector <G4String> deleteS;
    std::vector <G4int> deleteI;
    std::vector <G4double> deleteD;
    std::vector <G4int> ID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
