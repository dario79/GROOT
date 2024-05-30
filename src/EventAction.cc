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
/// \file hadronic/Hadr03/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 70759 2013-06-05 12:26:43Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "HistoManager.hh"
#include "EventActionMessenger.hh"

#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTypes.hh"
#include "G4VUserEventInformation.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
 fPrintModulo(100000),fEventMessenger(0)
{
    fEventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
    delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
    G4int evtNb = evt->GetEventID();

 //printing survey
    if (evtNb%fPrintModulo == 0)
        G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
   Energy.clear();
   Type.clear();
   PName.clear();
   VName.clear();
   Theta.clear();
   Phi.clear();
   ID.clear();

   swap(Energy, deleteD);
   deleteD = Energy;
   swap(Theta, deleteD);
   deleteD = Theta;
   swap(Phi, deleteD);
   deleteD = Phi;
   swap(ID, deleteI);
   deleteI = ID;
   swap(Type, deleteI);
   deleteI = Type;
   swap(PName, deleteS);
   deleteS = PName;
   swap(VName, deleteS);
   deleteS = VName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
    //std::ofstream outf;
    //outf.open("test.txt", std::ios::out | std::ios::app );

    if(outputType != "txt")
    {
        for(int i = 0; i < (int)Energy.size(); i++)
        {
            auto analysisManager = G4AnalysisManager::Instance();
            analysisManager->FillNtupleSColumn(2, 0, PName[i]);
            analysisManager->FillNtupleDColumn(2, 1, Energy[i]);
            //analysisManager->FillNtupleDColumn(2, 2, e_rec);
            analysisManager->FillNtupleDColumn(2, 2, Theta[i]);
            analysisManager->FillNtupleDColumn(2, 3, Phi[i]);
            analysisManager->FillNtupleSColumn(2, 4, VName[i]);
            analysisManager->FillNtupleIColumn(2, 5, evt->GetEventID());
            analysisManager->AddNtupleRow(2);

            //outf<<evt->GetEventID()<<' '<<PName[i]<<' '<<Energy[i]<<' '<<Theta[i]<<' '<<Phi[i]<<' '<<VName[i]<<'\n';
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4double edep, G4ThreeVector pos, G4String pname, G4String vname, G4String t, G4int type, G4int id)
{

    int sw = 0;

    for(int i = 0; i < (int)Energy.size(); i++)
    {
        if(Type[i] == 1 && ID[i] == id)
        {
            Energy[i] += edep;
            sw = 1;
        }
    }
    if(sw == 0)
    {
        PName.push_back(pname);
        Energy.push_back(edep);
        VName.push_back(vname);
        Type.push_back(type);
        theta = acos(pos.z());
        if(theta<=0) theta += twopi;
        phi =  atan2(pos.y(),pos.x());
        if(phi<=0) phi += twopi;
        Theta.push_back(theta);
        Phi.push_back(phi);
        ID.push_back(id);
    }
    outputType = t;
    //e_rec=0;
   // if (PName!="e-") e_rec = G4RandGauss::shoot(fEdepH.GetValue(),fEdepH.GetValue()/100.);
// acos z/dist 
//    theta2 = atan2(pos.y(),pos.z());
//    phi2 =  atan2(pos.x(),pos.z());
}


