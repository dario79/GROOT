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

#include "EventActionMessenger.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleTypes.hh"
#include "G4VUserEventInformation.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

void EventAction::EndOfEventAction(const G4Event* evt)
{


/*if () {
   G4EventManager* em = G4EventManager::GetEventManager();
   em->KeepTheCurrentEvent();
}*/

/* G4VUserEventInformation* eventInformation = (G4VUserEventInformation*) evt->GetUserInformation();
 std::ofstream outf;
 outf.open("camurria.txt", std::ios::out | std::ios::app );
 G4EventManager* yo = G4EventManager::GetEventManager();
 G4TrajectoryContainer* trajectoryContainer=evt->GetTrajectoryContainer();
 
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
 
  // extract the trajectories and draw them
     if (G4VVisManager::GetConcreteInstance()){
         for (G4int i=0; i<n_trajectories; i++){
             G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
             if(trj->GetParticleName()!="gamma" || trj->GetPDGEncoding() != 22){
//                   trj->SetForceDrawTrajectory(fForcedrawphotons);
//                   trj->SetForceNoDrawTrajectory(fForcenophotons);                  
                   yo->KeepTheCurrentEvent();
                   outf<<" Event is ... "<<trj->GetParticleName()<<" , code is "<<trj->GetPDGEncoding()<<G4endl;
//        delete yo;
                   G4cout<<" Event is ... "<<trj->GetParticleName()<<" , code is "<<trj->GetPDGEncoding()<<G4endl;
//          evt.ToBeKept();
             trj->DrawTrajectory();
             }
         }
     }

*/

/*       G4PrimaryParticle* particle = evt->GetPrimaryVertex()->GetPrimary();
//   G4String partName = particle->GetParticleName();
//       if (partname=='alpha' )
 std::ofstream outf;
 outf.open("camurria.txt", std::ios::out | std::ios::app );
 G4EventManager* yo = G4EventManager::GetEventManager();


        if(particle->GetPDGcode() != 22){  // NOT A GAMMA
        yo->KeepTheCurrentEvent();
 	outf<<" Event is ... "<<particle->GetG4code()->GetParticleName()<<" , code is "<<particle->GetPDGcode()<<G4endl;
        G4cout<<" Event is ... "<<particle->GetG4code()->GetParticleName()<<G4endl;
//          evt.ToBeKept();
        }
//        delete yo;
*/
//  outf.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


