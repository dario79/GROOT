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
/// \file hadronic/Hadr03/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 73011 2013-08-15 08:48:30Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "Run.hh"
#include "HistoManager.hh"
#include "SteppingMessenger.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct, EventAction* eventAction, PrimaryGeneratorAction* prim)
: G4UserSteppingAction(),
  fRunAction(RuAct),
  fEventAction(eventAction),
  fPrimary(prim)
{
    fSteppingMessenger = new SteppingMessenger(this);
    outputType = "both";
    track = "on";
    outputName = "GROOT";
    oldPName = "NONE";
    id = 0;
    oldEnergy = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::SelectOutputType(G4String type)
{
    outputType = type;
    fPrimary->SelectOutputType(type);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::SelectTrack(G4String t)
{
    track = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::SelectOutputName(G4String name)
{
    outputName = name;
    fPrimary->SelectOutputName(name);
    fRunAction->SelectOutputName(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

    G4String TheDeT = "NOT_IMPLEMENTED"; // TO DO : DEFINE STRING FOR ACTIVE VOLUMES
    Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    G4StepPoint* endPoint = aStep->GetPostStepPoint();
    G4VProcess* process = const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());

    std::ofstream outf;
    outf.open(outputName+".txt", std::ios::out | std::ios::app );
    G4EventManager* yo = G4EventManager::GetEventManager();

    // count processes
    //
    fRunAction->CountProcesses(process);

    // check that an real interaction occured (eg. not a transportation)
    G4StepStatus stepStatus = endPoint->GetStepStatus();
    G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
    if (transmit) return;

    //real processes : sum track length
    //
    G4double stepLength = aStep->GetStepLength();
    fRunAction->SumTrack(stepLength);

    //energy-momentum balance initialisation
    //
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4double Q             = - prePoint->GetKineticEnergy();
    G4ThreeVector Pbalance = - prePoint->GetMomentum();

    //initialisation of the nuclear channel identification
    //

    G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
    G4String partName = particle->GetParticleName();
    G4String nuclearChannel = partName;
    G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
    const G4Isotope* target = NULL;
    if (hproc) target = hproc->GetTargetIsotope();
    G4String targetName = "XXXX";  
    if (target) targetName = target->GetName();
    nuclearChannel += " + " + targetName + " --> ";
    if (targetName == "XXXX") run->SetTargetXXX(true);
    //scattered primary particle (if any)

//    G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    G4int ih = 3;
    G4int ih1 = ih+10;
    G4int ih2 = ih+20;
    G4int ih3 = ih+30;

    G4double posX = endPoint-> GetPosition().x() / CLHEP::cm ;
    G4double posY = endPoint-> GetPosition().y() / CLHEP::cm ;
    G4double posZ = endPoint-> GetPosition().z() / CLHEP::cm ;
    G4double momentumX = endPoint-> GetMomentumDirection().x();
    G4double momentumY = endPoint-> GetMomentumDirection().y();
    G4double momentumZ = endPoint-> GetMomentumDirection().z();
    G4double energy = endPoint->GetKineticEnergy();
    G4double posX0 = prePoint-> GetPosition().x() / CLHEP::cm ;
    G4double posY0 = prePoint-> GetPosition().y() / CLHEP::cm ;
    G4double posZ0 = prePoint-> GetPosition().z() / CLHEP::cm ;
    G4double momentumX0 = prePoint-> GetMomentumDirection().x();
    G4double momentumY0 = prePoint-> GetMomentumDirection().y();
    G4double momentumZ0 = prePoint-> GetMomentumDirection().z();
    G4double energy0 = prePoint->GetKineticEnergy();

    G4int trackid = aStep->GetTrack()->GetTrackID();
    G4int parentid = aStep->GetTrack()->GetParentID();
    G4double EdepStep = aStep->GetTotalEnergyDeposit();

    G4double t_i = prePoint->GetGlobalTime();
    G4double t_f = endPoint->GetGlobalTime();

    //  G4Step--> GetTotalEnergyDeposit()
//    analysis->FillH1(ih,energy);
//    analysis->FillH1(ih1,momentumX);
//    analysis->FillH1(ih2,momentumY);
//    analysis->FillH1(ih3,momentumZ);

    auto analysisManager = G4AnalysisManager::Instance();
    G4String volumeName =  prePoint->GetPhysicalVolume()->GetName();

    G4ThreeVector momentum = endPoint->GetMomentum();
    Q        += energy;
    Pbalance += momentum;

    nuclearChannel += partName + " + ";

    if(outputType != "txt" && track == "on")
    {
        analysisManager->FillNtupleDColumn(1, 0, posX0);
        analysisManager->FillNtupleDColumn(1, 1, posY0);
        analysisManager->FillNtupleDColumn(1, 2, posZ0);
        analysisManager->FillNtupleDColumn(1, 3, momentumX0);
        analysisManager->FillNtupleDColumn(1, 4, momentumY0);
        analysisManager->FillNtupleDColumn(1, 5, momentumZ0);
        analysisManager->FillNtupleDColumn(1, 6, energy0);
        analysisManager->FillNtupleDColumn(1, 7, posX);
        analysisManager->FillNtupleDColumn(1, 8, posY);
        analysisManager->FillNtupleDColumn(1, 9, posZ);
        analysisManager->FillNtupleDColumn(1, 10, momentumX);
        analysisManager->FillNtupleDColumn(1, 11, momentumY);
        analysisManager->FillNtupleDColumn(1, 12, momentumZ);
        analysisManager->FillNtupleDColumn(1, 13, energy);
        analysisManager->FillNtupleSColumn(1, 14, volumeName);
        analysisManager->FillNtupleIColumn(1, 15, yo->GetConstCurrentEvent()->GetEventID());
        analysisManager->FillNtupleSColumn(1, 16, partName);
        analysisManager->FillNtupleSColumn(1, 17, process->GetProcessName());
        analysisManager->FillNtupleIColumn(1, 18, trackid);
        analysisManager->FillNtupleIColumn(1, 19, parentid);
        analysisManager->FillNtupleDColumn(1, 20, EdepStep);
        analysisManager->FillNtupleDColumn(1, 21, t_i);
        analysisManager->FillNtupleDColumn(1, 22, t_f);


        analysisManager->AddNtupleRow(1);
    }

     //G4String volumeName =  prePoint->GetPhysicalVolume()->GetName();

    if(outputType != "root" && track == "on")
    {
        outf << " \t"
           << std::setprecision(5) << yo->GetConstCurrentEvent()->GetEventID()<< " \t"
           << partName << " \t"
           << std::setprecision(5) << energy0 << " \t"
           << std::setprecision(5) << posX0 << " \t"
           << std::setprecision(5) << posY0 << " \t"
           << std::setprecision(5) << posZ0 << " \t"
           << std::setprecision(5) << momentumX0 << " \t"
           << std::setprecision(5) << momentumY0 << " \t"
           << std::setprecision(5) << momentumZ0 << " \t"
           << std::setprecision(5) << energy << " \t"
           << std::setprecision(5) << posX << " \t"
           << std::setprecision(5) << posY << " \t"
           << std::setprecision(5) << posZ << " \t"
           << std::setprecision(5) << momentumX << " \t"
           << std::setprecision(5) << momentumY << " \t"
           << std::setprecision(5) << momentumZ << " \t"
           << process->GetProcessName() << " \t"
           << trackid << " \t"
           << parentid << " \t"
           << std::setprecision(5) << EdepStep << " \t"
           << volumeName << " \t" 
           << std::setprecision(5) << t_i << " \t"
           << std::setprecision(5) << t_f << " \t"
	   << G4endl;
    }

    G4double DIST = sqrt(posX*posX+posY*posY+posZ*posZ);
    G4ThreeVector thepos(posX/DIST,posY/DIST,posZ/DIST);

    if(oldPName != partName || oldEnergy != energy0)
    {
        id++;
    }
    oldPName = partName;
    oldEnergy = energy;

    fEventAction->AddEdep(energy0 - energy, thepos, partName, volumeName, outputType, 1, id);


  //secondaries
  //
    const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();

    for (size_t lp=0; lp<(*secondary).size(); lp++)
    {
        particle = (*secondary)[lp]->GetDefinition();
        G4String name   = particle->GetParticleName();
        G4String type   = particle->GetParticleType();
        G4double charge = particle->GetPDGCharge();
        energy = (*secondary)[lp]->GetKineticEnergy();
        momentumX = (*secondary)[lp] -> GetMomentumDirection().x();
        momentumY = (*secondary)[lp] -> GetMomentumDirection().y();
        momentumZ = (*secondary)[lp] -> GetMomentumDirection().z();
        posX = (*secondary)[lp] ->GetPosition().x()/ CLHEP::cm;
        posY = (*secondary)[lp] ->GetPosition().y()/ CLHEP::cm;
        posZ = (*secondary)[lp] ->GetPosition().z()/ CLHEP::cm;
        G4ThreeVector pos(posX/DIST,posY/DIST,posZ/DIST);

        fRunAction->ParticleCount(name,energy);
    //energy spectrum
/*        if (charge > 3.)  {ih = 4;  ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (particle == G4Gamma::Gamma())       {ih = 5; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (particle == G4Neutron::Neutron())   {ih = 6; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (particle == G4Proton::Proton())     {ih = 7; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (particle == G4Deuteron::Deuteron()) {ih = 8; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (particle == G4Alpha::Alpha())       {ih = 9; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (type == "nucleus")                  {ih = 10; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (type == "meson")                    {ih = 11; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        else if (type == "baryon")                   {ih = 12; ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
        analysis->FillH1(ih,energy);
        analysis->FillH1(ih1,momentumX);
        analysis->FillH1(ih2,momentumY);
        analysis->FillH1(ih3,momentumZ);
    //energy-momentum balance
        momentum = (*secondary)[lp]->GetMomentum();
        Q        += energy;
        Pbalance += momentum;
        volumeName = (*secondary)[lp]->GetVolume()->GetName();
*/
        if (0/*volumeName.compare(0,2,TheDeT)==0*/)
        {
            if(outputType != "root" && track == "on")
            {
                outf << " secondary: \t \t"
                    << std::setprecision(5) << yo->GetConstCurrentEvent()->GetEventID()<< " \t"
                    << name << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << (*secondary)[lp]->GetParentID() << " \t"
                    << std::setprecision(5) << energy << " \t"
                    << std::setprecision(5) << posX << " \t"
                    << std::setprecision(5) << posY << " \t"
                    << std::setprecision(5) << posZ << " \t"
                    << std::setprecision(5) << momentumX << " \t"
                    << std::setprecision(5) << momentumY << " \t"
                    << std::setprecision(5) << momentumZ << " \t"
                    << (*secondary)[lp]->GetCreatorProcess()->GetProcessName() << " \t"
                    << volumeName << G4endl;
            }

            if(outputType != "txt" && track == "on")
            {
                analysisManager = G4AnalysisManager::Instance();
                analysisManager->FillNtupleDColumn(1, 0, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 1, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 2, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 3, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 4, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 5, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 6, (*secondary)[lp]->GetParentID());
                analysisManager->FillNtupleDColumn(1, 7, posX);
                analysisManager->FillNtupleDColumn(1, 8, posY);
                analysisManager->FillNtupleDColumn(1, 9, posZ);
                analysisManager->FillNtupleDColumn(1, 10, momentumX);
                analysisManager->FillNtupleDColumn(1, 11, momentumY);
                analysisManager->FillNtupleDColumn(1, 12, momentumZ);
                analysisManager->FillNtupleDColumn(1, 13, energy);
                analysisManager->FillNtupleSColumn(1, 14, volumeName);
                analysisManager->FillNtupleIColumn(1, 15, yo->GetConstCurrentEvent()->GetEventID());
                analysisManager->FillNtupleSColumn(1, 16, partName);
                analysisManager->AddNtupleRow(1);
            }

            fEventAction->AddEdep(energy, pos, name, volumeName, outputType, 0, 0);
        }
    //particle flag
     fParticleFlag[particle]++;
    }
    outf.close();

  //energy-momentum balance
    G4double Pbal = Pbalance.mag();
    fRunAction->Balance(Pbal);
    ih = 2; // ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;
//    analysis->FillH1(ih,Q);

    ih = 1;  // ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;
 //   analysis->FillH1(ih,Q);


//    analysis->FillH1(0,0);

  // nuclear channel
    const G4int kMax = 16;
    const G4String conver[] = {"0","","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ",
                             "10 ","11 ","12 ","13 ","14 ","15 ","16 "};
    std::map<G4ParticleDefinition*,G4int>::iterator ip;
    for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++)
    {
        particle = ip->first;
        G4String name = particle->GetParticleName();
        G4int nb = ip->second;
        if (nb > kMax) nb = kMax;
        G4String Nb = conver[nb];
        if (particle == G4Gamma::Gamma())
        {
             fRunAction->CountGamma(nb);
             Nb = "N ";
        }
        if (ip != fParticleFlag.begin()) nuclearChannel += " + ";

        nuclearChannel += Nb + name;
    }

  ///G4cout << "\n nuclear channel: " << nuclearChannel << G4endl;
    fRunAction->CountNuclearChannel(nuclearChannel, Q);

    fParticleFlag.clear();
  // kill event after first interaction

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


