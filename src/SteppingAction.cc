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
#include "RunAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct)
: G4UserSteppingAction(),fRunAction(RuAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4String TheDeT = "YY";
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4VProcess* process = const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());

   std::ofstream outf;
   outf.open("GROOT.txt", std::ios::out | std::ios::app );
   G4EventManager* yo = G4EventManager::GetEventManager();

  // count processes
  // 
//  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
//  const G4VProcess* process   = endPoint->GetProcessDefinedStep();
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
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
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

/*  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
  G4String partName = particle->GetParticleName();
  G4String nuclearChannel = partName;
  G4HadronicProcess* hproc = (G4HadronicProcess*) process;
  const G4Isotope* target = hproc->GetTargetIsotope();
  G4String targetName; // = "XXXX";  
/*  if (target) { 
          targetName = target->GetName();
  }
  else {  
          targetName = "XXXX";
       }
  nuclearChannel += " + " + targetName + " --> ";
  if (targetName == "XXXX")  run->SetTargetXXX(true);  //G4cout<<"Could be a problem.."<<G4endl; //
  */
  //scattered primary particle (if any)
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int ih = 3;
  G4int ih1 = ih+10;
  G4int ih2 = ih+20;
  G4int ih3 = ih+30;
  if (1/* aStep->GetTrack()->GetTrackStatus() == fAlive*/) {
   G4double posX = aStep-> GetTrack()-> GetPosition().x() / CLHEP::cm ;
   G4double posY = aStep-> GetTrack()-> GetPosition().y() / CLHEP::cm ;
   G4double posZ = aStep-> GetTrack()-> GetPosition().z() / CLHEP::cm ;
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
//  G4Step--> GetTotalEnergyDeposit()
    analysis->FillH1(ih,energy);
    analysis->FillH1(ih1,momentumX);
    analysis->FillH1(ih2,momentumY);
    analysis->FillH1(ih3,momentumZ);
    //
    G4ThreeVector momentum = endPoint->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    //
    nuclearChannel += partName + " + ";

     G4String volumeName =  prePoint->GetPhysicalVolume()->GetName();
//    if( (volumeName=="CarbonLayerPhys" || volumeName=="MagnesiumLayerPhys") && (1.-momentumZ>0.0001) ) {
//     if ( volumeName=="VoidOrbPhys" && posZ>0. ){
     if (volumeName.compare(0,2,TheDeT)==0){
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
               << volumeName << G4endl;
     }
  }  
  
  //secondaries
  //
  const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();  
//   G4EventManager* yo = G4EventManager::GetEventManager();

  
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    particle = (*secondary)[lp]->GetDefinition(); 
    G4String name   = particle->GetParticleName();
    G4String type   = particle->GetParticleType();      
    G4double charge = particle->GetPDGCharge();
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    G4double momentumX = (*secondary)[lp] -> GetMomentumDirection().x();
    G4double momentumY = (*secondary)[lp] -> GetMomentumDirection().y();
    G4double momentumZ = (*secondary)[lp] -> GetMomentumDirection().z();
    G4double posX = (*secondary)[lp] ->GetPosition().x()/ CLHEP::cm;
    G4double posY = (*secondary)[lp] ->GetPosition().y()/ CLHEP::cm;
    G4double posZ = (*secondary)[lp] ->GetPosition().z()/ CLHEP::cm;
    

/*    if (particle == G4Alpha::Alpha()) G4cout<< momentumX << "\t"<< momentumY << "\t"<< momentumZ <<G4endl;
  	  G4ThreeVector momentumd = (*secondary)[lp] -> GetMomentumDirection();
    momentumX = momentumd.x();
    momentumY = momentumd.y();
    momentumZ = momentumd.z();
    if (particle == G4Alpha::Alpha()) G4cout<< momentumX << "\t"<< momentumY << "\t"<< momentumZ <<G4endl;
*/    fRunAction->ParticleCount(name,energy);
    //energy spectrum
    if (charge > 3.)  {ih = 4;  ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;}
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
    G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();
    Q        += energy;
    Pbalance += momentum;
    G4String volumeName = (*secondary)[lp]->GetVolume()->GetName();

//    if (particle != G4Gamma::Gamma()){  // NOT A GAMMA
/*    if (particle == G4Alpha::Alpha() && (*secondary).size()<5){  // ALPHA AND nucleus without outgoing gammas
//      if(G4String nextvol = (*secondary)[lp]->GetNextVolume()->GetName()) G4cout<<" jj "<<G4endl;
          yo->KeepTheCurrentEvent();
 	  outf<<" Event "<<yo->GetConstCurrentEvent()->GetEventID()<<" ... has a "<<name<<" , created by "
              <<(*secondary)[lp]->GetCreatorProcess()->GetProcessName()<<" in "<<(*secondary)[lp]->GetVolume()->GetName()<<" going "
//              <<(*secondary)[lp]->GetNextVolume()->GetName()
              <<" with other "<<(*secondary).size() -1<<" object(s)."<<G4endl;
//           exit(0);
//          delete yo;
//          G4cout<<" Event is ... "<<particle->GetG4code()->GetParticleName()<<G4endl;
//            evt.ToBeKept();
    }
    else {*/
//       if(volumeName=="physiQQ3" || volumeName=="physiSX3"){
//       if( (volumeName=="CarbonLayerPhys" || volumeName=="MagnesiumLayerPhys") && (1.-momentumZ>0.0001) ) {
//       if ( volumeName=="physiVoidOrb" && (momentumZ+1.>0.01) ){
//         if(momentumZ>0){     // particles going in the same direction of the beam -->

       
//        if( (volumeName=="physiQQ3" || volumeName=="physiSX3") 
/*        if ( volumeName!="WorldPhys"){ 
            G4String nextvolume = (*secondary)[lp]->GetNextVolume()->GetName();
          if(nextvolume=="physiQQ3" || nextvolume=="physiSX3") {
*/          
//            G4cout<<"DAMN IT!!!!!!!!"<<G4endl;
//     if ( volumeName=="VoidOrbPhys" && posZ>0. ){
      if (volumeName.compare(0,2,TheDeT)==0){
            outf << " \t \t"
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
//          }
      }
    //particle flag
    fParticleFlag[particle]++;
  }
  outf.close();  

  //energy-momentum balance
  G4double Pbal = Pbalance.mag();
  fRunAction->Balance(Pbal);
  ih = 2; // ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;
  analysis->FillH1(ih,Q);
//  analysis->FillH1(ih1,Pbal);
//  analysis->FillH1(ih2,Pbal);
//  analysis->FillH1(ih3,Pbal);
  ih = 1;  // ih1 = ih+10; ih2 = ih+20; ih3 = ih+30;
  analysis->FillH1(ih,Q);  
//  analysis->FillH1(ih1,Q);
//  analysis->FillH1(ih2,Q);
//  analysis->FillH1(ih3,Q);

  analysis->FillH1(0,0);

  // nuclear channel
  const G4int kMax = 16;  
  const G4String conver[] = {"0","","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ",
                             "10 ","11 ","12 ","13 ","14 ","15 ","16 "};
  std::map<G4ParticleDefinition*,G4int>::iterator ip;               
  for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++) {
    particle = ip->first;
    G4String name = particle->GetParticleName();      
    G4int nb = ip->second;
    if (nb > kMax) nb = kMax;   
    G4String Nb = conver[nb];    
    if (particle == G4Gamma::Gamma()) {
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
  //
//  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


