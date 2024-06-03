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
/// \file electromagnetic/TestEm1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// $Id: PrimaryGeneratorAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction* PrimaryGeneratorAction::instance = 0;
DetectorConstruction* det;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* DC)
:G4VUserPrimaryGeneratorAction(),
 fParticleGun(0),
 fRndmBeam(0),    
 fGunMessenger(0)                               
{
  // Definition of the General particle Source
//  fParticleGun = new G4GeneralParticleSource();

  fParticleGun  = new G4ParticleGun(1);
  SetDefaultKinematic(2);  // 2 = 2body, 3 = 3body  else if standard
    
  fRndmBeam = 0.;
  //create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetDefaultKinematic(G4int front=0)
{
  fuffa = front;
  targetphi = 0.*deg;// -35.*deg;
  g_Energy = 40.23*MeV; // MeV
//  counter = 1;
  SetCounter(1);
/*  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(11*MeV);
  G4double position = -1.*cm;
  if (front) position = -0.5*(fDetector->GetSize());
  fParticleGun->SetParticlePosition(G4ThreeVector(position,0.*cm,0.*cm));*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


/*G4SingleParticleSource* PrimaryGeneratorAction::GetParticleGun(){
   return fParticleGun->GetCurrentSource();
}
*/

/*void PrimaryGeneratorAction::GetGenMode(G4double sEne, G4ThreeVector sPosition, G4ThreeVector sDirection){
   energiaIniziale = sEne;
   oldPosition = sPosition;
   oldDirection = sDirection;
   counter++;
}*/

void PrimaryGeneratorAction::SetCounter(G4int VALUE=0){
  counter = VALUE;
}

void PrimaryGeneratorAction::AddCounter(G4int VALUE=0){
  counter += VALUE;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//   GetGenMode();
   std::ofstream outf;
   outf.open("Eli.txt", std::ios::out | std::ios::app );
//   targetphi = 35.*deg;

if(fuffa==2){
//  unita' di massa (rif. 12C)
   G4double root_u = 931.49432/1000.; //
   G4int root_out_particles = 2;     // n. output particles
////////// scattering N15,a up to 40.23*MeV /////////////////////////////
   G4double     root_mp=15.0001;    // massa proiettile             //
   G4double     root_mt=4.002603;   // massa target                 //
   G4double     root_m1=root_mp;    // massa prima part. uscente    //
   G4double     root_m2=root_mt;    // massa seconda part. uscente  //

//  (Momentum, Energy units are Gev/C, GeV)
   G4double root_masses[2] = {root_u*root_mp,root_u*root_mt};
   TGenPhaseSpace TheRootEvent;
   TRandom3* ran = new TRandom3(0);
   ran->SetSeed(0);

//        std::normal_distribution<double> distribution(/*mean=*/11.7, /*stddev=*/0.015);
//        std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015);
//        Double_t Gaus(mean,sigma);
//        double randomNumber = distribution(generator);

  //************************************************************ ****//
 ////////////////////////// EXTRACT ENERGY,Z HERE /////////////////// 
//*****************************************************************//
  if((counter%2)==0) {
// the 15N has been propagated, hopefully
//   counter++;
//  G4cout<<"COUNTER IS: "<<counter<<G4endl;
   AddCounter(1);
   // energia proiettile
//     G4double randomNumber = G4UniformRand()*(g_energy - 11.); // estimated energy loss in havar in MeV
  std::ifstream inshit;
  inshit.open("thatshit", std::ios::in);
// new lines will be skipped unless we stop it from happening:    
  inshit.unsetf(std::ios_base::skipws);
  inshit.clear(); // clear bad state after eof
  inshit.seekg(0);
  unsigned line_count = std::count(std::istream_iterator<char>(inshit), std::istream_iterator<char>(), '\n');
//  G4cout<<"line_count IS: "<<line_count<<G4endl;
  inshit.setf(std::ios_base::skipws);
  inshit.clear(); // clear bad state after eof
  inshit.seekg(0);

  unsigned randomline = std::rand() % line_count + 1;
//  G4cout<<"randomline IS: "<<randomline<<G4endl;
 for(unsigned yeah = 0; yeah < randomline+1; yeah++){
       inshit>>energiaIniziale>>oldposX>>oldposY>>oldposZ>>oldDirX>>oldDirY>>oldDirZ;
 }
   oldposX /= mm;
   oldposY /= mm;
   oldposZ /= mm;
   inshit.close();
// energia proiettile in GeV
   G4double root_ep=(energiaIniziale/1000.); // - root_u*root_m1);
// impulso proiettile in GeV/c
   G4double  root_pp = sqrt(2.*root_u*root_mp*root_ep); // root_ep/root_mp;
//   G4cout<<root_ep<<"\t"<<root_pp<<G4endl;
   TLorentzVector root_target(0.0, 0.0, 0.0, root_u*root_mt);
   TLorentzVector root_beam(root_pp*oldDirX, root_pp*oldDirY, root_pp*oldDirZ, root_ep + root_u*root_mp);
   TLorentzVector rootW = root_beam + root_target;
   G4int seed2 = (int) G4UniformRand()*10000 +1;
   gRandom->SetSeed(seed2);
   G4bool isThatOK = TheRootEvent.SetDecay(rootW, root_out_particles, root_masses);
//   G4cout<<"It's ok. "<<isThatOK<<G4endl;
 for(unsigned looping = 0; looping<1000000000; looping++){
    G4double root_weight = TheRootEvent.Generate();

    TLorentzVector *root_p1 = TheRootEvent.GetDecay(0);
    TLorentzVector *root_p2 = TheRootEvent.GetDecay(1);

    TVector3 root_pp1=root_p1->Vect();
    TVector3 root_pp2=root_p2->Vect();  
    TVector3 root_dir1=root_pp1.Unit();
    TVector3 root_dir2=root_pp2.Unit();

    G4ParticleDefinition* out_particle1 = G4IonTable::GetIonTable()->GetIon(7, 15, 0.);
    G4ParticleDefinition* out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("alpha");

    G4double E1 = 1000.*(root_p1->E() - root_u*root_m1);
    G4double E2 = 1000.*(root_p2->E() - root_u*root_m2);
    G4double theta2 = 180.*root_dir2.Theta()/pi; 
    G4double phi2 = 180.*root_dir2.Phi()/pi;
    G4double ThetaBin = 3.;
    G4double ThetaHit1 = 

/*    outf<<"THETA: "<<theta2<<"; PHI: "<<phi2<<G4endl;
    outf<<"X: "<<oldposX*mm<<"; Y: "<<oldposY*mm<<"; Z: "<<oldposZ*mm<<G4endl;
*/
 G4bool Hit_it = false;
 if(TMath::Abs(phi2)<0.1){ 
//      G4double tsquared = theta2*theta2;
      if(TMath::Abs(theta2-ThetaHit1)<ThetaBin) Hit_it = true;
      else if(TMath::Abs(theta2-ThetaHit2)<ThetaBin) Hit_it = true;
      else if(TMath::Abs(theta2-ThetaHit3)<ThetaBin) Hit_it = true;
      else if(TMath::Abs(theta2-ThetaHit4)<ThetaBin) Hit_it = true;
      else if(TMath::Abs(theta2-ThetaHit5)<ThetaBin) Hit_it = true;
 }
 if (root_weight>0. && Hit_it /*&& root_dir2(2)>0.1 && (((theta2*theta2)<=2500.) && (phi2*phi2)<0.01*/)){
  G4double t0 = (looping+1)/*(1.*s)*/;
  G4PrimaryVertex* vertexA = new G4PrimaryVertex(G4ThreeVector(oldposX*mm,oldposY*mm,oldposZ*mm), t0*s);
//  G4PrimaryParticle* particle1 = new G4PrimaryParticle(out_particle1);
  G4PrimaryParticle* particle2 = new G4PrimaryParticle(out_particle2);
//  particle1->SetMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));
//  particle1->SetKineticEnergy(E1*MeV);
  particle2->SetMomentumDirection(G4ThreeVector(root_dir2(0),root_dir2(1),root_dir2(2)));
  particle2->SetKineticEnergy(E2*MeV);
//  vertexA->SetPrimary(particle1);
  vertexA->SetPrimary(particle2);
  anEvent->AddPrimaryVertex(vertexA);
/*    fParticleGun->SetParticlePosition(G4ThreeVector(oldposX*mm,oldposY*mm,oldposZ*mm));
    fParticleGun->SetParticleDefinition(out_particle1);
    fParticleGun->SetParticleEnergy(E1*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));
//    G4cout<<root_p1->E()<<"\t"<<root_dir1.x()<<"\t"<<root_dir1.y()<<"\t"<<root_dir1.z()<<"\t"<<out_particle1->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);*/
/*         outf  << "Scattered particle 1 (not propagated): "
               << out_particle1->GetParticleName() << " \t" 
               << std::setprecision(5) << E1*MeV<< " \t"
               << std::setprecision(5) << oldposX*mm << " \t"
               << std::setprecision(5) << oldposY*mm << " \t"
               << std::setprecision(5) << oldposZ*mm << " \t"
               << std::setprecision(5) << root_dir1(0) << " \t"
               << std::setprecision(5) << root_dir1(1) << " \t"
               << std::setprecision(5) << root_dir1(2) << " \t"
               << G4endl;*/
         outf  << "Scattered particle 2: "
//              <<"THETA: "<<theta2<<"; PHI: "<<phi2<< " \t" 
               << out_particle2->GetParticleName() << " \t" 
               << std::setprecision(5) << E2*MeV<< " \t"
               << std::setprecision(5) << oldposX*mm << " \t"
               << std::setprecision(5) << oldposY*mm << " \t"
               << std::setprecision(5) << oldposZ*mm << " \t"
               << std::setprecision(5) << root_dir2(0) << " \t"
               << std::setprecision(5) << root_dir2(1) << " \t"
               << std::setprecision(5) << root_dir2(2) << " \t"
//               << std::setprecision(5) << t0<< " \t"
               << G4endl;
 }
 else {
/*       G4cout<<"Just a warning for root_weight <= 0. W = "<<root_weight<<G4endl;
       G4cout<<"counter:"<<counter<<"\n"<<energiaIniziale<<"\t"<<oldposX<<"\t"<<oldposY<<"\t"<<oldposZ<<"\t"<<oldDirX<<"\t"<<oldDirY<<"\t"<<oldDirZ<<G4endl;
       G4cout<<out_particle1->GetParticleName()<<"\t"<<E1*MeV<<"\t"<<oldposX<<"\t"<<oldposY<<"\t"<<oldposZ<<"\t"<<root_dir1(0)<<"\t"<<root_dir1(1)<<"\t"<<root_dir1(2)<<G4endl;
       G4cout<<out_particle2->GetParticleName()<<"\t"<<E2*MeV<<"\t"<<oldposX<<"\t"<<oldposY<<"\t"<<oldposZ<<"\t"<<root_dir2(0)<<"\t"<<root_dir2(1)<<"\t"<<root_dir2(2)<<G4endl;*/
 } // add "abort event" ?
//  }
// }
 }
 }
 else {
// 15N to calculate the energy loss 
//  counter++;
   AddCounter(1);
  std::ofstream outshit;
  outshit.open("thatshit", std::ios::out);
  outshit.close();
  G4double position = -1000.00215;
//  std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015*MeV);
  G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(7, 15, 0.);
  fParticleGun->SetParticleDefinition(particle);
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(g_Energy);
//  fParticleGun->SetParticlePosition(oldPosition);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,position*mm));
//  fParticleGun->SetParticleMomentumDirection(oldDirection);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << particle->GetParticleName() << " \t" 
               << std::setprecision(5) << g_Energy*MeV<< " \t"
               << std::setprecision(5) << 0. << " \t"
               << std::setprecision(5) << 0. << " \t"
               << std::setprecision(5) << position*mm << " \t"
               << std::setprecision(5) << 0. << " \t"
               << std::setprecision(5) << 0. << " \t"
               << std::setprecision(5) << 1. << " \t"
               << G4endl;
 }
}

else {
  G4cout<<" WHAT? "<<G4endl;
  G4ParticleDefinition* particle
//           = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
             = G4IonTable::GetIonTable()->GetIon(7, 15, 0.);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,-1.));
  fParticleGun->SetParticleEnergy(40.23*MeV);
  G4double position = -1.00000215; // chamber r + havar/2
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*m,0.*m,position*m));

/*    G4ThreeVector oldPosition = G4ThreeVector(0.,0.,0.);
    G4double x0 = oldPosition.x();
    G4double y0 = oldPosition.y() + (2*G4UniformRand()-1.)*0.5;
    G4double z0 = oldPosition.z() + (2*G4UniformRand()-1.)*0.5;
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
*/
    fParticleGun->GeneratePrimaryVertex(anEvent);
//    fParticleGun->SetParticlePosition(oldPosition);     


//    energiaIniziale = fParticleGun -> GetParticleEnergy();


}

outf.close();

}

