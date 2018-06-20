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
  det = DC; // get the pointer to the current DC
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
  g_Energy = 5.; // MeV
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

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

   std::ofstream outf;
   outf.open("GROOT.txt", std::ios::out | std::ios::app );
//   targetphi = 35.*deg;

if(fuffa==2){
//    G4cout<<"LET'S BEGIN!"<<G4endl;
//    gROOT.Reset();

////  ROOT GEN ////
// TApplication theApp("App", NULL, NULL);
//G4double root_E = 0.;
//G4ThreeVector root_P = G4ThreeVector(0.,0.,0.);
   G4double root_u = 931.49432/1000.; //
    // unita' di massa (rif. 12C)
   G4double phi = targetphi;
//11B 11.0093054
//2H 2.0141017778

////////// reazione 24Mg(g,a)20Ne  //////////////////////////////////
    G4int root_out_particles = 2;     // n. output particles
//    G4double     root_mp=0;           // massa proiettile             //
//    G4double     root_mt=23.985041;   // massa target                 //
//    G4double     root_m1=4.002603;    // massa prima part. uscente    //
//    G4double     root_m2=19.992440;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 27Al(g,a)23Na  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=26.981538;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=22.989769;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 27Al(g,p)26Mg  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=26.981538;   // massa target                 //
//    G4double     m1=1.007825;    // massa prima part. uscente    //
//    G4double     m2=25.982592;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 56Fe(g,a)52Cr  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=55.934937;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=51.940507;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 56Fe(g,p)55Mn  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=55.934937;   // massa target                 //
//    G4double     m1=1.007825;    // massa prima part. uscente    //
//    G4double     m2=54.938045;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 13C(g,a)9Be  ////////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=13.003355;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=9.012182;    // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 16O(g,a)12C  ////////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=15.994915;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=12.000000;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 28Si(g,a)24Mg  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=27.976926;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=23.985042;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 25Mg(g,a)21Ne  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=24.985837;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=20.993847;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
    
////////// reazione 26Mg(g,a)22Ne  //////////////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=25.982593;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=21.991385;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 96Ru(g,a)92Mo @ 9.3 MeV  ////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=95.907598;   // massa target                 //
//    G4double     m1=4.002603;    // massa prima part. uscente    //
//    G4double     m2=91.906811;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 74Se(g,p)73As @ 11.1 MeV ////////////////////////
//    G4double     mp=0;           // massa proiettile             //
//    G4double     mt=73.922476;   // massa target                 //
//    G4double     m1=1.007825;    // massa prima part. uscente    //
//    G4double     m2=72.923825;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 7Li(g,t)4He @ > 2.5 MeV ////////////////////////
//    G4double     root_mp=0;           // massa proiettile             //
//    G4double     root_mt=7.016004548;   // massa target                 //
//    G4double     root_m1=3.0160492;   // massa prima part. uscente    //
//    G4double     root_m2=4.002603;    // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////

////////// reazione 6Li(g,d)4He @ 5 MeV ////////////////////////
//    G4double     root_mp=0;           // massa proiettile             //
//    G4double     root_mt=6.015122795;   // massa target                 //
//    G4double     root_m1=4.002603;   // massa prima part. uscente    //
//    G4double     root_m2=2.01363;    // massa seconda part. uscente  //
////////////////////////////////////////////////////////////////////

////////// reazione 19F(g,a)15N @ > X MeV ////////////////////////
//    G4double     root_mp=0;           // massa proiettile             //
//    G4double     root_mt=18.998403;   // massa target                 //
//    G4double     root_m1=15.0001089;   // massa prima part. uscente    //
//    G4double     root_m2=4.002603;    // massa seconda part. uscente  //
////////////////////////////////////////////////////////////////////

////////// reazione 19F(g,p)18O @ > X MeV ////////////////////////
//    G4double     root_mp=0;           // massa proiettile             //
//    G4double     root_mt=;   // massa target                 //
//    G4double     root_m1=;   // massa prima part. uscente    //
//   G4double     root_m2=;    // massa seconda part. uscente  //
////////////////////////////////////////////////////////////////////

////////// reazione 112Sn(g,a)108Cd > 1.83034 MeV  ////////////////////////
    G4double     root_mp=0;           // massa proiettile             //
    G4double     root_mt=111.904826;   // massa target                 //
    G4double     root_m1=4.002603;    // massa prima part. uscente    //
    G4double     root_m2=107.90418;   // massa seconda part. uscente  //
/////////////////////////////////////////////////////////////////////
/*    G4double root_deg=180./3.141592654;
    G4double root_radius=50.;
    G4int root_geve=0;
*/  
    //(Momentum, Energy units are Gev/C, GeV)

    G4double root_masses[2] = {root_u*root_m1,root_u*root_m2};
    
    TGenPhaseSpace TheRootEvent;
    
 //   G4int seed = 0;
    TRandom3* ran = new TRandom3(0);
    ran->SetSeed(0);

        Double_t rnum = ran->Rndm();
        
        Int_t rnnum=rnum*10000000;
        
//        cout << "rnum= " << rnnum << "\n";
        
        std::default_random_engine generator;
        generator.seed(rnnum);
//        std::normal_distribution<double> distribution(/*mean=*/11.7, /*stddev=*/0.015);
        std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015);

//  Double_t Gaus(mean,sigma);
        double randomNumber = distribution(generator);
        
        Double_t root_ep=randomNumber/1000.;

//        Double_t     root_ep=11.7/1000.;

    // energia proiettile
    
//    cout << "Ebeam= " << ep << "\n";
    
    G4double     root_pp=root_ep;
    // impulso proiettile

///////   random interaction point here ?    
    TLorentzVector root_target(0.0, 0.0, 0.0, root_u*root_mt);
    TLorentzVector root_beam(0.0, 0.0, root_pp, root_ep);
    TLorentzVector rootW = root_beam + root_target;

//  I NEED TO SET THE RANDOM SEED FOR THE GENERATION SOMEWHERE HERE

    G4int seed2 = (int) G4UniformRand()*10000;
    gRandom->SetSeed(seed2);
  //  outf << "gRandom set to : " << seed2 << G4endl; 
//    gRandom->WriteRandom("random.dat");
  
    G4bool isThatOK = TheRootEvent.SetDecay(rootW, root_out_particles, root_masses);
//    G4cout<<"It's ok. "<<isThatOK<<G4endl;
    G4double root_weight = TheRootEvent.Generate();

    TLorentzVector *root_p1 = TheRootEvent.GetDecay(0);
    TLorentzVector *root_p2 = TheRootEvent.GetDecay(1);

    TVector3 root_pp1=root_p1->Vect();
    TVector3 root_pp2=root_p2->Vect();
    
    TVector3 root_dir1=root_pp1.Unit();
    TVector3 root_dir2=root_pp2.Unit();

    G4ParticleDefinition* out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
//    G4ParticleDefinition* out_particle1 = G4IonTable::GetIonTable()->GetIon(7, 15, 0.);
//    G4ParticleDefinition* out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
//    G4cout<<"Out part1: "<<out_particle1->GetParticleName()<<G4endl;

/*    static G4ParticleDefinition* part = 0;
    if(!part){
        part = G4IonTable::GetIonTable()->GetIon(10,20,0.);
        fParticleGun->SetParticleDefinition(part);
    }
    G4cout<<"Out part2 temp: "<<part->GetParticleName()<<G4endl;*/
    G4int TheA = 60, TheZ = 48;
//    G4ParticleDefinition* out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
//    G4ParticleDefinition* out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
//    G4ParticleDefinition* out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    G4ParticleDefinition* out_particle2 = G4IonTable::GetIonTable()->GetIon(TheZ, TheA, 0.);
//    G4cout<<"Out part2: "<<out_particle2->GetParticleName()<<G4endl;

//    gRandom->WriteRandom("random2.dat");

  //this function is called at the begining of event
  //
/*G4cout<<"==================================="<<G4endl;
   G4cout<<root_weight<<"\t"<<root_p1->E()<<"\t"<<root_p2->E()<<G4endl;
G4cout<<"==================================="<<G4endl;
*/
//  if (root_weight>0){

  G4double E1 = (1.e3)*(root_p1->E() - root_u*root_m1);
  G4double E2 = (1.e3)*(root_p2->E() - root_u*root_m2);

/*  G4cout<<root_p1->E() <<"\t"<<root_p2->E() <<G4endl;
  G4cout<<E1<<"  ------------------------------------  "<<E2<<G4endl;
  G4cout<<root_m1<<"\t"<<root_m2 <<G4endl;
  G4cout<<root_weight<<"\t"<<E1<<"\t"<<E2<<G4endl;
*/
  if (root_weight>0.){
    G4double t0 = 0.;
//    G4double X_target = DetectorConstruction::targetLayerXPosition;
//    G4cout<<"WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO: "<<X_target<<G4endl;
    G4ThreeVector oldPosition = G4ThreeVector(0.,0.,0.);
    G4double rbeam = 0.75*cm; // 0.5*mm; /* sin(targetrot*deg)* */ 
// (b*R*cos(2*pi*a/b), b*R*sin(2*pi*a/b)) CHECK!
//    G4double x0 = oldPosition.x() + b1*rbeam*cos(twopi*a1/b1);
//    G4double y0 = oldPosition.y() + b1*rbeam*sin(twopi*a1/b1);

//    G4double cosAlpha = 2.*G4UniformRand() - 1.;
//    G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
//    G4double x0 = oldPosition.x() + r*rbeam*cosAlpha;
//    G4double y0 = oldPosition.y() + r*rbeam*sinAlpha;
 
    G4double r = sqrt(G4UniformRand());
    G4double theta = 2.*pi*G4UniformRand()*rad;
    G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
    G4double y0 = oldPosition.y() + r*rbeam*sin(theta);

//    G4double y0 = oldPosition.y() + ( (a1>=0.) ? b1*rbeam*sin(acos(a1)) : -b1*rbeam*sin(acos(a1)));

//   std::ofstream yiyo;
//   yiyo.open("yiyo.txt", std::ios::out | std::ios::app );    
//   yiyo<<det->GetTargetThickness()<<"\t"<<det->GetTargetPosZ()<<G4endl;
//   yiyo.close();
//   G4cout<<det->GetTargetThickness()<<"\t"<<det->GetTargetPosZ()<<G4endl;
    G4double z0 = (det->GetTargetPosZ() +(2*G4UniformRand()-1.)*det->GetTargetThickness()); //<-- already in best units I hope
//   G4double z0 = oldPosition.z(); // + (2*G4UniformRand()-1.)*0.00057537 *mm;

/* target rotation phi */

      G4double phi = targetphi; // 0 - 35

//      z0 = - y0*sin(phi) + oldPosition.z()*cos(phi);
//      y0 = y0*cos(phi) + oldPosition.z()*sin(phi);
//    y0 = y0*(sin(phi) + cos(phi));
//      G4double x1 = x0*cos(phi) + z0*sin(phi);
//      G4double z1 = -x0*sin(phi) + z0*cos(phi);

    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

    fParticleGun->SetParticleDefinition(out_particle1);
    fParticleGun->SetParticleEnergy(E1*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));
//    G4cout<<root_p1->E()<<"\t"<<root_dir1.x()<<"\t"<<root_dir1.y()<<"\t"<<root_dir1.z()<<"\t"<<out_particle1->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << out_particle1->GetParticleName() << " \t" 
               << std::setprecision(5) << E1*MeV<< " \t"
               << std::setprecision(5) << x0/cm << " \t"
               << std::setprecision(5) << y0/cm << " \t"
               << std::setprecision(5) << z0/cm << " \t"
               << std::setprecision(5) << root_dir1(0) << " \t"
               << std::setprecision(5) << root_dir1(1) << " \t"
               << std::setprecision(5) << root_dir1(2) << " \t"
               << G4endl;

//    t0 = (1.E-3)*s;
    fParticleGun->SetParticleDefinition(out_particle2);
    fParticleGun->SetParticleEnergy(E2*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir2(0),root_dir2(1),root_dir2(2)));
//    G4cout<<root_p2->E()<<"\t"<<root_p2->X()<<"\t"<<root_p2->Y()<<"\t"<<root_p2->Z()<<"\t"<<out_particle2->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << out_particle2->GetParticleName() << " \t" 
               << std::setprecision(5) << E2*MeV<< " \t"
               << std::setprecision(5) << x0/cm << " \t"
               << std::setprecision(5) << y0/cm << " \t"
               << std::setprecision(5) << z0/cm << " \t"
               << std::setprecision(5) << root_dir2(0) << " \t"
               << std::setprecision(5) << root_dir2(1) << " \t"
               << std::setprecision(5) << root_dir2(2) << " \t"
               << G4endl;

//    fParticleGun->SetParticlePosition(oldPosition);     
//    energiaIniziale = fParticleGun -> GetParticleEnergy();
/*    root_p1->Delete();
    root_p2->Delete();
    TheRootEvent.Delete();
    root_target.Delete();
    root_beam.Delete();
    rootW.Delete();
*/
}
   else {G4cout<<"Just a warning for root_weight <= 0."<<G4endl;} // add "abort event" ?

}
else if(fuffa==3){
//    G4cout<<"LET'S BEGIN!"<<G4endl;
//    gROOT.Reset();

////  ROOT GEN ////
// TApplication theApp("App", NULL, NULL);
//G4double root_E = 0.;
//G4ThreeVector root_P = G4ThreeVector(0.,0.,0.);
   G4double      root_u=931.49432/1000.; //
    // unita' di massa (rif. 12C)

//11B 11.0093054
//2H 2.0141017778

    G4int root_out_particles = 3;     // n. output particles
    G4double     root_mp=0;           // massa proiettile             //
    G4double     root_mt=12.;         // massa target                 //
    G4double     root_m1=4.002603;    // massa prima part. uscente    //
    G4double     root_m2=4.002603;    // massa seconda part. uscente  //
    G4double     root_m3=4.002603;    // massa terza part. uscente  //
    //(Momentum, Energy units are Gev/C, GeV)
    G4double root_masses[3] = {root_u*root_m1,root_u*root_m2,root_u*root_m3};
    
    TGenPhaseSpace TheRootEvent;
    
 //   G4int seed = 0;
    TRandom3* ran = new TRandom3(0);
    ran->SetSeed(0);

        Double_t rnum = ran->Rndm();
        
        Int_t rnnum=rnum*10000000;
        
//        cout << "rnum= " << rnnum << "\n";
        
        std::default_random_engine generator;
        generator.seed(rnnum);
        std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015);

//  Double_t Gaus(mean,sigma);
        double randomNumber = distribution(generator);
        
        Double_t     root_ep=randomNumber/1000.;

//        Double_t     root_ep=11.7/1000.;

    // energia proiettile
    
//    cout << "Ebeam= " << ep << "\n";
    
    G4double     root_pp=root_ep;
    // impulso proiettile

///////   random interaction point here ?    
    TLorentzVector root_target(0.0, 0.0, 0.0, root_u*root_mt);
    TLorentzVector root_beam(0.0, 0.0, root_pp, root_ep);
    TLorentzVector rootW = root_beam + root_target;

//  I NEED TO SET THE RANDOM SEED FOR THE GENERATION SOMEWHERE HERE

    G4int seed2 = (int) G4UniformRand()*10000;
    gRandom->SetSeed(seed2);
  //  outf << "gRandom set to : " << seed2 << G4endl; 
//    gRandom->WriteRandom("random.dat");
  
    G4bool isThatOK = TheRootEvent.SetDecay(rootW, root_out_particles, root_masses);
//    G4cout<<"It's ok. "<<isThatOK<<G4endl;
    G4double root_weight = TheRootEvent.Generate();

    TLorentzVector *root_p1 = TheRootEvent.GetDecay(0);
    TLorentzVector *root_p2 = TheRootEvent.GetDecay(1);
    TLorentzVector *root_p3 = TheRootEvent.GetDecay(2);

    TVector3 root_pp1=root_p1->Vect();
    TVector3 root_pp2=root_p2->Vect();
    TVector3 root_pp3=root_p3->Vect();
    
    TVector3 root_dir1=root_pp1.Unit();
    TVector3 root_dir2=root_pp2.Unit();
    TVector3 root_dir3=root_pp3.Unit();

    G4ParticleDefinition* out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
//    G4cout<<"Out part1: "<<out_particle1->GetParticleName()<<G4endl;
    G4int TheA = 4, TheZ = 2;
//    G4ParticleDefinition* out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("20Ne");
    G4ParticleDefinition* out_particle2 = G4IonTable::GetIonTable()->GetIon(TheZ, TheA, 0.);
//    G4cout<<"Out part2: "<<out_particle2->GetParticleName()<<G4endl;
    G4ParticleDefinition* out_particle3 = G4IonTable::GetIonTable()->GetIon(TheZ, TheA, 0.);
//    gRandom->WriteRandom("random2.dat");

  //this function is called at the begining of event
  //
/*G4cout<<"==================================="<<G4endl;
   G4cout<<root_weight<<"\t"<<root_p1->E()<<"\t"<<root_p2->E()<<G4endl;
G4cout<<"==================================="<<G4endl;
*/
//  if (root_weight>0){

  G4double E1 = (1.e3)*(root_p1->E() - root_u*root_m1);
  G4double E2 = (1.e3)*(root_p2->E() - root_u*root_m2);
  G4double E3 = (1.e3)*(root_p3->E() - root_u*root_m3);

/*  G4cout<<root_p1->E() <<"\t"<<root_p2->E() <<G4endl;
  G4cout<<E1<<"  ------------------------------------  "<<E2<<G4endl;
  G4cout<<root_m1<<"\t"<<root_m2 <<G4endl;
  G4cout<<root_weight<<"\t"<<E1<<"\t"<<E2<<G4endl;
*/
  if (root_weight>0.){
    G4double t0 = (1.E-9)*s;

    G4ThreeVector oldPosition = G4ThreeVector(0.,0., (0.00057537+0.000266)/2.); //FIX THIS
    G4double rbeam = 0.5*mm; /* sin(targetrot*deg)* */ 
// (b*R*cos(2*pi*a/b), b*R*sin(2*pi*a/b)) CHECK!
//    G4double x0 = oldPosition.x() + b1*rbeam*cos(twopi*a1/b1);
//    G4double y0 = oldPosition.y() + b1*rbeam*sin(twopi*a1/b1);

//    G4double cosAlpha = 2.*G4UniformRand() - 1.;
//    G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
//    G4double x0 = oldPosition.x() + r*rbeam*cosAlpha;
//    G4double y0 = oldPosition.y() + r*rbeam*sinAlpha;
 
    G4double r = sqrt(G4UniformRand());
    G4double theta = 2.*pi*G4UniformRand()*rad;
    G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
    G4double y0 = oldPosition.y() + r*rbeam*sin(theta);

//    G4double y0 = oldPosition.y() + ( (a1>=0.) ? b1*rbeam*sin(acos(a1)) : -b1*rbeam*sin(acos(a1)));

    G4double z0 = oldPosition.z() ; // + (2*G4UniformRand()-1.)*0.00057537 *mm;
/* target rotation phi */

//    G4double phi = 35.*deg;
      G4double phi = targetphi; 
/*    x0 = x0*(cos(phi) - sin(phi));
//    y0 = y0*(sin(phi) + cos(phi));
      z0 = z0+ x0*(sin(phi) + cos(phi));
*/
      z0 = - x0*sin(phi) + oldPosition.z()*cos(phi);
      x0 = x0*cos(phi) + oldPosition.z()*sin(phi);
//    y0 = y0*(sin(phi) + cos(phi));


    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

    fParticleGun->SetParticleDefinition(out_particle1);
    fParticleGun->SetParticleEnergy(E1*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));
//    G4cout<<root_p1->E()<<"\t"<<root_dir1.x()<<"\t"<<root_dir1.y()<<"\t"<<root_dir1.z()<<"\t"<<out_particle1->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << out_particle1->GetParticleName() << " \t" 
               << std::setprecision(5) << E1*MeV<< " \t"
               << std::setprecision(5) << x0 << " \t"
               << std::setprecision(5) << y0 << " \t"
               << std::setprecision(5) << z0 << " \t"
               << std::setprecision(5) << root_dir1(0) << " \t"
               << std::setprecision(5) << root_dir1(1) << " \t"
               << std::setprecision(5) << root_dir1(2) << " \t"
               << G4endl;

//    t0 = (1.E-3)*s;
    fParticleGun->SetParticleDefinition(out_particle2);
    fParticleGun->SetParticleEnergy(E2*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir2(0),root_dir2(1),root_dir2(2)));
//    G4cout<<root_p2->E()<<"\t"<<root_p2->X()<<"\t"<<root_p2->Y()<<"\t"<<root_p2->Z()<<"\t"<<out_particle2->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << out_particle2->GetParticleName() << " \t" 
               << std::setprecision(5) << E2*MeV<< " \t"
               << std::setprecision(5) << x0 << " \t"
               << std::setprecision(5) << y0 << " \t"
               << std::setprecision(5) << z0 << " \t"
               << std::setprecision(5) << root_dir2(0) << " \t"
               << std::setprecision(5) << root_dir2(1) << " \t"
               << std::setprecision(5) << root_dir2(2) << " \t"
               << G4endl;

    fParticleGun->SetParticleDefinition(out_particle3);
    fParticleGun->SetParticleEnergy(E3*MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir3(0),root_dir3(1),root_dir3(2)));
//    G4cout<<root_p2->E()<<"\t"<<root_p2->X()<<"\t"<<root_p2->Y()<<"\t"<<root_p2->Z()<<"\t"<<out_particle2->GetParticleName()<<G4endl;
    fParticleGun->SetParticleTime(t0);
    fParticleGun->GeneratePrimaryVertex(anEvent);

         outf  << "Generated Primary: "
               << out_particle3->GetParticleName() << " \t" 
               << std::setprecision(5) << E3*MeV<< " \t"
               << std::setprecision(5) << x0 << " \t"
               << std::setprecision(5) << y0 << " \t"
               << std::setprecision(5) << z0 << " \t"
               << std::setprecision(5) << root_dir3(0) << " \t"
               << std::setprecision(5) << root_dir3(1) << " \t"
               << std::setprecision(5) << root_dir3(2) << " \t"
               << G4endl;

//    fParticleGun->SetParticlePosition(oldPosition);     
//    energiaIniziale = fParticleGun -> GetParticleEnergy();
/*    root_p1->Delete();
    root_p2->Delete();
    TheRootEvent.Delete();
    root_target.Delete();
    root_beam.Delete();
    rootW.Delete();
*/
}
   else {G4cout<<"Just a warning for root_weight <= 0."<<G4endl;} // add "abort event" ?

}

else {

  //    G4double X_target = DetectorConstruction::targetLayerXPosition;
//    G4cout<<"WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO: "<<X_target<<G4endl;
    G4ThreeVector oldPosition = G4ThreeVector(0.,0.,0.);
    G4double rbeam = 0.75*cm; // 0.5*mm; /* sin(targetrot*deg)* */ 
// (b*R*cos(2*pi*a/b), b*R*sin(2*pi*a/b)) CHECK!
//    G4double x0 = oldPosition.x() + b1*rbeam*cos(twopi*a1/b1);
//    G4double y0 = oldPosition.y() + b1*rbeam*sin(twopi*a1/b1);

//    G4double cosAlpha = 2.*G4UniformRand() - 1.;
//    G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
//    G4double x0 = oldPosition.x() + r*rbeam*cosAlpha;
//    G4double y0 = oldPosition.y() + r*rbeam*sinAlpha;
 
    G4double r = sqrt(G4UniformRand());
    G4double theta = 2.*pi*G4UniformRand()*rad;
    G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
    G4double y0 = oldPosition.y() + r*rbeam*sin(theta);


    G4double z0 = oldPosition.z()-200.*cm; // + (2*G4UniformRand()-1.)*0.00057537 *mm;

    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));


  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(g_Energy*MeV);
//  G4double position = 1.;
//  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position*cm));

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

