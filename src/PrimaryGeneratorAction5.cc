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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction3.cc
/// \brief Implementation of the PrimaryGeneratorAction3 class
//
//
// $Id: PrimaryGeneratorAction3.cc 83872 2014-09-20 22:23:50Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction5.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "SteppingAction.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction5::PrimaryGeneratorAction5(G4ParticleGun* gun)
: fParticleGun(gun)
{
  //solid angle
  //

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction5::~PrimaryGeneratorAction5()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction5::GeneratePrimaries(G4Event* anEvent, G4double g_Energy, G4String outputName, G4double rbeam)
{

    std::ofstream outf;
    outf.open(outputName+".txt", std::ios::out | std::ios::app );
    ////  ROOT GEN ////

       G4double      root_u=931.49432/1000.; //
    // unita' di massa (rif. 12C)
    //11B 11.0093054
    //2H 2.0141017778

       G4int root_out_particles = 3;     // n. output particles
       G4double     root_mp=25.98689186; // massa proiettile             //
       G4double     root_mt=2.01363;     // massa target                 //
       G4double     root_m1=1.007276466879;    // massa prima part. uscente    //
       G4double     root_m2=25.98259297;    // massa seconda part. uscente  //
//       G4double     root_m1=4.002603;    // massa prima part. uscente    //
//       G4double     root_m2=22.9897692820;    // massa seconda part. uscente  //

       G4double     root_m3=1.007276466879;    // massa terza part. uscente SPECT//
    //(Momentum, Energy units are Gev/C, GeV)

    G4double ExcEne1 = 0.;
    G4double ExcEne2 = 0.;
    G4double ExcEne3 = 0.;// 6.049; // <-- alfa01: 6.049 MeV alfa0g: 6.130 MeV  from MeV to GeV
   //(Momentum, Energy units are Gev/C, GeV)

       G4double root_masses[3] = {root_u*root_m1+ExcEne1,root_u*root_m2+ExcEne2,root_u*root_m3+ExcEne3};
       TGenPhaseSpace TheRootEvent;

       TRandom3* ran = new TRandom3(0);
       ran->SetSeed(0);

       Double_t rnum = ran->Rndm();
       Int_t rnnum=rnum*10000000;
       std::default_random_engine generator;
       generator.seed(rnnum);
       std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015);

       double randomNumber = distribution(generator);
    // energia proiettile
       Double_t     root_ep=randomNumber/1000.;

    // impulso proiettile
       G4double     root_pp=-1.;
       if(root_mp==0){
           root_pp=root_ep;
       }
       else {
           root_pp=sqrt(2.*root_u*root_mp*root_ep); // root_ep/299792458.;
       }


    //   random interaction point here ?
       TLorentzVector root_target(0.0, 0.0, 0.0, root_u*root_mt);
       TLorentzVector root_beam(0.0, 0.0, root_pp, root_ep+root_u*root_mp);
       TLorentzVector rootW = root_beam + root_target;

    //  I NEED TO SET THE RANDOM SEED FOR THE GENERATION SOMEWHERE HERE
   bool GO_ON=true;
   G4int countREIT = 0;
//   double prelmax = sqrt(2.*4.73*(4.*16.)/20.);
   G4double root_weight = -1000;
   TLorentzVector * root_p1, * root_p2, * root_p3;
   TVector3 root_pp1,root_pp2,root_pp3,root_dir1,root_dir2,root_dir3;
   while(GO_ON){
       if(!GO_ON) G4cout<<"Aspetta aspetta aspetta"<<G4endl;
       G4int seed2 = (int) G4UniformRand()*10000;
       gRandom->SetSeed(seed2);

       G4bool isThatOK = TheRootEvent.SetDecay(rootW, root_out_particles, root_masses);

        root_weight = TheRootEvent.Generate();

       root_p1 = TheRootEvent.GetDecay(0);
       root_p2 = TheRootEvent.GetDecay(1);
       root_p3 = TheRootEvent.GetDecay(2);

       root_pp1=root_p1->Vect();
       root_pp2=root_p2->Vect();
       root_pp3=root_p3->Vect();

       root_dir1=root_pp1.Unit();
       root_dir2=root_pp2.Unit();
       root_dir3=root_pp3.Unit();
//	G4double dotproduct = root_dir1.dot(root_dir2);//get the dot product between the beam and the spectator alpha to ensure QF kinematics
//	TODO: sample the relative momentum of the alpha and oxygen
//	G4ThreeVector prel = pprod1-pprod2;//relative momentum of alpha and oxygen
//	G4double prelt=prel.mag();//magnitude of prel
//	G4cout<<prelt<<"\t"<<prelmax<<G4endl;
        if((1.e3)*(root_p3->E() - root_u*root_m3) < 0.85){ // QF CONDITION : MARCO
           	GO_ON=false; break;
        }
//        else G4cout<<"Spectator Energy = "<< (1.e3)*(root_p3->E() - root_u*root_m3) << G4endl;
/*	if(dotproduct>0.5 && prelt<prelmax) {//ensure spectator is in QF kinematics by having the incoming beam and spectator have an angle of <60
		GO_ON=false;
	}
       if( abs(root_p1->Angle(root_p2->Vect())-root_p3->Angle(root_p2->Vect()))<0.001 ){
	  GO_ON=FALSE;
       }*/
       countREIT++;
   }

         out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("proton");
//         out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
         G4int TheA = 26, TheZ = 12;
//         G4int TheA = 23, TheZ = 11;
         out_particle2 = G4IonTable::GetIonTable()->GetIon(TheZ, TheA, 0.);
//         G4ParticleDefinition* out_particle3 = G4IonTable::GetIonTable()->GetIon(TheZ, TheA, 0.);
         G4ParticleDefinition * out_particle3 = G4ParticleTable::GetParticleTable()->FindParticle("proton");

         G4double E1 = (1.e3)*(root_p1->E() - root_u*root_m1);
         G4double E2 = (1.e3)*(root_p2->E() - root_u*root_m2);
         G4double E3 = (1.e3)*(root_p3->E() - root_u*root_m3);
   
         if (root_weight>0.){
           G4double t0 = (1.E-9)*s;

           G4ThreeVector oldPosition = G4ThreeVector(0.,0., (0.00057537+0.000266)/2.); //FIX THIS
           //rbeam = 0.5*mm; /* sin(targetrot*deg)* */

           G4double r = sqrt(G4UniformRand());
           G4double theta = 2.*pi*G4UniformRand()*rad;
           G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
           G4double y0 = oldPosition.y() + r*rbeam*sin(theta);

           G4double z0 = oldPosition.z() ; // + (2*G4UniformRand()-1.)*0.00057537 *mm;
       /* target rotation phi */
           G4double phi = targetphi;

           z0 = - x0*sin(phi) + oldPosition.z()*cos(phi);
           x0 = x0*cos(phi) + oldPosition.z()*sin(phi);

           fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

           fParticleGun->SetParticleDefinition(out_particle1);
           fParticleGun->SetParticleEnergy(E1*MeV);
           fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));

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
//	   G4cout<<"After "<<countREIT<<" tries."<<G4endl;
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

       }
       else
       {
           G4cout<<"Just a warning for root_weight <= 0. YOU MESSED UP! "<<G4endl;
       } // add "abort event" ?
      
       outf.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
