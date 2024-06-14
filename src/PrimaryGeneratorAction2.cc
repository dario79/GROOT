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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction2.cc
/// \brief Implementation of the PrimaryGeneratorAction2 class
//
//
// $Id: PrimaryGeneratorAction2.cc 83872 2014-09-20 22:23:50Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction2.hh"
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

PrimaryGeneratorAction2::PrimaryGeneratorAction2(G4ParticleGun* gun)
: fParticleGun(gun)
{
  //solid angle
  //

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction2::~PrimaryGeneratorAction2()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction2::GeneratePrimaries(G4Event* anEvent, G4double state, G4double g_Energy, G4String outputType, G4String outputName,
                                                G4int Part1A, G4int Part1Z, G4int Part2A, G4int Part2Z, G4double Part1Mass, G4double Part2Mass,
                                                G4double ProjectileMass, G4double TargetMass, G4double rbeam, DetectorConstruction* det)
{

    std::ofstream outf;
    outf.open(outputName+".txt", std::ios::out | std::ios::app );

   ////  ROOT GEN ////
   // TApplication theApp("App", NULL, NULL);

       G4double root_u = 931.49432/1000.; //
       G4int root_out_particles = 2;   // n. output particles
   // unita' di massa (rif. 12C)
   //G4double phi = targetphi;
   //11B 11.0093054
   //2H 2.0141017778

    ////////// reazione 24Mg(g,a)20Ne  //////////////////////////////////
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
    //    G4double     root_mp=0;           // massa proiettile             //
    //    G4double     root_mt=111.904826;   // massa target                 //
    //    G4double     root_m1=4.002603;    // massa prima part. uscente    //
    //    G4double     root_m2=107.90418;   // massa seconda part. uscente  //
    /////////////////////////////////////////////////////////////////////

    ////////// reazione 6Li(g,d)4He @ 5 MeV ////////////////////////
    //    G4double     root_mp=0;           // massa proiettile             //
    //    G4double     root_mt=6.015122795;   // massa target                 //
    //   G4double     root_m1=4.002603;   // massa prima part. uscente    //
    //    G4double     root_m2=2.01363;    // massa seconda part. uscente  //
    ////////////////////////////////////////////////////////////////////

    //(Momentum, Energy units are Gev/C, GeV)

    G4double ExcEne1 = 0.;
    G4double ExcEne2 = state/1000.;// 6.049; // <-- alfa01: 6.049 MeV alfa0g: 6.130 MeV  from MeV to GeV
   //(Momentum, Energy units are Gev/C, GeV)

   G4double root_masses[2] = {root_u*(Part1Mass)+ExcEne1,root_u*(Part2Mass)+ExcEne2};
//    G4cout<<Part1Mass<<"\t"<<root_masses[1]<<"\t"<<state<<G4endl;
//G4double root_masses[2] = {root_u*Part1Mass,root_u*Part2Mass};
       TGenPhaseSpace TheRootEvent;

       TRandom3* ran = new TRandom3(0);
       ran->SetSeed(0);
       Double_t rnum = ran->Rndm();
       Int_t rnnum=rnum*10000000;

       std::default_random_engine generator;
       generator.seed(rnnum);
       std::normal_distribution<double> distribution(/*mean=*/g_Energy, /*stddev=*/0.015);
       double randomNumber = distribution(generator);

       Double_t root_ep=randomNumber/1000.;
       G4double     root_pp=-1.;
//      impulso proiettile
       if(ProjectileMass==0){
           root_pp=root_ep;
       }
       else {
           root_pp=sqrt(2.*root_u*ProjectileMass*root_ep); // root_ep/299792458.;
       }

   ///////   random interaction point here ?

        TLorentzVector root_target(0.0, 0.0, 0.0, root_u*TargetMass);
        TLorentzVector root_beam(0.0, 0.0, root_pp, root_ep +ProjectileMass*root_u);
        TLorentzVector rootW = root_beam + root_target;

   //  I NEED TO SET THE RANDOM SEED FOR THE GENERATION SOMEWHERE HERE

        G4int seed2 = (int) G4UniformRand()*10000;
        gRandom->SetSeed(seed2);

        G4bool isThatOK = TheRootEvent.SetDecay(rootW, root_out_particles, root_masses);
        G4double root_weight = TheRootEvent.Generate();

        TLorentzVector *root_p1 = TheRootEvent.GetDecay(0);
        TLorentzVector *root_p2 = TheRootEvent.GetDecay(1);

        TVector3 root_pp1=root_p1->Vect();
        TVector3 root_pp2=root_p2->Vect();

        TVector3 root_dir1=root_pp1.Unit();
        TVector3 root_dir2=root_pp2.Unit();

       // out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle(Part1_Name);

        if(Part1A==0) {
             out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        }
        else  if(Part1A==1 && Part1Z==0) {
             out_particle1 = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        }
        else {
             out_particle1 = G4IonTable::GetIonTable()->GetIon(Part1Z, Part1A, 0.);
        }
        if(Part2A==0) {
             out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        }
        else  if(Part2A==1 && Part2Z==0) {
             out_particle2 = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        }
        else {
             out_particle2 = G4IonTable::GetIonTable()->GetIon(Part2Z, Part2A, 0.);
        }

        G4double E1 = (1.e3)*(root_p1->E() - root_u*Part1Mass);
        G4double E2 = (1.e3)*(root_p2->E() - root_u*Part2Mass);


        if (root_weight>0.)
        {
           G4double t0 = 0.;
           G4ThreeVector oldPosition = G4ThreeVector(0.,0.,0.);
           G4double r = sqrt(G4UniformRand());
           G4double theta = 2.*pi*G4UniformRand()*rad;
           G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
           G4double y0 = oldPosition.y() + r*rbeam*sin(theta);
           G4double z0 = (det->GetTargetPosZ()) +(2*G4UniformRand()-1.)*(det->GetTargetThickness())/2.; //<-- already in best units I hope
       /* target rotation phi */
       //G4double phi = targetphi; // 0 - 35

           fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

           fParticleGun->SetParticleDefinition(out_particle1);
           fParticleGun->SetParticleEnergy(E1*MeV);
           fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir1(0),root_dir1(1),root_dir1(2)));

           fParticleGun->SetParticleTime(t0);
           fParticleGun->GeneratePrimaryVertex(anEvent);

           if(outputType != "root")
           {

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
           }

           if(outputType != "txt")
           {
               G4String name = out_particle1->GetParticleName();
               auto analysisManager = G4AnalysisManager::Instance();
               analysisManager->FillNtupleSColumn(0, 0, name);
               analysisManager->FillNtupleDColumn(0, 1, E1);
               analysisManager->FillNtupleDColumn(0, 2, x0);
               analysisManager->FillNtupleDColumn(0, 3, y0);
               analysisManager->FillNtupleDColumn(0, 4, z0);
               analysisManager->FillNtupleDColumn(0, 5, root_dir1(0));
               analysisManager->FillNtupleDColumn(0, 6, root_dir1(1));
               analysisManager->FillNtupleDColumn(0, 7, root_dir1(2));
               analysisManager->AddNtupleRow(0);
           }

       //    t0 = (1.E-3)*s;
           fParticleGun->SetParticleDefinition(out_particle2);
           fParticleGun->SetParticleEnergy(E2*MeV);
           fParticleGun->SetParticleMomentumDirection(G4ThreeVector(root_dir2(0),root_dir2(1),root_dir2(2)));

           fParticleGun->SetParticleTime(t0);
           fParticleGun->GeneratePrimaryVertex(anEvent);

           if(outputType != "root")
           {
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
           }

           if(outputType != "txt")
           {

               G4String name = out_particle2->GetParticleName();
               auto analysisManager = G4AnalysisManager::Instance();
               analysisManager->FillNtupleSColumn(0, 0, name);
               analysisManager->FillNtupleDColumn(0, 1, E2);
               analysisManager->FillNtupleDColumn(0, 2, x0);
               analysisManager->FillNtupleDColumn(0, 3, y0);
               analysisManager->FillNtupleDColumn(0, 4, z0);
               analysisManager->FillNtupleDColumn(0, 5, root_dir2(0));
               analysisManager->FillNtupleDColumn(0, 6, root_dir2(1));
               analysisManager->FillNtupleDColumn(0, 7, root_dir2(2));
               analysisManager->AddNtupleRow(0);
           }

       }
       else
       {
            G4cout<<"Just a warning for root_weight <= 0."<<G4endl;
       } // add "abort event" ?
        outf.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
