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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction4.cc
/// \brief Implementation of the PrimaryGeneratorAction4 class
//
//
// $Id: PrimaryGeneratorAction4.cc 83872 2014-09-20 22:23:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction4.hh"
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

PrimaryGeneratorAction4::PrimaryGeneratorAction4(G4ParticleGun* gun)
: fParticleGun(gun)
{
  //solid angle
  //

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction4::~PrimaryGeneratorAction4()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction4::GeneratePrimaries(G4Event* anEvent, G4double rbeam, G4double g_Energy, 
           G4String outputType, G4String outputName, DetectorConstruction* det)
{

    std::ofstream outf;
    outf.open(outputName+".txt", std::ios::out | std::ios::app );

//    G4ThreeVector oldPosition = G4ThreeVector(0.*cm,-38.*cm,0.); // MOVE THE SOURCE (p-a)
    G4ThreeVector oldPosition = G4ThreeVector(0.*cm,0.*cm,0.); // MOVE THE SOURCE (p-a)
    //rbeam = 0.75*cm; // 0.5*mm; /* sin(targetrot*deg)* */

    G4double r = sqrt(G4UniformRand());
    G4double theta = 2.*pi*G4UniformRand()*rad;
    G4double x0 = oldPosition.x(); //+ r*rbeam*cos(theta);
    G4double y0 = oldPosition.y(); //+ r*rbeam*sin(theta);

    G4double z0 = oldPosition.z(); // + (2*G4UniformRand()-1.)*0.00057537 *mm;

   G4double t0 = 0.1*ns;

   G4double MinTheta = 0*deg;//-10*deg;//0.*deg;//; TMath::ATan2(-4.,52.); // 0.07692307692307692308
   G4double MaxTheta = 180*deg;//10*deg;//180.*deg; //  TMath::ATan2(4.,52.); // 0.07692307692307692308
//  G4double rndm = G4UniformRand();
//  G4double costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
//  G4double sintheta = std::sqrt(1. - costheta*costheta);

//  rndm = G4UniformRand();
/*   G4double phi = 2.*pi*G4UniformRand()*rad; // MinPhi + (MaxPhi - MinPhi) * rndm;
   G4double pz = sin(theta) * cos(phi);
   G4double px = sin(theta) * sin(phi);//py = sin(theta) * sin(phi);
   G4double py = -cos(theta);//px = cos(theta);
*/

// from Hadr06
  G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double px = sinTheta*std::cos(phi),
           py = sinTheta*std::sin(phi),
           pz = cosTheta;

   px = -1; py = 0, pz = 0; // towards the Scionix

//  G4double Eg1 = 0.; // g_Energy*MeV; // *G4UniformRand();
// MAXWELLIAN EXTRACTION WITH T= g_Energy = 1.592 MeV ; n(E) = 2.146*sqrt(E)*exp(-E/1.592); --> normalizzata a 3.82 neutroni (integrale = 3.82) Cf source!
// x = -1.18*T*ln(1-r^(2/3)) -------- (-1.18*1.592*TMath::Log(1-TMath::Power(xx/1000.,2./3.)));

//   Eg1 =(-1.18*1.592*log(1-pow(G4UniformRand(),2./3.)));
   G4double  Eg1 =(-1.18*1.592*log(1-pow(G4UniformRand(),2./3.)));
//   G4double  Eg1 = 100.*keV;

   G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
   fParticleGun->SetParticleDefinition(particle);
   fParticleGun->SetParticleTime(t0);
   fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
   fParticleGun->SetParticleEnergy(Eg1*MeV);
   fParticleGun->GeneratePrimaryVertex(anEvent);

            if(outputType != "root")
            {

                outf  << "Generated Primary: "
                    << particle->GetParticleName() << " \t"
                    << std::setprecision(5) << Eg1*MeV<< " \t"
                    << std::setprecision(5) << x0 << " \t"
                    << std::setprecision(5) << y0 << " \t"
                    << std::setprecision(5) << z0 << " \t"
                    << std::setprecision(5) << px << " \t"
                    << std::setprecision(5) << py << " \t"
                    << std::setprecision(5) << pz << " \t"
                    << G4endl;
            }

            if(outputType != "txt")
            {
                G4String name = particle->GetParticleName();
                auto analysisManager = G4AnalysisManager::Instance();
                analysisManager->FillNtupleSColumn(0, 0, name);
                analysisManager->FillNtupleDColumn(0, 1, Eg1*MeV);
                analysisManager->FillNtupleDColumn(0, 2, x0);
                analysisManager->FillNtupleDColumn(0, 3, y0);
                analysisManager->FillNtupleDColumn(0, 4, z0);
                analysisManager->FillNtupleDColumn(0, 5, px);
                analysisManager->FillNtupleDColumn(0, 6, py);
                analysisManager->FillNtupleDColumn(0, 7, pz);
                analysisManager->AddNtupleRow(0);
            }

//  rndm = G4UniformRand();
/*   phi = 2.*pi*G4UniformRand()*rad; // MinPhi + (MaxPhi - MinPhi) * rndm;

//  G4double sinphi = std::sin(Phi);
//  G4double cosphi = std::cos(Phi);

   theta = MinTheta - G4UniformRand()*(MinTheta - MaxTheta);
//   phi = 2.*TMath::Pi()*G4UniformRand();

   px = sin(theta) * cos(phi*rad);
   py = sin(theta) * sin(phi);
   pz = cos(theta);
*/

   px = 1; py = 0, pz = 0; // towards the BaF

   G4double Eg2 = 10.;

   t0 = 0.*ns;
   G4ParticleDefinition* particle2 = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
   fParticleGun->SetParticleDefinition(particle2);
   fParticleGun->SetParticleTime(t0);
   fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
   fParticleGun->SetParticleEnergy(Eg2*MeV);
   fParticleGun->GeneratePrimaryVertex(anEvent);

            if(outputType != "root")
            {

                outf  << "Generated Primary: "
                    << particle2->GetParticleName() << " \t"
                    << std::setprecision(5) << Eg2*MeV<< " \t"
                    << std::setprecision(5) << x0 << " \t"
                    << std::setprecision(5) << y0 << " \t"
                    << std::setprecision(5) << z0 << " \t"
                    << std::setprecision(5) << px << " \t"
                    << std::setprecision(5) << py << " \t"
                    << std::setprecision(5) << pz << " \t"
                    << G4endl;
            }

            if(outputType != "txt")
            {
                G4String name = particle2->GetParticleName();
                auto analysisManager = G4AnalysisManager::Instance();
                analysisManager->FillNtupleSColumn(0, 0, name);
                analysisManager->FillNtupleDColumn(0, 1, Eg2*MeV);
                analysisManager->FillNtupleDColumn(0, 2, x0);
                analysisManager->FillNtupleDColumn(0, 3, y0);
                analysisManager->FillNtupleDColumn(0, 4, z0);
                analysisManager->FillNtupleDColumn(0, 5, px);
                analysisManager->FillNtupleDColumn(0, 6, py);
                analysisManager->FillNtupleDColumn(0, 7, pz);
                analysisManager->AddNtupleRow(0);
            }

   G4double Eg3 = 10.;
   px = -1; py = 0, pz = 0; // towards the BaF
   t0 = 0.*ns;

   G4ParticleDefinition* particle3 = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
   fParticleGun->SetParticleDefinition(particle3);
   fParticleGun->SetParticleTime(t0);
   fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
//  fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
   fParticleGun->SetParticleEnergy(Eg3*MeV);
   fParticleGun->GeneratePrimaryVertex(anEvent);

            if(outputType != "root")
            {

                outf  << "Generated Primary: "
                    << particle2->GetParticleName() << " \t"
                    << std::setprecision(5) << Eg3*MeV<< " \t"
                    << std::setprecision(5) << x0 << " \t"
                    << std::setprecision(5) << y0 << " \t"
                    << std::setprecision(5) << z0 << " \t"
                    << std::setprecision(5) << px << " \t"
                    << std::setprecision(5) << py << " \t"
                    << std::setprecision(5) << pz << " \t"
                    << G4endl;
            }

            if(outputType != "txt")
            {
                G4String name = particle2->GetParticleName();
                auto analysisManager = G4AnalysisManager::Instance();
                analysisManager->FillNtupleSColumn(0, 0, name);
                analysisManager->FillNtupleDColumn(0, 1, Eg3*MeV);
                analysisManager->FillNtupleDColumn(0, 2, x0);
                analysisManager->FillNtupleDColumn(0, 3, y0);
                analysisManager->FillNtupleDColumn(0, 4, z0);
                analysisManager->FillNtupleDColumn(0, 5, px);
                analysisManager->FillNtupleDColumn(0, 6, py);
                analysisManager->FillNtupleDColumn(0, 7, pz);
                analysisManager->AddNtupleRow(0);
            }


   outf.close();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
