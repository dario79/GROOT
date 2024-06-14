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
/// \file eventgenerator/particleGun/src/PrimaryGeneratorAction0.cc
/// \brief Implementation of the PrimaryGeneratorAction0 class
//
//
// $Id: PrimaryGeneratorAction0.cc 83872 2014-09-20 22:23:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction0.hh"
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

PrimaryGeneratorAction0::PrimaryGeneratorAction0(G4ParticleGun* gun)
: fParticleGun(gun)
{
  //solid angle
  //

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction0::~PrimaryGeneratorAction0()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction0::GeneratePrimaries(G4Event* anEvent, G4double rbeam, G4double g_Energy, 
           G4String outputType, G4String outputName, DetectorConstruction* det)
{

    std::ofstream outf;
    outf.open(outputName+".txt", std::ios::out | std::ios::app );

       G4ThreeVector oldPosition = G4ThreeVector(0.,0.,0.);
        //rbeam = 0.75*cm; // 0.5*mm; /* sin(targetrot*deg)* */

        G4double r = sqrt(G4UniformRand());
        G4double theta = 2.*pi*G4UniformRand()*rad;
        G4double x0 = oldPosition.x() + r*rbeam*cos(theta);
        G4double y0 = oldPosition.y() + r*rbeam*sin(theta);
        G4double z0 = oldPosition.z()-50.*cm; // + (2*G4UniformRand()-1.)*0.00057537 *mm;
        G4double px = 0.;
        G4double py = 0.;
        G4double pz = 1.;

        fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

        G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
        fParticleGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
        fParticleGun->SetParticleEnergy(g_Energy*MeV);


   fParticleGun->GeneratePrimaryVertex(anEvent);


            if(outputType != "root")
            {

                outf  << "Generated Primary: "
                    << particle->GetParticleName() << " \t"
                    << std::setprecision(5) << g_Energy*MeV<< " \t"
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
                analysisManager->FillNtupleDColumn(0, 1, g_Energy*MeV);
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
