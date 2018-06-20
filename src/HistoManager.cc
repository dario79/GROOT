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
/// \file hadronic/Hadr03/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 72245 2013-07-12 08:51:51Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Hadr03")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  //  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
  // Define histograms start values
//  const G4int kMaxHisto = 13;
  const G4int kMaxHisto = 43;  
  const G4String id[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
                         "20","21","22","23","24","25","26","27","28","29",  
                         "30","31","32","33","34","35","36","37","38","39",
                         "40","41","42"};
  const G4String title[] = 
                { "dummy",                                              //0
                  "Q = Ekin out - Ekin in",                            	//1
                  "Pbalance = mag(P_out - P_in)",              		//2
                  "kinetic energy of scattered primary particle",   	//3
                  "kinetic energy of recoil nuclei",                	//4
                  "kinetic energy of gamma",                        	//5
                  "kinetic energy of neutrons",                     	//6
                  "kinetic energy of protons",                      	//7
                  "kinetic energy of deuterons",                    	//8
                  "kinetic energy of alphas",                       	//9
                  "kinetic energy of all others ions",              	//10
                  "kinetic energy of all others mesons",            	//11
                  "kinetic energy of all others baryons",          	//12
                  "X momentum of scattered primary particle",   	//13
                  "X momentum of recoil nuclei",            		//14
                  "X momentum of gamma",                        	//15
                  "X momentum of neutrons",                     	//16
                  "X momentum of protons",                      	//17
                  "X momentum of deuterons",                    	//18
                  "X momentum of alphas",                       	//19
                  "X momentum of all others ions",              	//20
                  "X momentum of all others mesons",            	//21
                  "X momentum of all others baryons",           	//22
                  "Y momentum of scattered primary particle",   	//23
                  "Y momentum of recoil nuclei",                	//24
                  "Y momentum of gamma",                        	//25
                  "Y momentum of neutrons",                     	//26
                  "Y momentum of protons",                      	//27
                  "Y momentum of deuterons",                    	//28
                  "Y momentum of alphas",                       	//29
                  "Y momentum of all others ions",              	//30
                  "Y momentum of all others mesons",            	//31
                  "Y momentum of all others baryons",           	//32
                  "Z momentum of scattered primary particle",   	//33
                  "Z momentum of recoil nuclei",                	//34
                  "Z momentum of gamma",                        	//35
                  "Z momentum of neutrons",                     	//36
                  "Z momentum of protons",                      	//37
                  "Z momentum of deuterons",                    	//38
                  "Z momentum of alphas",                       	//39
                  "Z momentum of all others ions",              	//40
                  "Z momentum of all others mesons",            	//41
                  "Z momentum of all others baryons",           	//42
                 };    
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
//    G4int ih1 = analysisManager->CreateH1(id[k]+100, title1[k], nbins, vmin, vmax);
//    G4int ih2 = analysisManager->CreateH1(id[k]+200, title2[k], nbins, vmin, vmax);
//    G4int ih3 = analysisManager->CreateH1(id[k]+300, title3[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
//    analysisManager->SetH1Activation(ih1, false);
//    analysisManager->SetH1Activation(ih2, false);
//    analysisManager->SetH1Activation(ih3, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
