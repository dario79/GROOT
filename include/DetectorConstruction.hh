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
/// \file electromagnetic/TestEm1/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 66586 2012-12-21 10:48:39Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

// CADMESH //
#include "CADMesh.hh"


class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     virtual G4VPhysicalVolume* Construct();
     
     void SetDetThickness(G4double);
     void SetDetMaterial (G4String);
     void SetMagField (G4double);
     void SetPath     (G4String);

     void SetBackingThickness(G4double);
     void SetBackingMaterial (G4String);
     void SetTargetPhi (G4double);
     void SetBackingPos(G4String);
     void SetReactionType(G4int);

     void UpdateGeometry();
    
    G4Material *targetLayerMaterial;
    G4Material *backingLayerMaterial;
    G4Material *worldMaterial;
    G4Material *ChamberMaterial;
    G4Material  *aluminium;
    G4Material *DetectorMaterial;
    G4Material *PTMaterial;
    G4Material* windowMaterial;
    G4Material* AirMaterial; 

  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return fPBox;};           
                    
     G4double           GetDetSize()       {return fDetThickness;};
     G4double           GetTargetPhi()     {return targetPhi;}
     G4int           GetReactionType()     {return fuffa;}
     G4Material*        GetDetMaterial()   {return fMaterial;};
     G4double           GetBackSize()       {return fBackingThickness;};
     G4Material*        GetBackMaterial()   {return fBackingMaterial;};
     
     void               PrintParameters();
     G4double targetLayerXPosition, targetLayerYPosition, targetLayerZPosition, targetLayerXDimension;
     G4double GetTargetPosZ(); 
     G4double GetTargetThickness();
     G4String detPath;

  private:

    G4VPhysicalVolume*    fPBox;
    G4LogicalVolume*      fLBox;
    G4LogicalVolume*      logictargetLayer;
     
    G4double              fDetThickness;
    G4double              fBackingThickness;
    G4double              targetPhi;
    G4int                 fuffa;
    G4String              backingPos;
    G4Material*           fMaterial;
    G4Material*           fBackingMaterial;
    G4UniformMagField*    fMagField;
     
    DetectorMessenger* fDetectorMessenger;
//    G4VSolid * cad_solid, *cad_solid2, *cad_solid3, *cad_solidb, *cad_solid4, *cad_solid5;
//    G4LogicalVolume * cad_logical, *cad_logical2, *cad_logical3, *cad_logicalb, *cad_logical4, *cad_logical5;
//    G4VPhysicalVolume * cad_physical, *cad_physical2, *cad_physical3, *cad_physicalb, *cad_physical4, *cad_physical5;

//    G4VSolid /* *SX3_solid[140],*/ * QQ3_solid[500];
    G4LogicalVolume /* *SX3_logical[140],*/ * QQ3_logical[500];
    G4VPhysicalVolume /* *SX3_physical[140],*/ * QQ3_physical[500];
//    CADMesh /* *SX3mesh[140],*/ * QQ3mesh[500];
//    CADMesh * mesh2, * mesh3, * mesh4, * mesh5;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

