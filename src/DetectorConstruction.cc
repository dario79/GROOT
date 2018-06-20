
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
/// \file electromagnetic/TestEm1/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 68245 2013-03-19 18:51:00Z maire $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"


  int verblevel = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fPBox(0), fLBox(0), fMaterial(0), fMagField(0)
{
    fBoxSize = 0.1*mm;
    DefineMaterials();
    SetMaterial("G4_Galactic");
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    //
    // define Elements
    //
    G4double z,a;

//    G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
//  G4Element* C  = new G4Element("Carbonium" ,"C" , z= 6., a=  12.00*g/mole);
    G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
    G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
//    G4Element* Mg = new G4Element("Magnesium"   ,"Mg" , z= 12., a=  24.305*g/mole);
//    G4Element* Ge = new G4Element("Germanium","Ge", z=32., a=  72.59*g/mole);
//    G4Element* Bi = new G4Element("Bismuth"  ,"Bi", z=83., a= 208.98*g/mole);
  
    //
    // define materials
    //
    G4double density;
    G4int ncomponents, natoms;
    G4double fractionmass;
    
    G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
    Air->AddElement(N, fractionmass=70.*perCent);
    Air->AddElement(O, fractionmass=30.*perCent);
/*    
    G4Material* H2l =
    new G4Material("H2liquid", density= 70.8*mg/cm3, ncomponents=1);
    H2l->AddElement(H, fractionmass=1.);
    
    G4Material* H2O =
    new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);
    ///H2O->SetChemicalFormula("H_2O");
    H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
    
    density = 0.001*mg/cm3;
    G4Material* CO2 = new G4Material("CO2", density, ncomponents=2);
    CO2->AddElement(C, natoms=1);
    CO2->AddElement(O, natoms=2);
*/   
//    new G4Material("D2_gas", z=2., a= 2.0141*g/mole, density= 0.036*mg/cm3);
    
//    new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
    
      aluminium = new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);
    
      G4Material  *silicon = new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
    
//    new G4Material("Chromium"   , z=24., a= 51.99*g/mole, density= 7.140*g/cm3);
    
//    new G4Material("Germanium"  , z=32., a= 72.61*g/mole, density= 5.323*g/cm3);
    
    G4Material *carbonMaterial = new G4Material("Carbonium", z=6., a=12.017*g/mole, density= 2.253*g/cm3);  // ATTENZIONE AL FATTORE 1000 (g->kg)
    G4Material* magnesiumMaterial = new G4Material("Magnesium", z=12., a=24.305*g/mole, density= 1.737*g/cm3);  // NELLA DENSITA'
    G4Material* lithiumMaterial = new G4Material("lithium", z=3., a=6.94*g/mole, density= 0.534*g/cm3);
    G4Element* elO  = new G4Element("Oxygen", "O", z= 8., a=16.00*g/mole);

    G4Material* Sn112Material = new G4Material("Sn112", z=50., a=111.904818*g/mole, density= 7.28*g/cm3);

//    G4Material* carbonMaterial = new G4Material("Ruthenium", z=44.,  a=101.07*g/mole, density= 12.1 *g/cm3);  

/*    G4Material* BGO = new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
    BGO->AddElement(O , natoms=12);
    BGO->AddElement(Ge, natoms= 3);
    BGO->AddElement(Bi, natoms= 4);
*/    
//    G4Material  *iron = new G4Material("Iron", z=26., a= 55.85*g/mole, density= 7.870*g/cm3);
    
//    new G4Material("Tungsten"   , z=74., a=183.85*g/mole, density= 19.30*g/cm3);
    
    G4Material *AuMaterial = new G4Material("Gold", z=79., a=196.97*g/mole, density= 19.32*g/cm3);
    
    G4Material *PbMaterial = new G4Material("Lead", z=82., a=207.19*g/mole, density= 11.35*g/cm3);
    
//    new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);
    
    
//    G4Material* argonGas = new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,kStateGas, 273.15*kelvin, 1*atmosphere);
    
/*    G4Material* butane = new G4Material("Isobutane",density= 2.42*mg/cm3, ncomponents=2,kStateGas,273.15*kelvin, 1*atmosphere);
    butane->AddElement(C, natoms=4);
    butane->AddElement(H, natoms=10);
    
    G4Material* ArButane =
    new G4Material("ArgonButane", density= 1.835*mg/cm3, ncomponents=2,
                   kStateGas,273.15*kelvin,1.*atmosphere);
    ArButane->AddMaterial(argonGas, fractionmass=70*perCent);
    ArButane->AddMaterial(butane ,  fractionmass=30*perCent);

*/


    // MATERIAL DEFINIOTION VIA THE NIST MATERIAL BUILDER
    G4bool isotopes = false;
//    G4Material* Graphite = G4NistManager::Instance()->FindOrBuildMaterial("G4_GRAPHITE", isotopes);
//    G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
// Vacuum
  G4int nElem;
  density = 6.078e-12 *g/cm3; ;//2.376e-15 *g/cm3; 
  G4double temperature = 300. *kelvin;
  G4double pressure = 1.e-8 *bar;//; 2.e-7 *bar; 
  G4Material* Vacuum = new G4Material("Vacuum", density, nElem=1,kStateGas,temperature,pressure);
  Vacuum->AddMaterial(Air,1.);

        G4int zz;
        G4NistManager* pMatMan = G4NistManager::Instance();
        G4Element* elCr = pMatMan->FindOrBuildElement(zz=24);
//	G4Element* elMn = G4NistManager::Instance()->FindOrBuildElement(zz=25);
	G4Element* elFe = pMatMan->FindOrBuildElement(zz=26);
//	G4Element* elCo = G4NistManager::Instance()->FindOrBuildElement(zz=27);
	G4Element* elNi = pMatMan->FindOrBuildElement(zz=28);
        G4Element* elC =  pMatMan->FindOrBuildElement(zz=6);



//	G4double a,z;
//       G4double weightRatio;

/* 	G4String name,symbol;

        G4Element* elCr= new G4Element(name="Chromium", symbol="Cr", z=24., a = 51.9961*g/mole);
	G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z=25., a = 54.94*g/mole);
	G4Element* elFe= new G4Element(name="Iron", symbol="Fe", z=26., a = 55.845*g/mole);
	G4Element* elNi= new G4Element(name="Nickel", symbol="Ni", z=28., a = 58.6934*g/mole);
	G4Element* elCu= new G4Element(name="Copper", symbol="Cu", z=29., a = 63.55*g/mole);
	G4Element* elMo = new G4Element(name="Molydb", symbol="Mo", z=42., a = 95.96*g/mole);
*/
	//StainlessSteel, Material Names : stainless steel
	//Material : Fe-Cr-Ni-Mo,  93xx: 3.25% Ni,1.2%Cr,0.12%Mo,

/* OLD	G4Material* stainlesssteel = new G4Material("StainlessSteel", density = 7.85 *g/cm3, nElem=3);
	stainlesssteel->AddElement(elFe, weightRatio=0.9550);
	stainlesssteel->AddElement(elCr, weightRatio=0.0330);
	stainlesssteel->AddElement(elNi, weightRatio=0.0120);

C: 0,04%
Cr: 18,2%
Ni: 8,1%
Fe: 73,66%

densità 8,00 g/cm3

*/    

       G4Material*  stainlesssteel = new G4Material("StainlessSteel", density = 8.0 *g/cm3, nElem=4);
	    stainlesssteel->AddElement(elFe, 0.7366);
	    stainlesssteel->AddElement(elCr, 0.1820);
	    stainlesssteel->AddElement(elNi, 0.0810);
	    stainlesssteel->AddElement(elC, 0.0004);


    G4Element* Zn = pMatMan->FindOrBuildElement(30);
    G4Element* Cu = pMatMan->FindOrBuildElement(29);
    G4Material* Brass = new G4Material("Brass", density = 8.70*g/cm3, nElem=2);
    Brass->AddElement(Zn, 3);
    Brass->AddElement(Cu, 7);

    G4Element* elLi = pMatMan->FindOrBuildElement(zz=3);
    G4Element* elH = pMatMan->FindOrBuildElement(zz=1);
    G4Material* LiH = new G4Material("LiH", density = 0.82*g/cm3, nElem=2);
    LiH->AddElement(elLi, 1);
    LiH->AddElement(elH, 1);
    G4Element* elF = pMatMan->FindOrBuildElement(zz=9);
    G4Material* LiF = new G4Material("LiF", density = 2.64*g/cm3, nElem=2);
    LiF->AddElement(elLi, 1);
    LiF->AddElement(elF, 1);

// Mylar

   density = 1.39*g/cm3;
   G4Material* Mylar = new G4Material("Mylar", density, nElem=3);
   Mylar->AddElement(elO,2);
   Mylar->AddElement(elC,5);
   Mylar->AddElement(elH,4);

// window? 

   G4Element* elSi = pMatMan->FindOrBuildElement(14);
   density = 2.200*g/cm3;
   G4Material* SiO2 = new G4Material("quartz", density, nElem=2);
   SiO2->AddElement(elSi,1);
   SiO2->AddElement(elO,2);

//Havar
  //
  G4Element* elCo = pMatMan->FindOrBuildElement(zz=27);
  G4Element* elW  = pMatMan->FindOrBuildElement(zz=74);

  G4Material* Havar = new G4Material("Havar", density= 8.3*g/cm3, nElem=5);
  Havar->AddElement(elCr, 0.1785);
  Havar->AddElement(elFe, 0.1822);
  Havar->AddElement(elCo, 0.4452);
  Havar->AddElement(elNi, 0.1310);
  Havar->AddElement(elW , 0.0631);

   G4Element* elTi = pMatMan->FindOrBuildElement(zz=22);

   G4Material* Titanium = new G4Material("Titanium", z=22, a=47.867*g/mole, density= 4.506*g/cm3);

// Kapton Dupont de Nemur (density: 1.396-1.430, get middle )
 
  G4Element* elN  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  density = 1.413*g/cm3;
  G4Material* Kapton = new G4Material("Kapton", density, nElem=4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

    PTMaterial = Brass;
//    targetLayerMaterial = magnesiumMaterial;
    targetLayerMaterial = Sn112Material;//; LiF; // lithiumMaterial; // LiH; // lithiumMaterial;
//    targetLayerMaterial = PbMaterial;
    backingLayerMaterial =  Mylar; // AuMaterial; // carbonMaterial;
    worldMaterial = Vacuum;
    AirMaterial = Air;
    ChamberMaterial = aluminium; //stainlesssteel;
//    ChamberMaterial = stainlesssteel;
    windowMaterial = aluminium; // tainlesssteel; // Kapton; //aluminium; // Titanium;// SiO2; // Havar; // 
//    VoidOrbMaterial = Vacuum;
/*    SX3Material = silicon;*/
//    G4cout << " Now.."<<G4endl;
    DetectorMaterial = silicon;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
//    G4cout << " .. and then."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

//    G4Colour  cyan (0.0, 1.0, 1.0);
  G4VisAttributes *theWhite = new G4VisAttributes( G4Colour(255/255., 255/255., 255/255. ));
    theWhite -> SetVisibility(true);
    theWhite -> SetForceSolid(true);
//    theWhite -> SetForceWireframe(true);

 
    G4VisAttributes *darkGreen = new G4VisAttributes( G4Colour(0/255., 100/255., 0/255. ));
    darkGreen -> SetVisibility(true);
    darkGreen -> SetForceSolid(true);
//    darkGreen -> SetForceWireframe(true);

/*    G4VisAttributes *red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
    red -> SetVisibility(true);
    red-> SetForceSolid(true);
//    red -> SetForceWireframe(true);
*/
    G4VisAttributes *cyan = new G4VisAttributes(G4Colour(0/255., 255/255., 255/255.));
    cyan -> SetVisibility(true);
//    cyan-> SetForceSolid(true);
    cyan -> SetForceWireframe(true);
    cyan->SetForceAuxEdgeVisible(true);

    G4VisAttributes *yellow = new G4VisAttributes(G4Colour(255/255., 255/255. ,0/255.));
    yellow -> SetVisibility(true);
//    yellow-> SetForceSolid(true);
    yellow -> SetForceWireframe(true);
    yellow->SetForceAuxEdgeVisible(true);
    yellow-> SetForceLineSegmentsPerCircle(64); //  16 segments of 4 lines each


    G4VisAttributes *red = new G4VisAttributes(G4Colour(255/255., 0/255., 0/255.));
    red -> SetVisibility(true);
//    red-> SetForceSolid(true);
    red -> SetForceWireframe(true);
    red->SetForceAuxEdgeVisible(true);
    red-> SetForceLineSegmentsPerCircle(64); //  16 segments of 4 lines each

    G4double worldVolumeXDimension = 3.0 *m;
    //G4double worldVolumeXDimension = 10.0 *mm;
    G4double worldVolumeYDimension = 3.0 *m;
    G4double worldVolumeZDimension = 5.0 *m;

//    G4double chamberRmin = 200. *mm; // 1000
//    G4double chamberRmax = 210. *mm; // 1010


    G4double gapBeforeTarget = 0.0001 *mm; // 0.01 *mm;
    
// Au 200.*nm; // 60ug/cm2 di 12C 0.000266 *mm; 
    G4double backingLayerXDimension = 0.0013*mm; //200.*nm;
//    G4double backingLayerXDimension = 0.60 *mm; // *density/section 
//    G4double backingLayerYDimension = 10.0 *mm;
//    G4double backingLayerZDimension = 10.0 *mm;

// 200ug/cm2 LiH 0.0024390243902439024 *mm; // 100ug/cm2 24Mg 0.000575 *mm; // 0.1g/cm2 Pb 0.088105726872246696 *mm // 100ug/cm2 7Li 0.0018726592 *mm
// 200ug/cm2 LiF 0.00075757575.. *mm; 
// 112Sn 350ug/cm2 --> 4.8e-5 m = 48um
/*    G4double */ targetLayerXDimension = 50.*um; // 0.0018726592 *mm;
//    G4double targetLayerXDimension = 0.100 *mm; // *density/section 
//    G4double targetLayerYDimension = 10.0 *mm;
//    G4double targetLayerZDimension = 10.0 *mm;

/*    G4double */ targetLayerXPosition = 0. *mm;
/*    G4double */ targetLayerYPosition = 0. *mm;
/*    G4double */ targetLayerZPosition = /* gapBeforeTarget + */ 0. *mm ;

    G4double backingLayerXPosition = 0. *mm;
    G4double backingLayerYPosition = 0. *mm;    
    G4double backingLayerZPosition = gapBeforeTarget +  backingLayerXDimension/2. + targetLayerXDimension/2.;
    G4double backingLayerZPosition2 = -gapBeforeTarget - backingLayerXDimension/2. -targetLayerXDimension/2.;

    // WORLD VOLUME
    G4Box* sBox = new G4Box("Container",
                            worldVolumeXDimension/2,worldVolumeYDimension/2,worldVolumeZDimension/2);          
    
    
    fLBox = new G4LogicalVolume(sBox,
                                worldMaterial,            
                                "WorldLogical");
    
    fPBox = new G4PVPlacement(0,                          //no rotation
                              G4ThreeVector(0,0,0),            //at (0,0,0)
                              fLBox,                       //its logical volume
                              "WorldPhys",
                              0,                           //its mother  volume
                              false,                       //no boolean operation
                              0);                          //copy number
    
   fLBox->SetVisAttributes(G4VisAttributes::Invisible);

   G4RotationMatrix* rotm  = new G4RotationMatrix();
//   rotm->rotateY(90.*deg);
//   rotm->rotateY(2.*M_PI);
/*
   // Chamber ? al
   G4Sphere* solidChamber =  new G4Sphere("Chamber",chamberRmin,chamberRmax,0.*deg,360.*deg,5.*deg,170.*deg);
   G4LogicalVolume* logicChamber = new G4LogicalVolume(solidChamber,ChamberMaterial,"ChamberLogic");
   G4VPhysicalVolume* physiChamber = new G4PVPlacement(rotm,G4ThreeVector(0.,0.,0.),logicChamber,"ChamberPhys",fLBox,false,0);
   // Visualisation of the External part
   logicChamber-> SetVisAttributes(cyan);
*/
   // void inside Chamber? vuoto
/*   G4Sphere* solidVoidOrb =  new G4Sphere("VoidOrb",0.1*m,0.7*m,0.0*deg,360.*deg,0.*deg,180.*deg);
   G4LogicalVolume* logicVoidOrb = new G4LogicalVolume(solidVoidOrb,VoidOrbMaterial,"VoidOrbLogic");
   G4VPhysicalVolume* physiVoidOrb = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(0.,0.,0.)),logicVoidOrb,"VoidOrbPhys",fLBox,false,0);
 // new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicVoidOrb,"VoidOrbPhys",fLBox,false,0);
   // Visualisation of the External part
   logicVoidOrb-> SetVisAttributes(yellow);
//   logicVoidOrb->SetVisAttributes(G4VisAttributes::Invisible);
*/

//  Chamber "a botte"

//    G4double ChamberOffSet = 10. *cm;

// CILINDRO "pieno di vuoto"
//    G4Tubs* solidVoidOrb = new G4Tubs("VoidOrb",0.2*m,1.*m,1.*m,0.*deg,360.*deg);
//    G4LogicalVolume* logicVoidOrb = new G4LogicalVolume(solidVoidOrb,VoidOrbMaterial,"VoidOrbLogic");
//    G4VPhysicalVolume* physiVoidOrb = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicVoidOrb,"VoidOrbPhys",fLBox,false,0);
//    logicVoidOrb-> SetVisAttributes(yellow);
//    logicVoidOrb->SetVisAttributes(G4VisAttributes::Invisible);



// 8.8cm radius (cathetus) ~4.2cm (cathetus)  angle =+/-
// ip = 8.8./cos t = 4.2/sin t   t = arctg 4.2/8.8  = 25.4
   G4double  targetrot = 0.; // 35.; // in deg
//   rotm->rotateY(90.*deg);

    // TARGET: 12C Layer
//    G4Box* solidbackingLayer = new G4Box("backingLayer",backingLayerXDimension/2,backingLayerYDimension/2,backingLayerZDimension/2);
/*    G4Tubs* solidbackingLayer = new G4Tubs("backingLayer",0.,0.99*cm,backingLayerXDimension/2.,0.*deg,360.*deg);
    G4LogicalVolume* logicbackingLayer = new G4LogicalVolume(solidbackingLayer,backingLayerMaterial,"backingLayerLogic");
    G4VPhysicalVolume* physibackingLayer = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,backingLayerZPosition)),logicbackingLayer,"backingLayerPhys",fLBox,false,0);
//    G4VPhysicalVolume* physibackingLayer = new G4PVPlacement(0,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,backingLayerZPosition),logicbackingLayer,"backingLayerPhys",fLBox,false,0);
     G4VisAttributes * backVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.,1.));
    backVisAtt -> SetForceSolid(true);
    logicbackingLayer -> SetVisAttributes(backVisAtt);
*/
/*    G4Tubs* solidbackingLayer2 = new G4Tubs("backingLayer2",0.,0.99*cm,backingLayerXDimension/2.,0.*deg,360.*deg);
    G4LogicalVolume* logicbackingLayer2 = new G4LogicalVolume(solidbackingLayer2,backingLayerMaterial,"backingLayerLogic2");
    G4VPhysicalVolume* physibackingLayer2 = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,backingLayerZPosition2)),logicbackingLayer2,"backingLayerPhys2",fLBox,false,0);
     G4VisAttributes * backVisAtt2 = new G4VisAttributes(G4Colour(0.,0.,1.,1.));
    backVisAtt2 -> SetForceSolid(true);
    logicbackingLayer2 -> SetVisAttributes(backVisAtt2);
*/
    // TARGET: 27Mg Layer    
//    G4Box* solidtargetLayer = new G4Box("targetLayer",targetLayerXDimension/2,targetLayerYDimension/2,targetLayerZDimension/2);
    G4Tubs* solidtargetLayer = new G4Tubs("targetLayer",0.,0.99*cm,targetLayerXDimension/2.,0.*deg,360.*deg);
    G4LogicalVolume* logictargetLayer = new G4LogicalVolume(solidtargetLayer,targetLayerMaterial,"targetLayerLogic");
//    G4VPhysicalVolume* physitargetLayer = new G4PVPlacement(0,G4ThreeVector(targetLayerXPosition,targetLayerYPosition,targetLayerZPosition),logictargetLayer,"targetLayerPhys",fLBox,false,0);
    G4VPhysicalVolume* physitargetLayer = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(targetLayerXPosition,targetLayerYPosition,targetLayerZPosition)),logictargetLayer,"targetLayerPhys",fLBox,false,0);
    logictargetLayer -> SetVisAttributes(theWhite);


/*
    G4Tubs* collimator = new G4Tubs("collimator",0.5*cm,10.*cm,25.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logiccollimator = new G4LogicalVolume(collimator,aluminium,"logiccollimator");
    G4VPhysicalVolume* physicollimator = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,100.*cm)),logiccollimator,"physicollimator",fLBox,false,0);
//     G4VisAttributes * collimatorVisAtt2 = new G4VisAttributes(G4Colour(0.,0.,1.,1.));
//    logiccollimator -> SetForceSolid(true);
    logiccollimator ->  SetVisAttributes(theWhite); // SetVisAttributes(G4Color(1., 1., 1., 0.6));
*/

    G4Tubs* window = new G4Tubs("window",0.*cm,2.*cm,/*0.1*cm/*0.0025*cm */1.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicwindow = new G4LogicalVolume(window,windowMaterial,"logicwindow");
    G4VPhysicalVolume* physiwindow = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,-127.5*cm)),logicwindow,"physiwindow",fLBox,false,0);
    logicwindow ->  SetVisAttributes(theWhite); // SetVisAttributes(G4Color(1., 1., 1., 0.6));

    G4Tubs* flange = new G4Tubs("flange",2.*cm,12.*cm,3.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicflange = new G4LogicalVolume(flange, aluminium,"logicflange");
    G4VPhysicalVolume* physiflange = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,-127.5*cm)),logicflange,"physiflange",fLBox,false,0);
    logicflange ->  SetVisAttributes(red); // SetVisAttributes(G4Color(1., 1., 1., 0.6));

    G4Tubs* theair = new G4Tubs("theair",0.*cm,5.*cm,5.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicair = new G4LogicalVolume(theair, AirMaterial,"logicair");
    G4VPhysicalVolume* physiair = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,-140.*cm)),logicair,"physiair",fLBox,false,0);
    logicair ->  SetVisAttributes(theWhite); // SetVisAttributes(G4Color(1., 1., 1., 0.6));



   // CAD model rotation. 
    G4RotationMatrix * rot = new G4RotationMatrix();
//    rot->rotateX(90*deg);
    rot->rotateY(0.*deg); // 90 or -90?

    G4ThreeVector Theoffset;
    // Load CAD file as tessellated solid //

// IN STP FILE: distance from reference frame origin to the center of the bottom target holder circle ( =====|O|O|O|O|0| ) is (524.5*mm, 5.*mm, 0.2*mm), with target holder along Z-axis and X the beam direction.

    G4double X_offset = 0.*cm; // 52.45*cm;       OLD -> //  -15.363*cm; // DET mid  -15.103*cm; SENSITIVE DET mid  -14.363*cm; 
    G4double Y_offset = 0.*cm; // -0.5*cm;
    G4double Z_offset =  0.*cm;
    Theoffset.set(Z_offset, Y_offset, X_offset); // FIX THIS
    
    // Note that offset is applied to the points in mesh directly before placement.

 G4int a=1, b=1;
/*
 for (G4int i = 0; i<140; i++){ 

    if(i>0 && i%4==0) { 
          if(a==12)a++; 
       a++;
       b = 1;
    }
    else if(i>0){
       b++;
    }
    G4String stripfilename = "/home/lattuadad/geom_import/M/SX3";
    stripfilename.append("-");
    stripfilename.append(std::to_string(a));
    stripfilename.append("-");
    stripfilename.append(std::to_string(b));
    stripfilename.append(".stl");

    G4cout<<stripfilename.c_str()<<G4endl;
    G4String stripID = "SX3-";
    stripID.append(std::to_string(a));
    stripID.append("-");
    stripID.append(std::to_string(b));

    G4String stripIDlog = stripID;
    stripIDlog.append("_logical");

//    G4String stripIDphys = stripID;
    G4String stripIDphys = "SX3-";
    G4String TheDetID = std::to_string(a);
    if(TheDetID.size()==1) {
         stripIDphys.append("0");
    }
    stripIDphys.append(TheDetID);
    stripIDphys.append("-");
    G4String TheStripID = std::to_string(b);
    if(TheStripID.size()==1) {
         stripIDphys.append("0");
    }
    stripIDphys.append(TheStripID);

    stripIDphys.append("_physical");

    SX3mesh[i] = new CADMesh( (char*) stripfilename.c_str(), (char*) "STL");
    SX3mesh[i]->SetScale(mm);
//    SX3mesh->SetOffset(offset);
    SX3mesh[i]->SetReverse(false);

    SX3_solid[i] = SX3mesh[i]->TessellatedMesh();
    SX3_logical[i] = new G4LogicalVolume(SX3_solid[i], DetectorMaterial, stripIDlog, 0, 0, 0);
//    SX3_physical[i] = new G4PVPlacement(rot, Theoffset, SX3_logical[i],
    SX3_physical[i] = new G4PVPlacement (G4Transform3D(*rot,Theoffset), SX3_logical[i],
                                     stripIDphys, fLBox, false, 0);
    SX3_logical[i]->SetVisAttributes(G4Color(0., 0.4, 0.4, 1.));
 }
*/
/*  
    CADMesh * mesh = new CADMesh((char*) "/home/lattuadad/geom_import/Final/SX3-1-1.stl", (char*) "STL");
    mesh->SetScale(mm);
//    mesh->SetOffset(offset);
    mesh->SetReverse(false);

    cad_solid = mesh->TessellatedMesh();
    cad_logical = new G4LogicalVolume(cad_solid, DetectorMaterial, "cad_logical", 0, 0, 0);
    cad_physical = new G4PVPlacement(rot, Theoffset, cad_logical,
                                     "cad_physical", fLBox, false, 0);
    cad_logical->SetVisAttributes(G4Color(0.0, 1.0, 0.0, 1.));

    CADMesh * meshb = new CADMesh((char*) "/home/lattuadad/geom_import/Final/SX3-1-2.stl", (char*) "STL");
    meshb->SetScale(mm);
//    meshb->SetOffset(offset);
    meshb->SetReverse(false);

    cad_solidb = meshb->TessellatedMesh();
    cad_logicalb = new G4LogicalVolume(cad_solidb, DetectorMaterial, "cad_logicalb", 0, 0, 0);
    cad_physicalb = new G4PVPlacement(rot, Theoffset, cad_logicalb,
                                     "cad_physicalb", fLBox, false, 0);
    cad_logicalb->SetVisAttributes(G4Color(0.0, 0.5, 0.5, 0.5));
*/
//    offset = G4ThreeVector(-300*cm, -300*cm, -300*cm);

 a = 1, b = 1;

    G4RotationMatrix  rot3[12]; // = new G4RotationMatrix(); // rot; // 
    G4ThreeVector Theoffset3(0.,0.,-8.8*cm); // there's a Z offset in the stl file.. I had to correct for that
   
 for (G4int i = 0; i<12; i++){ 
/*
    G4String stripfilename = "/home/lattuadad/geom_import/M/YY1";
    stripfilename.append("-");
    stripfilename.append(std::to_string(i));
//    stripfilename.append("-");
//    stripfilename.append(std::to_string(b));
    stripfilename.append(".stl");

    G4cout<<stripfilename.c_str()<<G4endl;
    G4String stripID = "YY1-";
    stripID.append(std::to_string(i));
//    stripID.append("-");
//    stripID.append(std::to_string(b));

    G4String stripIDlog = stripID;
    stripIDlog.append("_logical");

//    G4String stripIDphys = stripID;
    G4String stripIDphys = "YY1-";
    G4String TheDetID = std::to_string(i);
    if(TheDetID.size()==1) {
         stripIDphys.append("0");
    }
    stripIDphys.append(TheDetID);
/*    stripIDphys.append("-");
    G4String TheStripID = std::to_string(b);
    if(TheStripID.size()==1) {
         stripIDphys.append("0");
    }
    stripIDphys.append(TheStripID);
*/ /*
    stripIDphys.append("_physical");   

    if(i<6) {
         Theoffset3.set(Z_offset, Y_offset, X_offset-10.*cm);
    }
    else {
         Theoffset3.set(Z_offset, Y_offset, X_offset+10.*cm);
    }
    rot3[i].rotateY(90.*deg); // 90 or -90?
    rot3[i].rotateZ(i*60.*deg);
//    rot3[i].rotateX(47.*deg);

*/
    G4String stripfilename = "/home/lattuadad/geom_import/SIDAR/SIDAR";
    stripfilename.append(std::to_string(i+1));
//    stripfilename.append("-");
//    stripfilename.append(std::to_string(b));
    stripfilename.append(".stl");

    G4cout<<stripfilename.c_str()<<G4endl;
    G4String stripID = "YY1-";
    stripID.append(std::to_string(i));
//    stripID.append("-");
//    stripID.append(std::to_string(b));

    G4String stripIDlog = stripID;
    stripIDlog.append("_logical");

//    G4String stripIDphys = stripID;
    G4String stripIDphys = "YY1-";
    G4String TheDetID = std::to_string(i);
    if(TheDetID.size()==1) {
         stripIDphys.append("0");
    }
    stripIDphys.append(TheDetID);
    stripIDphys.append("_physical");
//    Theoffset3.set(Z_offset, Y_offset, X_offset);
    

//    G4cout<<stripfilename.c_str()<<G4endl;
    QQ3mesh[i] = new CADMesh( (char*) stripfilename.c_str(), (char*) "STL");
    QQ3mesh[i]->SetScale(mm);
//    QQ3mesh->SetOffset(offset);
    QQ3mesh[i]->SetReverse(false);

    QQ3_solid[i] = QQ3mesh[i]->TessellatedMesh();
    QQ3_logical[i] = new G4LogicalVolume(QQ3_solid[i], DetectorMaterial, stripIDlog, 0, 0, 0);
//    QQ3_physical[i] = new G4PVPlacement(rot, Theoffset, QQ3_logical[i],
    QQ3_physical[i] = new G4PVPlacement (G4Transform3D(rot3[i],Theoffset3), QQ3_logical[i],
                                     stripIDphys, fLBox, false, 0);
    QQ3_logical[i]->SetVisAttributes(G4Color(0., 1., 0., 1.));
 }

/*
     G4Colour  white   ()              ;  // white
     G4Colour  white   (1.0, 1.0, 1.0) ;  // white
     G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
     G4Colour  black   (0.0, 0.0, 0.0) ;  // black
     G4Colour  red     (1.0, 0.0, 0.0) ;  // red
     G4Colour  green   (0.0, 1.0, 0.0) ;  // green
     G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
     G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
     G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
     G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
*/
    // Note that offset is applied to the points in mesh directly before placement.
/*    CADMesh * */
/*    mesh3 = new CADMesh((char*) "/home/lattuadad/geom_import/M/YY1_1b.stl", (char*) "STL");
    mesh3->SetScale(mm);
//    mesh3->SetOffset(offset);
    mesh3->SetReverse(false);

    G4ThreeVector Theoffset3;
    Theoffset3.set(0.,0.,0.); // FIX THIS
//    Theoffset2.set(Z_offset, -0.5*cm-52.4*cm*(sin(targetrot*deg)), 52.4*cm*(cos(targetrot*deg))); // FIX THIS

    G4RotationMatrix * rot3 = rot; // new G4RotationMatrix();
    rot3->rotateY(90.*deg); // 90 or -90?

    Theoffset3.set(Z_offset, Y_offset, X_offset+10.*cm);
    cad_solid3 = mesh3->TessellatedMesh();
    cad_logical3 = new G4LogicalVolume(cad_solid3, DetectorMaterial, "cad_logical3", 0, 0, 0);
    cad_physical3 = new G4PVPlacement(// rot, Theoffset,
                                     G4Transform3D(*rot3,Theoffset3), cad_logical3,
                                     "cad_physical3", fLBox, false, 0);
    cad_logical3->SetVisAttributes(G4Color(0., 1., 0., 1.)); // ->SetVisAttributes(G4Color(1., 1., 1., 0.6));
*/
/*
    Theoffset.set(0., 0., 60.*cm +Z_offset);  // FIX THIS
    G4RotationMatrix * rot2 = new G4RotationMatrix();
    rot2->rotateX(90*deg);  // -90, 90, 90  and flip..
    rot2->rotateY(90*deg);
    rot2->rotateZ(90*deg);
*/

    // Note that offset is applied to the points in mesh directly before placement.
/*    CADMesh * */
    mesh2 = new CADMesh((char*) "/home/lattuadad/geom_import/M/Chamber2.stl", (char*) "STL");
    mesh2->SetScale(mm);
//    mesh2->SetOffset(offset);
    mesh2->SetReverse(false);

    cad_solid2 = mesh2->TessellatedMesh();
    cad_logical2 = new G4LogicalVolume(cad_solid2, aluminium, "cad_logical2", 0, 0, 0);
    cad_physical2 = new G4PVPlacement(G4Transform3D(*rot,Theoffset), cad_logical2,
                                     "cad_physical2", fLBox, false, 0);
    cad_logical2->SetVisAttributes(G4Color(1., 1., 1., 0.3));

/*    CADMesh * */ 
/*    mesh4 = new CADMesh((char*) "/home/lattuadad/geom_import/M/PCB.stl", (char*) "STL");
    mesh4->SetScale(mm);
//    mesh4->SetOffset(offset);
    mesh4->SetReverse(false);

    cad_solid4 = mesh4->TessellatedMesh();
    cad_logical4 = new G4LogicalVolume(cad_solid4, DetectorMaterial, "cad_logical4", 0, 0, 0);   // FIX THE MATERIAL? Fiberglass?
    cad_physical4 = new G4PVPlacement(//rot, Theoffset, cad_logical4,
                                     G4Transform3D(*rot,Theoffset), cad_logical4,
                                     "cad_physical4", fLBox, false, 0);
    cad_logical4->SetVisAttributes(G4Color(1., 1., 0., 1.));
*/
/*
//    Theoffset.set(-0.44*cm, (12.5+4.)*cm, 0.85*cm);// 15.*cm -Z_offset);  // FIX THIS
//    Theoffset.set(-7.75*cm, -10.4*cm, 0.*cm);// 15.*cm -Z_offset);  // FIX THIS
    Theoffset.set((-7.75+2.45*(sin(targetrot*deg)))*cm, -10.4*cm, sin(targetrot*deg)*7.75*cm);
//targetLayerXPosition+cos(targetrot*deg)*1.*cm), targetLayerYPosition, targetLayerZPosition-(sin(targetrot*deg)*1.*cm)
   G4RotationMatrix * rot3 = new G4RotationMatrix();
    rot3->rotateX(90*deg);
    rot3->rotateZ(90.*deg);
    rot3->rotateZ(-targetrot*deg);
*/
// /control/execute vis.mac



// QQ3 HERE


   // CAD model ladder (co)rotation. 
    G4ThreeVector Theoffset2;
    Theoffset2.set(0.,0.,0.); // FIX THIS
//    Theoffset2.set(Z_offset, -0.5*cm-52.4*cm*(sin(targetrot*deg)), 52.4*cm*(cos(targetrot*deg))); // FIX THIS

    G4RotationMatrix * rot2 = rot; // new G4RotationMatrix();
//    rot2->rotateY(90*deg); // 90 or -90?
//    rot2->rotateZ(90.*deg); // 90 or -90? 
    rot2->rotateY(targetrot*deg); // 90 or -90?


//    CADMesh 
    mesh5 = new CADMesh((char*) "/home/lattuadad/geom_import/M/TargetLadder.stl", (char*) "STL");
    mesh5->SetScale(mm);
//    mesh->SetOffset(offset);
    mesh5->SetReverse(false);
    cad_solid5 = mesh5->TessellatedMesh();
    cad_logical5 = new G4LogicalVolume(cad_solid5, aluminium, "cad_logical5", 0, 0, 0);
    cad_physical5 = new G4PVPlacement(//rot2, Theoffset, cad_logical5,
                                      G4Transform3D(*rot2,Theoffset2), cad_logical5,
                                     "cad_physical5", fLBox, false, 0);
    cad_logical5->SetVisAttributes(G4Color(1., 1., 1., 0.25));

   G4cout<<"Detector ready."<<G4endl;   

//
    PrintParameters();
   
    //always return the root volume
    //
    return fPBox;
//
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
    << " of " << fMaterial->GetName() <<" but who cares? "<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMaterial(G4String materialChoice)
{
    // search the material by its name
    ////G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
    
    if (pttoMaterial) { fMaterial = pttoMaterial;
    } else {
        G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
        << materialChoice << " not found" << G4endl;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSize(G4double value)
{
  if(verblevel >1)    G4cout << "Dentro Detector Contruction --> STAMPO IL VALORE DELLO SPESSORE DELLA FETTA ===>  " << value<< G4endl;
    fBoxSize = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
    //apply a global uniform magnetic field along Z axis
    G4FieldManager* fieldMgr
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
    if (fMagField) delete fMagField;        //delete the existing magn field
    
    if (fieldValue!=0.)                        // create a new one if non nul
    {
        fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
        fieldMgr->SetDetectorField(fMagField);
        fieldMgr->CreateChordFinder(fMagField);
    }
    else
    {
        fMagField = 0;
        fieldMgr->SetDetectorField(fMagField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}



G4double DetectorConstruction::GetTargetPosZ()
{
return targetLayerZPosition;
}


G4double DetectorConstruction::GetTargetThickness()
{
return targetLayerXDimension;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
