
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
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fPBox(0), fLBox(0), fMaterial(0), fMagField(0)
{
    fDetThickness = 25*um;
    DefineMaterials();
//    SetDetMaterial("G4_Galactic");
    fDetectorMessenger = new DetectorMessenger(this);
    detPath = "/home/groot/work/geom_import/SIDAR/";
    fuffa = 2;
    targetPhi = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
// define Elements

    G4double z,a;
//    G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
//    G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);

    G4int zz;
    G4int nElem;
    G4NistManager* pMatMan = G4NistManager::Instance();
    G4Element* elCr = pMatMan->FindOrBuildElement(zz=24);
    G4Element* elFe = pMatMan->FindOrBuildElement(zz=26);
    G4Element* elNi = pMatMan->FindOrBuildElement(zz=28);
    G4Element* elC =  pMatMan->FindOrBuildElement(zz=6);
    G4Element* elZn = pMatMan->FindOrBuildElement(30);
    G4Element* elCu = pMatMan->FindOrBuildElement(29);
    G4Element* elLi = pMatMan->FindOrBuildElement(zz=3);
    G4Element* elH = pMatMan->FindOrBuildElement(zz=1);
    G4Element* elF = pMatMan->FindOrBuildElement(zz=9);
    G4Element* elCo = pMatMan->FindOrBuildElement(zz=27);
    G4Element* elW  = pMatMan->FindOrBuildElement(zz=74);
    G4Element* elTi = pMatMan->FindOrBuildElement(zz=22);
    G4Element* elN  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
    G4Element* elO  = new G4Element("Oxygen", "O", z= 8., a=16.00*g/mole);
    G4Element* elD  = new G4Element("Deuterium","D",z=1.,2.01*g/mole);
    G4Element* elB = pMatMan->FindOrBuildElement(zz=5);

// define materials

    G4double density;
    G4int ncomponents;
    G4double fractionmass;

    G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
    Air->AddElement(elN, fractionmass= 70.*perCent);
    Air->AddElement(elO, fractionmass= 30.*perCent);

    aluminium = new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);
    G4Material*  silicon = new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
    G4Material *carbonMaterial = new G4Material("Carbonium", z=6., a=12.017*g/mole, density= 2.253*g/cm3);  // ATTENZIONE AL FATTORE 1000 (g->kg)
    G4Material* magnesiumMaterial = new G4Material("Magnesium", z=12., a=24.305*g/mole, density= 1.737*g/cm3);  // NELLA DENSITA'
    G4Material* lithiumMaterial = new G4Material("lithium", z=3., a=6.94*g/mole, density= 0.534*g/cm3);
    
    G4Material* Sn112Material = new G4Material("Sn112", z=50., a=111.904818*g/mole, density= 7.28*g/cm3);
    G4Material *AuMaterial = new G4Material("Gold", z=79., a=196.97*g/mole, density= 19.32*g/cm3);
    G4Material *PbMaterial = new G4Material("Lead", z=82., a=207.19*g/mole, density= 11.35*g/cm3);
    G4Material* Titanium = new G4Material("Titanium", z=22, a=47.867*g/mole, density= 4.506*g/cm3);
    G4Material* Cadmium = new G4Material("Cadmium", z=48, a=112.411*g/mole, density= 8.65*g/cm3);

    G4Material* Kapton = new G4Material("Kapton", density = 1.413*g/cm3, nElem=4); // Kapton Dupont de Nemur (density: 1.396-1.430, get middle )
    Kapton->AddElement(elO,5);
    Kapton->AddElement(elC,22);
    Kapton->AddElement(elN,2);
    Kapton->AddElement(elH,10);

    G4Material* Brass = new G4Material("Brass", density = 8.70*g/cm3, nElem=2);
    Brass->AddElement(elZn, 3);
    Brass->AddElement(elCu, 7);

    G4Material*  stainlesssteel = new G4Material("StainlessSteel", density = 8.0 *g/cm3, nElem=4);
    stainlesssteel->AddElement(elFe, 0.7366);
    stainlesssteel->AddElement(elCr, 0.1820);
    stainlesssteel->AddElement(elNi, 0.0810);
    stainlesssteel->AddElement(elC, 0.0004);

    G4Material* LiH = new G4Material("LiH", density = 0.82*g/cm3, nElem=2);
    LiH->AddElement(elLi, 1);
    LiH->AddElement(elH, 1);

    G4Material* LiF = new G4Material("LiF", density = 2.64*g/cm3, nElem=2);
    LiF->AddElement(elLi, 1);
    LiF->AddElement(elF, 1);

// CH2
// Polypropelene - polyethylene
//   G4Material* CH2 = new G4Material ("CH2" , 0.91*g/cm3, 2);
   density = 0.9*g/cm3; // 0.935*g/cm3; // 0.857 g/cm3
   G4Material* CH2 = new G4Material("CH2",density, nElem=2);
   CH2->AddElement(elC,1);
   CH2->AddElement(elH,2);

   density = 1.06*g/cm3;
   G4Material* CD2 = new G4Material("CD2",density, nElem=2);
   CD2->AddElement(elC,1);
   CD2->AddElement(elD,2);

   G4Material* Paraffin = new G4Material("Paraffin",  density= 0.94*g/cm3, nElem=2);
//    Paraffin->AddElement(elC,52);
//    Paraffin->AddElement(elH,25);
    Paraffin->AddElement(elC,  85.6*perCent);
    Paraffin->AddElement(elH,  14.4*perCent);

// Mylar
    G4Material* Mylar = new G4Material("Mylar", density = 1.39*g/cm3, nElem=3);
    Mylar->AddElement(elO,2);
    Mylar->AddElement(elC,5);
    Mylar->AddElement(elH,4);

// MATERIAL DEFINIOTION VIA THE NIST MATERIAL BUILDER
// window?

    G4Element* elSi = pMatMan->FindOrBuildElement(14);
    G4Material* SiO2 = new G4Material("quartz", density = 2.200*g/cm3, nElem=2);
    SiO2->AddElement(elSi,1);
    SiO2->AddElement(elO,2);

//Havar

    G4Material* Havar = new G4Material("Havar", density= 8.3*g/cm3, nElem=5);
    Havar->AddElement(elCr, 0.1785);
    Havar->AddElement(elFe, 0.1822);
    Havar->AddElement(elCo, 0.4452);
    Havar->AddElement(elNi, 0.1310);
    Havar->AddElement(elW , 0.0631);


   density = 1.2*g/cm3;
   G4Material* Epoxy = new G4Material("Epoxy", density, nElem=2);
   Epoxy->AddElement(elH, 2);
   Epoxy->AddElement(elC, 2);
  
//FR4 (Glass + Epoxy)
   density = 1.86*g/cm3;
   G4Material* FR4 = new G4Material("FR4", density, nElem=2);
   FR4->AddMaterial(SiO2, fractionmass=0.528);
   FR4->AddMaterial(Epoxy, fractionmass=0.472);

// Vacuum
    density = 6.078e-12 *g/cm3; ;//2.376e-15 *g/cm3;
    G4double temperature = 300. *kelvin;
    G4double pressure = 1.e-8 *bar;//; 2.e-7 *bar;
    G4Material* Vacuum = new G4Material("Vacuum", density, nElem=1,kStateGas,temperature,pressure);
    Vacuum->AddMaterial(Air,1.);

    G4Element* elTl = new G4Element("Thallium", "Tl", z=81., a=204.383*g/mole);
    G4Element* elLa = new G4Element("Lanthanum", "La", z=57.,a=138.90547*g/mole);
    G4Element* elBr = new G4Element("Bromium", "Br", z=35., a=79.904*g/mole);
    G4Element* elCe = new G4Element("Cerium", "Tl", z=58., a=140.116*g/mole);
    G4Element* elGe = new G4Element("Germanium", "Ge", z=32, a=72.630*g/mole);

    //LaBr3
    G4Material* LaBr3 = new G4Material("LaBr3", density = 5.07*g/cm3, ncomponents=2);
    LaBr3->AddElement(elLa, nElem=1);
    LaBr3->AddElement(elBr, nElem=3);

    //LaBr3_Ce
    G4Material* LaBr3_Ce = new G4Material("LaBr3_Ce", density = 5.08*g/cm3, ncomponents=2);
    LaBr3_Ce->AddMaterial(LaBr3, 99.5*perCent);
    LaBr3_Ce->AddElement(elCe, 0.5*perCent);

// CeBr3
    G4Material* CeBr3 = new G4Material("CeBr3", density = 5.1*g/cm3, ncomponents=2);
    CeBr3->AddElement(elCe, nElem=1);
    CeBr3->AddElement(elBr, nElem=3);


//He3 def from geant4 forum:
  //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/materials/122.html?inline=-1
  G4int protons=2, neutrons=1, nucleons=protons+neutrons;
  G4double atomicMass = 3.016*g/mole;
  G4Isotope* he3 = new G4Isotope("He3", protons, nucleons, atomicMass);
  G4Element* He3 = new G4Element("Helium3", "He3", 1);
  He3->AddIsotope(he3, 100*perCent);
  pressure = 4*bar;
  temperature = 293*kelvin;
  G4double molar_constant = CLHEP::Avogadro*CLHEP::k_Boltzmann;  //from clhep
  density = (atomicMass*pressure)/(temperature*molar_constant);
// fix it 
//  G4cout<<"3He Tube density: " << density << G4endl; // volume  l:660mm  d:24mm , 298577 mm3 // pressure: 4bar mass: 144.4g; density is 0.48 g/cm3 //ATTENZ: mm3 -->cm3
  density = 0.48*g/cm3; // 177.43*mg/cm3;  // 0.0013881658858234*g/cm3;  (0.1346*mg/cm3) <-- at 5bar? //  
  G4cout<<"3He Tube density: " << density << G4endl;
  G4Material* Helium3 = new G4Material("Helium3", density, 1, kStateGas, temperature, pressure);
  Helium3->AddElement(He3, 100*perCent);

//Diborano
  atomicMass = 27.7*g/mole;
  pressure = 52.5*bar; // 300*0.07*bar;
  temperature = 86.*kelvin;
//  molar_constant = CLHEP::Avogadro*CLHEP::k_Boltzmann;  //from clhep
  density = (atomicMass*pressure)/(temperature*molar_constant);
 density = 0.00118 *g/cm3; 

//  density = 1.e21 // what about the expansion?
  G4Material* Diborane = new G4Material("Diborane", density, 2, kStateGas, temperature, pressure);
//    G4Material* Diborane = new G4Material("Diborane", density = 0.00118 *g/cm3, ncomponents=2);
    Diborane->AddElement(elB, nElem=2);
    Diborane->AddElement(elH, nElem=6);


  G4Material* SiC = new G4Material("SiC", density= 3.22*g/cm3, ncomponents=2);
    SiC->AddElement(elC, nElem=1);
    SiC->AddElement(elSi, nElem=1);


// Eljen EJ301 Density [g/cc] 0.874 #H atoms/cc 4.82 x 10 22 Ratio H:C 1.212

  G4Material* EJ301 = new G4Material("EJ301", density= 0.874*g/cm3, 2, kStateLiquid);
	G4double nC = 4.82/1.212; // No. of C Atoms per cm3 (x10^22)
	G4double nH = 4.82; // No. of H Atoms per cm3 (x10^22)
	EJ301->AddElement(elH, nH/(nH+nC) );
	EJ301->AddElement(elC, nC/(nH+nC) );

  G4Element* Ba = new G4Element("Barium", "Ba", z=56, a= 137.327*g/mole);

  G4Material* BaF2 = new G4Material("BaF2", density= 4.89*g/cm3, 2, kStateSolid);
	BaF2->AddElement(Ba, nElem=1);
	BaF2->AddElement(elF, nElem=2);

// Te = a1 [1- exp(-a2 T_p ^a3)] + a4 Tp  

    PTMaterial = Brass;
    targetLayerMaterial = Sn112Material;//; LiF; // lithiumMaterial; // LiH; // lithiumMaterial;
    backingLayerMaterial =  Mylar; // AuMaterial; // carbonMaterial;
    worldMaterial = Vacuum; // POLYCUBE // Vacuum; //Air
    AirMaterial = Air;
    ChamberMaterial = aluminium; // stainlesssteel;
    windowMaterial = aluminium; // tainlesssteel; // Kapton; //aluminium; // Titanium;// SiO2; // Havar; //
    DetectorMaterial = Vacuum;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    targetLayerXDimension = 1.*nm; // 0.0018726592 *mm;
    fDetThickness = targetLayerXDimension/2;

    G4double backingLayerXDimension = 0.0013*mm; //200.*nm;
    fBackingThickness = backingLayerXDimension/2.;
    backingPos = "before";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
// Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

// G4Colour  cyan (0.0, 1.0, 1.0);
  G4VisAttributes *theWhite = new G4VisAttributes( G4Colour(255/255., 255/255., 255/255. ));
    theWhite -> SetVisibility(true);
    theWhite -> SetForceSolid(true);

    G4VisAttributes *darkGreen = new G4VisAttributes( G4Colour(0/255., 100/255., 0/255. ));
    darkGreen -> SetVisibility(true);
    darkGreen -> SetForceSolid(true);

    G4VisAttributes *cyan = new G4VisAttributes(G4Colour(0/255., 255/255., 255/255.));
    cyan -> SetVisibility(true);
    cyan -> SetForceWireframe(true);
    cyan->SetForceAuxEdgeVisible(true);

    G4VisAttributes *yellow = new G4VisAttributes(G4Colour(255/255., 255/255. ,0/255.));
    yellow -> SetVisibility(true);
    yellow -> SetForceWireframe(true);
    yellow->SetForceAuxEdgeVisible(true);
    yellow-> SetForceLineSegmentsPerCircle(64); //  16 segments of 4 lines each

    G4VisAttributes *red = new G4VisAttributes(G4Colour(255/255., 0/255., 0/255.));
    red -> SetVisibility(true);
    red -> SetForceWireframe(true);
    red->SetForceAuxEdgeVisible(true);
    red-> SetForceLineSegmentsPerCircle(64); //  16 segments of 4 lines each

    G4double worldVolumeXDimension = 5. *m;
    G4double worldVolumeYDimension = 1. *m;
    G4double worldVolumeZDimension = 1. *m;

    G4double gapBeforeTarget = 0.0001 *mm;
    G4double gapAfterTarget = 0.0001 *mm;// 0.01 *mm;

    targetLayerXPosition = 0. *mm;
    targetLayerYPosition = 0. *mm;
    targetLayerZPosition = /* gapBeforeTarget + */ 0. *mm ;

    G4double backingLayerXPosition = 0. *mm;
    G4double backingLayerYPosition = 0. *mm;

    G4double backingLayerZPosition;

    if(backingPos == "before")
        backingLayerZPosition = gapBeforeTarget +  fBackingThickness/2. + fDetThickness/2.;
    else
        backingLayerZPosition = -gapAfterTarget -  fBackingThickness/2. - fDetThickness/2.;
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

   fLBox->SetVisAttributes(G4VisAttributes::GetInvisible());

   G4RotationMatrix* rotm  = new G4RotationMatrix();
   //rotm->rotateY(targetPhi*deg);

// TARGET: 27Mg Layer
    G4Tubs* solidtargetLayer = new G4Tubs("targetLayer",0.98*cm,0.99*cm,fDetThickness/2.,0.*deg,360.*deg);
    logictargetLayer = new G4LogicalVolume(solidtargetLayer,targetLayerMaterial,"targetLayerLogic");
    G4VPhysicalVolume* physitargetLayer = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(targetLayerXPosition,targetLayerYPosition,targetLayerZPosition)),logictargetLayer,"targetLayerPhys",fLBox,false,0);
    logictargetLayer -> SetVisAttributes(theWhite);
/*
    G4Tubs* solidbackingLayer = new G4Tubs("backingLayer",0.,0.99*cm,fBackingThickness/2.,0.*deg,360.*deg);
    G4LogicalVolume* logicbackingLayer = new G4LogicalVolume(solidbackingLayer,backingLayerMaterial,"backingLayerLogic");
    G4VPhysicalVolume* physibackingLayer = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,backingLayerZPosition)),logicbackingLayer,"backingLayerPhys",fLBox,false,0);
    logicbackingLayer -> SetVisAttributes(red);
    rotm->rotateY(-targetPhi*deg);

    G4VisAttributes * backVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.,1.));
   backVisAtt -> SetForceSolid(true);
    logicbackingLayer -> SetVisAttributes(backVisAtt);*/
/*
    G4Tubs* window = new G4Tubs("window",0.*cm,2.*cm,1.*cm,0.*deg,360.*deg);
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


    G4Tubs* apizza = new G4Tubs("pizza",1.*cm,13.*cm,8.*um,0.*deg,360.*deg);
    G4LogicalVolume* logicpizza = new G4LogicalVolume(apizza, windowMaterial,"logicpizza");
    G4VPhysicalVolume* physipizza = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,10.*cm)),logicpizza,"physipizza",fLBox,false,0);
    logicpizza->SetVisAttributes(theWhite); // SetVisAttributes(G4Color(1., 1., 1., 0.6));

    G4Tubs* FCup = new G4Tubs("FCup",2.*cm,3.*cm,26.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicFCup = new G4LogicalVolume(FCup, PTMaterial,"logicFCup"); // BRASS OR COPPER?
    G4VPhysicalVolume* physiFCup = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,(99+26)*cm)),logicFCup,"physiFCup",fLBox,false,0);
    logicFCup ->  SetVisAttributes(theWhite);
    G4Tubs* FCup2 = new G4Tubs("FCup2",0.*cm,2.*cm,2.*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicFCup2 = new G4LogicalVolume(FCup2, PTMaterial,"logicFCup2"); // BRASS OR COPPER?
    G4VPhysicalVolume* physiFCup2 = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(backingLayerXPosition,backingLayerYPosition,(101.+52.)*cm)),logicFCup2,"physiFCup2",fLBox,false,0);
    logicFCup2 ->  SetVisAttributes(red);
*/
// CAD model rotation.
/*    G4RotationMatrix * rot = new G4RotationMatrix();
    rot->rotateY(0.*deg); // 90 or -90?

    G4ThreeVector Theoffset;
    // Load CAD file as tessellated solid //

// IN STP FILE: distance from reference frame origin to the center of the bottom target holder circle ( =====|O|O|O|O|0| ) is (524.5*mm, 5.*mm, 0.2*mm), with target holder along Z-axis and X the beam direction.

    G4double X_offset = 0.*cm; // 52.45*cm;       OLD -> //  -15.363*cm; // DET mid  -15.103*cm; SENSITIVE DET mid  -14.363*cm;
    G4double Y_offset = 0.*cm; // -0.5*cm;
    G4double Z_offset =  0.*cm;
    Theoffset.set(Z_offset, Y_offset, X_offset); // FIX THIS

// Note that offset is applied to the points in mesh directly before placement.

// offset on ELISSAPI files --> 2.4412 is distance between the centroid of the small circle and the origin
    //G4ThreeVector Theoffset3(0.,0.,+(20.8+2.4412)*cm);   not needed anymore
*/

/*
    G4Tubs* LeadFloor = new G4Tubs("LeadFloor",0.,1.5*m,2.5*cm,0.*deg,360.*deg);
    G4LogicalVolume* logicLeadFloor = new G4LogicalVolume(LeadFloor,DetectorMaterial,"LeadFloorLogic");
    G4VPhysicalVolume* physiLeadFloor = new G4PVPlacement(G4Transform3D(*rotm,G4ThreeVector(0.,0.,-10.5*cm)),logicLeadFloor,"LeadFloorPhys",fLBox,false,0);
    logicLeadFloor -> SetVisAttributes(G4Color(0.3, 0.3, 0.3, 1.));
*/
//    rotm->rotateY(-targetPhi*deg);



    std::ifstream infile;
    infile.open ("inputFileDetectors.txt");
    G4String name, material;
    G4int i = 0;

    while(infile >> name >> material)
    {

        G4double z, y, x, rotx, roty, rotz;

        infile>>x>>y>>z>>rotx>>roty>>rotz;

        G4ThreeVector position(x*cm,y*cm,z*cm);

        G4RotationMatrix * rot3 = new G4RotationMatrix(); // [12]; // = new G4RotationMatrix(); // rot; //
        if(roty!=0) rot3->rotateY(roty*deg); // 90 or -90?
        if(rotx!=0) rot3->rotateX(rotx*deg);
        if(rotz!=0) rot3->rotateZ(rotz*deg);

        G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
        if (pttoMaterial)
        {
             DetectorMaterial = pttoMaterial;
        }
        else
        {
            G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
            << material << " not found" << G4endl;
            exit(0);
        }

        G4String stripfilename = detPath+name;
        //out<<stripfilename<<' ';
        //stripfilename.append(std::to_string(i+1));
        //out<<stripfilename<<' ';
        //stripfilename.append(".stl");
        //out<<stripfilename<<' ';

        G4cout<<stripfilename.c_str()<<G4endl;

// REMOVE STL?
        G4String stripID = name;  /* "YY1-";  
        stripID.append(std::to_string(i));*/
        G4String stripIDlog = stripID;
        stripIDlog.append("_logical");

        G4String stripIDphys = name; // "YY1-";
/*        G4String TheDetID = std::to_string(i);
        if(TheDetID.size()==1)
        {
            stripIDphys.append("0");
        }
        stripIDphys.append(TheDetID);*/
        stripIDphys.append("_physical");
/*
auto bunny_mesh = CADMesh::TessellatedMesh::FromSTL("./bunny.stl");

        auto bunny_logical = new G4LogicalVolume( bunny_mesh->GetSolid() 
                                                 , water
                                                 , "logical"
                                                 , 0, 0, 0
*/
      if(name.substr(name.find_last_of(".") +1) == "stl" || name.substr(name.find_last_of(".") +1) == "STL") { // fn.substr(fn.find_last_of(".") + 1) == "conf")
          auto the_mesh = CADMesh::TessellatedMesh::FromSTL( (char*) stripfilename.c_str());
          the_mesh->SetScale(mm);
          the_mesh->SetReverse(false);
          QQ3_logical[i] = new G4LogicalVolume(the_mesh->GetSolid(), DetectorMaterial, stripIDlog, 0, 0, 0);
      }
      else {
          auto the_mesh = CADMesh::TessellatedMesh::FromPLY( (char*) stripfilename.c_str());
          the_mesh->SetScale(mm);
          the_mesh->SetReverse(false);
          QQ3_logical[i] = new G4LogicalVolume(the_mesh->GetSolid(), DetectorMaterial, stripIDlog, 0, 0, 0);

      }

    //    QQ3_physical[i] = new G4PVPlacement(rot, Theoffset, QQ3_logical[i],
          QQ3_physical[i] = new G4PVPlacement (G4Transform3D(*rot3/*[i]*/,position), QQ3_logical[i], stripIDphys, fLBox, false, 0);
          QQ3_logical[i]->SetVisAttributes(G4Color(0.5, 0.5, 0.5, 1.));

        i++;
 }

    // Note that offset is applied to the points in mesh directly before placement.


    G4cout<<"Detector ready."<<G4endl;

//
    PrintParameters();
    return fPBox; // return the root volume
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    G4cout << "\n Maybe it's a lie: The Backing is " << G4BestUnit(fBackingThickness,"Length")
    << " of " << backingLayerMaterial->GetName() <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial)
    {
         targetLayerMaterial = pttoMaterial;
    }
    else
    {
        G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
        << materialChoice << " not found" << G4endl;
    }

}\

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBackingMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial)
    {
         backingLayerMaterial = pttoMaterial;
    }
    else
    {
        G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
        << materialChoice << " not found" << G4endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBackingThickness(G4double value)
{
    fBackingThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetThickness(G4double value)
{
    fDetThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPath(G4String value)
{
    detPath = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
    //apply a global uniform magnetic field along Z axis
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
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

void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetPosZ()
{
    return targetLayerZPosition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetThickness()
{
    return fDetThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetPhi(G4double phi)
{
    targetPhi = phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBackingPos(G4String pos)
{
    backingPos = pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetReactionType(G4int type)
{
    fuffa = type;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
