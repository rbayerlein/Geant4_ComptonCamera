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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "RunControl.hh"
#include "SGDetector.hh"
#include "PMTDetector.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <CLHEP/Units/PhysicalConstants.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ 
  useBGO = false; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

	abs_en.clear();
	ifstream in2;
	in2.open("abs_energy.txt");
	string dummy3;
	if(!in2) cout << "abs_energy: file does not exist" << endl;
	else{
		while(in2 >> dummy3){
			abs_en.push_back(dummy3);
		}
	}	
  in2.close();

	abs_le.clear();
	ifstream in3;
	in3.open("abs_length.txt");
	string dummy4;
	if(!in3) cout << "abs_length: file does not exist" << endl;
	else{
		while(in3 >> dummy4){
			abs_le.push_back(dummy4);
		}
	}	
  in3.close();

RunControl * RC= RunControl::GetInstance();
stringstream ss_volumeList;
ss_volumeList << RC->RAWfileName << "_VolumeList.txt";
string s_vL = ss_volumeList.str();
ofstream o_volumes;
o_volumes.open(s_vL.c_str());


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CHERENKOV GEOMETRICAL PARAMETERS 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //SiPM Parameters

 int n = RC->arraySideLength_x; 				//number of channels x
 int v = RC->arraySideLength_y;					//number of channels y
 int p = RC->number_of_arrays_x;
 int w = RC->number_of_arrays_y;
 int z1 = RC->number_of_arrays_z;							//number of arrays sides
 int z2 = RC->arraySideLength_z;							//number of channels sides

 G4double space = RC->space; 					//space between individual SiPM channels
 G4double spacem = RC->spacem;        //space between SiPM arrays

	//SiPM_ref
  G4double ref_thick = RC->ref_thick;

   G4double SiPMchannel_size_xy = RC->SiPMchannel_size_xy;  //xy-size SiPM
   G4double SiPM_thick = RC->SiPM_thick;
   G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_Galactic");

 //PMMA Parameters

 G4double space_x = RC->space_x;        //overlap PMMA above SiPM x
 G4double space_y = RC->space_y;        //overlap PMMA above SiPM y
  PMMA_size_x = p*(n*SiPMchannel_size_xy+space*(n-1))+(p-1)*spacem+2*space_x; // declared in header file
  PMMA_size_y = w*(v*SiPMchannel_size_xy+space*(v-1))+(w-1)*spacem+2*space_y; // declared in header file
 G4double PMMA_thick = RC->PMMA_thick;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MATERIALS 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4Material* PMMA_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

	 int natoms=0;
	 int z=0;
	 int ncomp=0;
	 G4String symbol;
 	 G4double a = 1.01*g/mole; G4Element* elH = new G4Element("Hydrogen",symbol="H",z=1.,a); 
 	 a = 16.00*g/mole; 
 	 G4Element* elO  =  new G4Element("Oxygen",symbol="O",z=8.,a); 
	 G4double density = 1.000*g/cm3; 
	 G4Material* ref_mat = new G4Material("Water",density, ncomp=2); 
 	 ref_mat->AddElement(elH, natoms=2); 
	 ref_mat->AddElement(elO, natoms=1);



 //Envelope parameters
  
  G4double env_sizeXY = 10*PMMA_size_x;
  G4double env_sizeZ = 10*PMMA_size_x;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

 //World Parameters

  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  //G4Material* SG_mat = nist->FindOrBuildMaterial("G4_PLEXIGLASS");

  //Bismuth

  G4double M_Bi = 209*g/mole;
  G4int Z_Bi = 83;
  G4String symbol_Bi = "Bi";
  G4Element* elBi = new G4Element("Bismuth", symbol_Bi, Z_Bi, M_Bi);

  //Germanium
  G4double M_Ge = 74*g/mole;
  G4int Z_Ge = 32;
  G4String symbol_Ge = "Ge";
  G4Element* elGe = new G4Element("Germanium", symbol_Ge, Z_Ge, M_Ge);

  //Lanthanum
  G4double M_La = 138.905*g/mole;
  G4int Z_La = 57;
  G4String symbol_La = "La";
  G4Element* elLa = new G4Element("Lanthanum", symbol_La, Z_La, M_La);

  //Bromine
  G4double M_Br = 79.904*g/mole;
  G4int Z_Br = 35;
  G4String symbol_Br = "Br";
  G4Element* elBr = new G4Element("Bromine", symbol_Br, Z_Br, M_Br);

  // BGO 
  int nComp_BGO = 3;
  double density_BGO = 7.12*g/cm3;

  // LaBr3
  int nComp_LaBr3 = 2;
  double density_LaBr3 = 5.08*g/cm3;

  G4Material* SG_mat;
  if(useBGO){
    SG_mat = new G4Material("BGO", density_BGO, nComp_BGO);
    SG_mat->AddElement(elBi, 4);
    SG_mat->AddElement(elGe, 3);
    SG_mat->AddElement(elO, 12);
  }else{
    SG_mat = new G4Material("LaBr3", density_LaBr3, nComp_LaBr3);
    SG_mat->AddElement(elLa, 1);
    SG_mat->AddElement(elBr, 3);
  }

  // Option to switch on/off checking of volume overlaps
  G4bool checkOverlaps = true;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BUILD VOLUMES AND GEOMETRIES
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //     
  // World
  //
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name

  G4VPhysicalVolume* physEnv =              
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

G4SDManager* SDManager = G4SDManager::GetSDMpointer();

SGDetector* sensitiveGammaDetector;
sensitiveGammaDetector = new SGDetector("scintillator");

PMTDetector* sensitivePMTDetector;
sensitivePMTDetector = new PMTDetector("PMTDetector");

SDManager->AddNewDetector(sensitiveGammaDetector);
SDManager->AddNewDetector(sensitivePMTDetector);

//Scattered Gamma Detector Matrices 

G4Box *solidGammaDetector_Backside = new G4Box("solidGammaDetector_Backside",RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_z/2);
G4Box *solidGammaDetector_Sides = new G4Box("solidGammaDetector_Sides",RC->SG_ChannelSize_z/2, RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_xy/2);
G4Box *solidGammaDetector_TopBottom = new G4Box("solidGammaDetector_TopBottom",RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_z/2, RC->SG_ChannelSize_xy/2);

G4Box *solidPMTs_Backside = new G4Box("solidPMTs_Backside", RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_xy/2, RC->PMTs_ChannelSize_z/2);
G4Box *solidPMTs_Sides = new G4Box("solidPMTs_Sides", RC->PMTs_ChannelSize_z/2, RC->SG_ChannelSize_xy/2, RC->SG_ChannelSize_xy/2);
G4Box *solidPMTs_TopBottom = new G4Box("solidPMTs_TopBottom", RC->SG_ChannelSize_xy/2, RC->PMTs_ChannelSize_z/2, RC->SG_ChannelSize_xy/2);

G4VisAttributes * SG_visAtt = new G4VisAttributes(G4Colour(0.1, 0.5, 1, 0.4));
G4VisAttributes * PMTs_visAtt = new G4VisAttributes(G4Colour(0.2, 1, 0.6, 0.3));

G4OpticalSurface* SG_OptSurface = new G4OpticalSurface("SG_OptSurface");  // for internal reflection
G4OpticalSurface* SG_OptWrapper = new G4OpticalSurface("SG_OptWrapper");  // for photons coming from outside. 
G4OpticalSurface* PMT_InnerSurface = new G4OpticalSurface("PMT_InnerSurface");
G4OpticalSurface* PMT_TowardsScintSurface = new G4OpticalSurface("PMT_TowardsScintSurface");
G4OpticalSurface* PMT_TowardsPMTSurface = new G4OpticalSurface("PMT_TowardsPMTSurface");


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PLACEMENT OF SCINTILLATORS AND PMTs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int SG_SideLength_x= RC-> SG_SideLength_x;
int SG_SideLength_y= RC-> SG_SideLength_y;
int SG_SideLength_z= RC-> SG_SideLength_z;

// BACKSIDE

G4LogicalVolume* SG_logic_Backside[100][100];
G4VPhysicalVolume* SG_physical_Backside[100][100];
G4LogicalBorderSurface* SG_Surface_Backside[100][100];  // for internal reflection
G4LogicalBorderSurface* SG_Wrapper_Backside[100][100];  // for photons coming from outside. 

G4LogicalVolume* PMTs_logic_Backside[100][100];
G4VPhysicalVolume* PMTs_physical_Backside[100][100];
G4LogicalBorderSurface* PMT_Inner_Surface_Backside[100][100];
G4LogicalBorderSurface* PMT_TowardsScintSurface_Backside[100][100];
G4LogicalBorderSurface* PMT_TowardsPMTSurface_Backside[100][100];

for (int i = 0; i < SG_SideLength_x; ++i)
{
  for (int j = 0; j < SG_SideLength_y; ++j)
  { 
  char name[512];
  sprintf(name, "SG_Back_%d_%d", i, j);
  SG_logic_Backside[i][j] = new G4LogicalVolume(solidGammaDetector_Backside, SG_mat, name);
  SG_logic_Backside[i][j]->SetVisAttributes(SG_visAtt);
  SG_logic_Backside[i][j]->SetSensitiveDetector(sensitiveGammaDetector);

  char name_PMTs[512];
  sprintf(name_PMTs, "PMTs_Back_%d_%d", i, j);  
  PMTs_logic_Backside[i][j] = new G4LogicalVolume(solidPMTs_Backside, SG_mat, name_PMTs);
  PMTs_logic_Backside[i][j]->SetVisAttributes(PMTs_visAtt);
  PMTs_logic_Backside[i][j]->SetSensitiveDetector(sensitivePMTDetector);
  // here set sensitive detector for PMTs

  G4double pos_x = (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_x/2 - (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 - i* (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_y = (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_y/2 - (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 - j* (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_z = RC->SG_DetectorPosition_z+0.5*RC->SG_ChannelSize_z;
  G4ThreeVector position = G4ThreeVector(pos_x, pos_y, pos_z);
  o_volumes << name_PMTs << "\t" << pos_x << "\t" << pos_y << "\t" << pos_z << endl;

  SG_physical_Backside[i][j] = new G4PVPlacement(
    0,                          // rotation 
    position,                   // position
    SG_logic_Backside[i][j],    // logical volume
    name,                       // name
    logicEnv,                   // logical volume of mother
    false,                      // boolean operation
    0,                          // copy number
    checkOverlaps);             // check for overlaps
  
  G4double pos_z_PMT = RC->SG_DetectorPosition_z+RC->SG_ChannelSize_z + 0.5*RC->PMTs_ChannelSize_z;
  G4ThreeVector position_PMT = G4ThreeVector(pos_x, pos_y, pos_z_PMT);


  PMTs_physical_Backside[i][j] = new G4PVPlacement(
    0,
    position_PMT,
    PMTs_logic_Backside[i][j],
    name_PMTs,
    logicEnv,
    false,
    0,
    checkOverlaps);

  char surf_name[512];
  sprintf(surf_name, "SG_Surf_Back_%d_%d", i, j);
  SG_Surface_Backside[i][j] = new G4LogicalBorderSurface(surf_name,SG_physical_Backside[i][j],physEnv,SG_OptSurface);

  char wrapper_name[512];
  sprintf(wrapper_name, "SG_Wrapper_Back_%d_%d", i, j);
  SG_Wrapper_Backside[i][j] = new G4LogicalBorderSurface(wrapper_name,physEnv, SG_physical_Backside[i][j], SG_OptWrapper);

  char PMT_surf_name[512];
  sprintf(PMT_surf_name, "PMT_InnerSurf_Back_%d_%d", i ,j);
  PMT_Inner_Surface_Backside[i][j] = new G4LogicalBorderSurface(PMT_surf_name, PMTs_physical_Backside[i][j], physEnv, PMT_InnerSurface);

  char PMT_toScint_surf_name[512];
  sprintf(PMT_toScint_surf_name, "PMT_ToScint_Surf_Back_%d_%d", i ,j);
  PMT_TowardsScintSurface_Backside[i][j] = new G4LogicalBorderSurface(PMT_toScint_surf_name, PMTs_physical_Backside[i][j], SG_physical_Backside[i][j], PMT_TowardsScintSurface);

  char PMT_toPMT_surf_name[512];
  sprintf(PMT_toPMT_surf_name, "PMT_ToPMT_Surf_Back_%d_%d", i ,j);
  PMT_TowardsPMTSurface_Backside[i][j] = new G4LogicalBorderSurface(PMT_toPMT_surf_name, SG_physical_Backside[i][j], PMTs_physical_Backside[i][j], PMT_TowardsPMTSurface);

  // here set surfaces for PMTs
  }
}

// LEFT SIDE

G4LogicalVolume* SG_logic_Left[100][100];
G4VPhysicalVolume* SG_physical_Left[100][100];
G4LogicalBorderSurface* SG_Surface_Left[100][100];
G4LogicalBorderSurface* SG_Wrapper_Left[100][100];  // for photons coming from outside. 

G4LogicalVolume* PMTs_logic_Left[100][100];
G4VPhysicalVolume* PMTs_physical_Left[100][100];
G4LogicalBorderSurface* PMT_Inner_Surface_Left[100][100];
G4LogicalBorderSurface* PMT_TowardsScintSurface_Left[100][100];
G4LogicalBorderSurface* PMT_TowardsPMTSurface_Left[100][100];

for (int i = 0; i < SG_SideLength_z; ++i)
{
  for (int j = 0; j < SG_SideLength_y; ++j)
  {
  char name[512];
  sprintf(name, "SG_Left_%d_%d", i, j);
  SG_logic_Left[i][j] = new G4LogicalVolume(solidGammaDetector_Sides, SG_mat, name);
  SG_logic_Left[i][j]->SetVisAttributes(SG_visAtt);
  SG_logic_Left[i][j]->SetSensitiveDetector(sensitiveGammaDetector);

  char name_PMTs[512];
  sprintf(name_PMTs, "PMTs_Left_%d_%d", i, j);  
  PMTs_logic_Left[i][j] = new G4LogicalVolume(solidPMTs_Sides, SG_mat, name_PMTs);
  PMTs_logic_Left[i][j]->SetVisAttributes(PMTs_visAtt);
  PMTs_logic_Left[i][j]->SetSensitiveDetector(sensitivePMTDetector);

  G4double pos_x = (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_x/2 + 0.5* RC->SG_ChannelSize_z;
  G4double pos_y = - (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_y/2 + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + j*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_z = RC->SG_DetectorPosition_z-SG_SideLength_z*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + i*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);

  G4ThreeVector position = G4ThreeVector(pos_x, pos_y, pos_z);
  o_volumes << name_PMTs << "\t" << pos_x << "\t" <<  pos_y << "\t" << pos_z << endl;

  SG_physical_Left[i][j] = new G4PVPlacement(
  0,                          // rotation 
  position,                   // position
  SG_logic_Left[i][j],    // logical volume
  name,                       // name
  logicEnv,                   // logical volume of mother
  false,                      // boolean operation
  0,                          // copy number
  checkOverlaps);             // check for overlaps

  G4double pos_x_PMT = pos_x + 0.5*RC->SG_ChannelSize_z + 0.5*RC->PMTs_ChannelSize_z;
  G4ThreeVector position_PMT = G4ThreeVector(pos_x_PMT,pos_y,pos_z);

  PMTs_physical_Left[i][j] = new G4PVPlacement(
    0,
    position_PMT,
    PMTs_logic_Left[i][j],
    name_PMTs,
    logicEnv,
    false,
    0,
    checkOverlaps);

  char surf_name[512];
  sprintf(surf_name, "SG_Surf_Left_%d_%d", i, j);
  SG_Surface_Left[i][j] = new G4LogicalBorderSurface(surf_name,SG_physical_Left[i][j],physEnv,SG_OptSurface);

  char wrapper_name[512];
  sprintf(wrapper_name, "SG_Wrapper_Left_%d_%d", i, j);
  SG_Wrapper_Left[i][j] = new G4LogicalBorderSurface(wrapper_name,physEnv, SG_physical_Left[i][j], SG_OptWrapper);

  char PMT_surf_name[512];
  sprintf(PMT_surf_name, "PMT_InnerSurf_Left_%d_%d", i ,j);
  PMT_Inner_Surface_Left[i][j] = new G4LogicalBorderSurface(PMT_surf_name, PMTs_physical_Left[i][j], physEnv, PMT_InnerSurface);

  char PMT_toScint_surf_name[512];
  sprintf(PMT_toScint_surf_name, "PMT_ToScint_Surf_Left_%d_%d", i ,j);
  PMT_TowardsScintSurface_Left[i][j] = new G4LogicalBorderSurface(PMT_toScint_surf_name, PMTs_physical_Left[i][j], SG_physical_Left[i][j], PMT_TowardsScintSurface);

  char PMT_toPMT_surf_name[512];
  sprintf(PMT_toPMT_surf_name, "PMT_ToPMT_Surf_Left_%d_%d", i ,j);
  PMT_TowardsPMTSurface_Left[i][j] = new G4LogicalBorderSurface(PMT_toPMT_surf_name, SG_physical_Left[i][j], PMTs_physical_Left[i][j], PMT_TowardsPMTSurface);

  }
}


// RIGHT SIDE

G4LogicalVolume* SG_logic_Right[100][100];
G4VPhysicalVolume* SG_physical_Right[100][100];
G4LogicalBorderSurface* SG_Surface_Right[100][100];
G4LogicalBorderSurface* SG_Wrapper_Right[100][100];  // for photons coming from outside. 

G4LogicalVolume* PMTs_logic_Right[100][100];
G4VPhysicalVolume* PMTs_physical_Right[100][100];
G4LogicalBorderSurface* PMT_Inner_Surface_Right[100][100];
G4LogicalBorderSurface* PMT_TowardsScintSurface_Right[100][100];
G4LogicalBorderSurface* PMT_TowardsPMTSurface_Right[100][100];

for (int i = 0; i < SG_SideLength_z; ++i)
{
  for (int j = 0; j < SG_SideLength_y; ++j)
  {
  char name[512];
  sprintf(name, "SG_Right_%d_%d", i, j);
  SG_logic_Right[i][j] = new G4LogicalVolume(solidGammaDetector_Sides, SG_mat, name);
  SG_logic_Right[i][j]->SetVisAttributes(SG_visAtt);
  SG_logic_Right[i][j]->SetSensitiveDetector(sensitiveGammaDetector);

  char name_PMTs[512];
  sprintf(name_PMTs, "PMTs_Right_%d_%d", i, j);  
  PMTs_logic_Right[i][j] = new G4LogicalVolume(solidPMTs_Sides, SG_mat, name_PMTs);
  PMTs_logic_Right[i][j]->SetVisAttributes(PMTs_visAtt);
  PMTs_logic_Right[i][j]->SetSensitiveDetector(sensitivePMTDetector);

  G4double pos_x = -(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_x/2 - 0.5* RC->SG_ChannelSize_z;
  G4double pos_y = - (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_y/2 + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + j*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_z = RC->SG_DetectorPosition_z-SG_SideLength_z*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + i*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);

  G4ThreeVector position = G4ThreeVector(pos_x, pos_y, pos_z);
  o_volumes << name_PMTs << "\t" << pos_x << "\t" <<  pos_y << "\t" << pos_z << endl;

  SG_physical_Right[i][j] = new G4PVPlacement(
  0,                          // rotation 
  position,                   // position
  SG_logic_Right[i][j],    // logical volume
  name,                       // name
  logicEnv,                   // logical volume of mother
  false,                      // boolean operation
  0,                          // copy number
  checkOverlaps);             // check for overlaps

  G4double pos_x_PMT = pos_x - 0.5*RC->SG_ChannelSize_z - 0.5*RC->PMTs_ChannelSize_z;
  G4ThreeVector position_PMT = G4ThreeVector(pos_x_PMT,pos_y,pos_z);

  PMTs_physical_Right[i][j] = new G4PVPlacement(
    0,
    position_PMT,
    PMTs_logic_Right[i][j],
    name_PMTs,
    logicEnv,
    false,
    0,
    checkOverlaps);

  char surf_name[512];
  sprintf(surf_name, "SG_Surf_Right_%d_%d", i, j);
  SG_Surface_Right[i][j] = new G4LogicalBorderSurface(surf_name,SG_physical_Right[i][j],physEnv,SG_OptSurface);

  char wrapper_name[512];
  sprintf(wrapper_name, "SG_Wrapper_Right_%d_%d", i, j);
  SG_Wrapper_Right[i][j] = new G4LogicalBorderSurface(wrapper_name,physEnv, SG_physical_Right[i][j], SG_OptWrapper);

  char PMT_surf_name[512];
  sprintf(PMT_surf_name, "PMT_InnerSurf_Right_%d_%d", i ,j);
  PMT_Inner_Surface_Right[i][j] = new G4LogicalBorderSurface(PMT_surf_name, PMTs_physical_Right[i][j], physEnv, PMT_InnerSurface);

  char PMT_toScint_surf_name[512];
  sprintf(PMT_toScint_surf_name, "PMT_ToScint_Surf_Right_%d_%d", i ,j);
  PMT_TowardsScintSurface_Right[i][j] = new G4LogicalBorderSurface(PMT_toScint_surf_name, PMTs_physical_Right[i][j], SG_physical_Right[i][j], PMT_TowardsScintSurface);

  char PMT_toPMT_surf_name[512];
  sprintf(PMT_toPMT_surf_name, "PMT_ToPMT_Surf_Right_%d_%d", i ,j);
  PMT_TowardsPMTSurface_Right[i][j] = new G4LogicalBorderSurface(PMT_toPMT_surf_name, SG_physical_Right[i][j], PMTs_physical_Right[i][j], PMT_TowardsPMTSurface);

  }
}

// BOTTOM SIDE

G4LogicalVolume* SG_logic_Bottom[100][100];
G4VPhysicalVolume* SG_physical_Bottom[100][100];
G4LogicalBorderSurface* SG_Surface_Bottom[100][100];
G4LogicalBorderSurface* SG_Wrapper_Bottom[100][100];  // for photons coming from outside. 

G4LogicalVolume* PMTs_logic_Bottom[100][100];
G4VPhysicalVolume* PMTs_physical_Bottom[100][100];
G4LogicalBorderSurface* PMT_Inner_Surface_Bottom[100][100];
G4LogicalBorderSurface* PMT_TowardsScintSurface_Bottom[100][100];
G4LogicalBorderSurface* PMT_TowardsPMTSurface_Bottom[100][100];

for (int i = 0; i < SG_SideLength_x; ++i)
{
  for (int j = 0; j < SG_SideLength_z; ++j)
  {
  char name[512];
  sprintf(name, "SG_Bottom_%d_%d", i, j);
  SG_logic_Bottom[i][j] = new G4LogicalVolume(solidGammaDetector_TopBottom, SG_mat, name);
  SG_logic_Bottom[i][j]->SetVisAttributes(SG_visAtt);  
  SG_logic_Bottom[i][j]->SetSensitiveDetector(sensitiveGammaDetector);

  char name_PMTs[512];
  sprintf(name_PMTs, "PMTs_Bottom_%d_%d", i, j);  
  PMTs_logic_Bottom[i][j] = new G4LogicalVolume(solidPMTs_TopBottom, SG_mat, name_PMTs);
  PMTs_logic_Bottom[i][j]->SetVisAttributes(PMTs_visAtt);
  PMTs_logic_Bottom[i][j]->SetSensitiveDetector(sensitivePMTDetector);

  G4double pos_x =  (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_x/2 -0.5*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) - i*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_y = -(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_y/2 -0.5*RC->SG_ChannelSize_z;
  G4double pos_z = RC->SG_DetectorPosition_z-SG_SideLength_z*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + j*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);

  G4ThreeVector position = G4ThreeVector(pos_x, pos_y, pos_z);
  o_volumes << name_PMTs << "\t" << pos_x << "\t" <<  pos_y << "\t" << pos_z << endl;

  SG_physical_Bottom[i][j] = new G4PVPlacement(
  0,                          // rotation 
  position,                   // position
  SG_logic_Bottom[i][j],    // logical volume
  name,                       // name
  logicEnv,                   // logical volume of mother
  false,                      // boolean operation
  0,                          // copy number
  checkOverlaps);             // check for overlaps

  G4double pos_y_PMT = pos_y - 0.5*RC->SG_ChannelSize_z - 0.5*RC->PMTs_ChannelSize_z;
  G4ThreeVector position_PMT = G4ThreeVector(pos_x,pos_y_PMT,pos_z);

  PMTs_physical_Bottom[i][j] = new G4PVPlacement(
    0,
    position_PMT,
    PMTs_logic_Bottom[i][j],
    name_PMTs,
    logicEnv,
    false,
    0,
    checkOverlaps);


  char surf_name[512];
  sprintf(surf_name, "SG_Surf_Bottom_%d_%d", i, j);
  SG_Surface_Bottom[i][j] = new G4LogicalBorderSurface(surf_name,SG_physical_Bottom[i][j],physEnv,SG_OptSurface);

  char wrapper_name[512];
  sprintf(wrapper_name, "SG_Wrapper_Bottom_%d_%d", i, j);
  SG_Wrapper_Bottom[i][j] = new G4LogicalBorderSurface(wrapper_name,physEnv, SG_physical_Bottom[i][j], SG_OptWrapper);

  char PMT_surf_name[512];
  sprintf(PMT_surf_name, "PMT_InnerSurf_Bottom_%d_%d", i ,j);
  PMT_Inner_Surface_Bottom[i][j] = new G4LogicalBorderSurface(PMT_surf_name, PMTs_physical_Bottom[i][j], physEnv, PMT_InnerSurface);

  char PMT_toScint_surf_name[512];
  sprintf(PMT_toScint_surf_name, "PMT_ToScint_Surf_Right_%d_%d", i ,j);
  PMT_TowardsScintSurface_Bottom[i][j] = new G4LogicalBorderSurface(PMT_toScint_surf_name, PMTs_physical_Bottom[i][j], SG_physical_Bottom[i][j], PMT_TowardsScintSurface);

  char PMT_toPMT_surf_name[512];
  sprintf(PMT_toPMT_surf_name, "PMT_ToPMT_Surf_Bottom_%d_%d", i ,j);
  PMT_TowardsPMTSurface_Bottom[i][j] = new G4LogicalBorderSurface(PMT_toPMT_surf_name, SG_physical_Bottom[i][j], PMTs_physical_Bottom[i][j], PMT_TowardsPMTSurface);

  }
}

// TOP SIDE

G4LogicalVolume* SG_logic_Top[100][100];
G4VPhysicalVolume* SG_physical_Top[100][100];
G4LogicalBorderSurface* SG_Surface_Top[100][100];
G4LogicalBorderSurface* SG_Wrapper_Top[100][100];  // for photons coming from outside. 

G4LogicalVolume* PMTs_logic_Top[100][100];
G4VPhysicalVolume* PMTs_physical_Top[100][100];
G4LogicalBorderSurface* PMT_Inner_Surface_Top[100][100];
G4LogicalBorderSurface* PMT_TowardsScintSurface_Top[100][100];
G4LogicalBorderSurface* PMT_TowardsPMTSurface_Top[100][100];


for (int i = 0; i < SG_SideLength_x; ++i)
{
  for (int j = 0; j < SG_SideLength_z; ++j)
  {
  char name[512];
  sprintf(name, "SG_Top_%d_%d", i, j);
  SG_logic_Top[i][j] = new G4LogicalVolume(solidGammaDetector_TopBottom, SG_mat, name);
  SG_logic_Top[i][j]->SetVisAttributes(SG_visAtt);  
  SG_logic_Top[i][j]->SetSensitiveDetector(sensitiveGammaDetector);

  char name_PMTs[512];
  sprintf(name_PMTs, "PMTs_Top_%d_%d", i, j);  
  PMTs_logic_Top[i][j] = new G4LogicalVolume(solidPMTs_TopBottom, SG_mat, name_PMTs);
  PMTs_logic_Top[i][j]->SetVisAttributes(PMTs_visAtt);
  PMTs_logic_Top[i][j]->SetSensitiveDetector(sensitivePMTDetector);

  G4double pos_x =  (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_x/2 -0.5*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) - i*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);
  G4double pos_y = +(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)*SG_SideLength_y/2 +0.5*RC->SG_ChannelSize_z;
  G4double pos_z = RC->SG_DetectorPosition_z-SG_SideLength_z*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy) + (RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy)/2 + j*(RC->SG_ChannelSize_xy+RC->SG_ChannelGap_xy);

  G4ThreeVector position = G4ThreeVector(pos_x, pos_y, pos_z);
  o_volumes << name_PMTs << "\t" << pos_x << "\t" <<  pos_y << "\t" << pos_z << endl;

  SG_physical_Top[i][j] = new G4PVPlacement(
  0,                          // rotation 
  position,                   // position
  SG_logic_Top[i][j],    // logical volume
  name,                       // name
  logicEnv,                   // logical volume of mother
  false,                      // boolean operation
  0,                          // copy number
  checkOverlaps);             // check for overlaps

  G4double pos_y_PMT = pos_y + 0.5*RC->SG_ChannelSize_z + 0.5*RC->PMTs_ChannelSize_z;
  G4ThreeVector position_PMT = G4ThreeVector(pos_x,pos_y_PMT,pos_z);

    PMTs_physical_Top[i][j] = new G4PVPlacement(
    0,
    position_PMT,
    PMTs_logic_Top[i][j],
    name_PMTs,
    logicEnv,
    false,
    0,
    checkOverlaps);

  char surf_name[512];
  sprintf(surf_name, "SG_Surf_Top_%d_%d", i, j);
  SG_Surface_Top[i][j] = new G4LogicalBorderSurface(surf_name,SG_physical_Top[i][j],physEnv,SG_OptSurface);

  char wrapper_name[512];
  sprintf(wrapper_name, "SG_Wrapper_Top_%d_%d", i, j);
  SG_Wrapper_Top[i][j] = new G4LogicalBorderSurface(wrapper_name,physEnv, SG_physical_Top[i][j], SG_OptWrapper);

  char PMT_surf_name[512];
  sprintf(PMT_surf_name, "PMT_InnerSurf_Top_%d_%d", i ,j);
  PMT_Inner_Surface_Top[i][j] = new G4LogicalBorderSurface(PMT_surf_name, PMTs_physical_Top[i][j], physEnv, PMT_InnerSurface);

  char PMT_toScint_surf_name[512];
  sprintf(PMT_toScint_surf_name, "PMT_ToScint_Surf_Top_%d_%d", i ,j);
  PMT_TowardsScintSurface_Top[i][j] = new G4LogicalBorderSurface(PMT_toScint_surf_name, PMTs_physical_Top[i][j], SG_physical_Top[i][j], PMT_TowardsScintSurface);

  char PMT_toPMT_surf_name[512];
  sprintf(PMT_toPMT_surf_name, "PMT_ToPMT_Surf_Top_%d_%d", i ,j);
  PMT_TowardsPMTSurface_Top[i][j] = new G4LogicalBorderSurface(PMT_toPMT_surf_name, SG_physical_Top[i][j], PMTs_physical_Top[i][j], PMT_TowardsPMTSurface);

  }
}

  SDManager->ListTree();







//####################################################################################//####################################################################################
// SURFACES
//####################################################################################//####################################################################################


// ############ Scintillator Wrapper Surface (outside) ############

  SG_OptWrapper -> SetModel(unified);
  // SG_OptWrapper -> SetFinish(polished);
  SG_OptWrapper -> SetFinish(groundbackpainted);
  //SG_OptWrapper -> SetFinish(polishedbackpainted);
  SG_OptWrapper -> SetType(dielectric_dielectric);
  G4MaterialPropertiesTable* wrapperProperty = new G4MaterialPropertiesTable();

  const G4int wrapper_NUM = 2;
  G4double wrapper_pp[wrapper_NUM] = {1.0*eV, 7.5*eV};
  G4double wrapper_reflectivity[wrapper_NUM] = {0, 0.};
  G4double wrapper_efficiency[wrapper_NUM] = {1.0, 1.0};
  G4double wrapper_Rindex[wrapper_NUM] = {0,0};   // makes sure, no photon enters the scintillator channel from outside.

  wrapperProperty->AddProperty("REFLECTIVITY", wrapper_pp, wrapper_reflectivity, wrapper_NUM);
  wrapperProperty->AddProperty("EFFICIENCY", wrapper_pp, wrapper_efficiency, wrapper_NUM);
  wrapperProperty->AddProperty("RINDEX", wrapper_pp, wrapper_Rindex, wrapper_NUM);

  SG_OptWrapper->SetMaterialPropertiesTable(wrapperProperty);



// ############ Scintillator Inner Surface (inside) ############
// for more Info see Janecek and Moses 2010: 10.1109/TNS.2010.2042731

  SG_OptSurface -> SetModel(LUT);
  SG_OptSurface -> SetFinish(groundteflonair);
  SG_OptSurface -> SetType(dielectric_LUT);


// ############ PMT towards Scintillator ############


  PMT_TowardsScintSurface->SetType(dielectric_dielectric);
  PMT_TowardsScintSurface->SetModel(unified);
  PMT_TowardsScintSurface->SetFinish(groundbackpainted);
  const G4int entranceSurf_NUM = 2;
  G4double entranceSurf_pp[entranceSurf_NUM] = {1.0*eV, 7.5*eV};
  G4double entranceSurf_reflectivity[entranceSurf_NUM] = {0, 0.};
  G4double entranceSurf_efficiency[entranceSurf_NUM] = {0.0, 0.0};
  G4double entranceSurf_Rind[entranceSurf_NUM] = {0,0};
  G4MaterialPropertiesTable* PMT_TowardsScintSurface_MPT = new G4MaterialPropertiesTable();
  PMT_TowardsScintSurface_MPT->AddProperty("RINDEX", entranceSurf_pp, entranceSurf_Rind, entranceSurf_NUM);
  PMT_TowardsScintSurface_MPT->AddProperty("REFLECTIVITY", entranceSurf_pp, entranceSurf_reflectivity, entranceSurf_NUM);
  PMT_TowardsScintSurface_MPT->AddProperty("EFFICIENCY", entranceSurf_pp, entranceSurf_efficiency, entranceSurf_NUM);

  PMT_TowardsScintSurface->SetMaterialPropertiesTable(PMT_TowardsScintSurface_MPT);


// ############ PMT - from PMT to scintillator surface ############

    // not necessary as long as its the same material as the scintillator
/*  
  PMT_TowardsPMTSurface->SetType(dielectric_dielectric);
  PMT_TowardsPMTSurface->SetModel(unified);
  PMT_TowardsPMTSurface->SetFinish(groundbackpainted);
  G4MaterialPropertiesTable* PMT_TowardsPMTSurface_MPT = new G4MaterialPropertiesTable();
  PMT_TowardsPMTSurface_MPT->AddProperty("REFLECTIVITY", entranceSurf_pp, entranceSurf_reflectivity, entranceSurf_NUM);
  PMT_TowardsPMTSurface->SetMaterialPropertiesTable(PMT_TowardsPMTSurface_MPT);

*/

// ############ PMT internal surface ############

  PMT_InnerSurface -> SetModel(unified);
  // PMT_InnerSurface -> SetFinish(polished);
  PMT_InnerSurface -> SetFinish(groundbackpainted);
  //PMT_InnerSurface -> SetFinish(polishedbackpainted);
  PMT_InnerSurface -> SetType(dielectric_dielectric);
  G4MaterialPropertiesTable* PMT_InnerSurface_MPT = new G4MaterialPropertiesTable();

  const G4int innerSurf_NUM = 2;
  G4double PMT_innerSurf_pp[innerSurf_NUM] = {1.0*eV, 7.5*eV};
  G4double PMT_innerSurf_reflectivity[innerSurf_NUM] = {0, 0.};
  G4double PMT_innerSurf_efficiency[innerSurf_NUM] = {1.0, 1.0};
  G4double PMT_innerSurf_Rindex[innerSurf_NUM] = {0,0};   // makes sure, no photon enters the scintillator channel from outside.

  PMT_InnerSurface->SetMaterialPropertiesTable(PMT_InnerSurface_MPT);









//####################################################################################//####################################################################################
// MATERIAL - optical properties
//####################################################################################//####################################################################################

// read in BGO emission properties
  ifstream in_emission;
  ifstream in_transm;
  if(useBGO) in_emission.open("BGO_scint_emission_data.csv");
  else in_emission.open("LaBr3_scint_emission_data.csv");
  
  std::vector<double> wl;
  std::vector<double> intensity;
  std::vector<double> wl_transm;
  std::vector<double> transm;
  double dummy1, dummy2;
  double integrated_intensity=0;

  while(!in_emission.eof()){
    in_emission >> dummy1 >> dummy2;
    integrated_intensity+=dummy2;
    wl.push_back(dummy1);
    intensity.push_back(dummy2);
  }
  cout << "read in emission data from file sucessfully." << endl;
  in_emission.close();
 
  double dummyX, dummyY;
  if(useBGO) in_transm.open("BGO_scint_transmission_data.csv");
  else in_transm.open("LaBr3_scint_attenuation_data.csv");

  while(!in_transm.eof()){
    in_transm >> dummyX >> dummyY;
    wl_transm.push_back(dummyX);
    transm.push_back(dummyY);
  }
  cout << "read in transmission data from file sucessfully." << endl;
  in_transm.close();

  double scint_normalized_intensity[1000];
  double scint_photon_energies[1000];
  double scint_abs_length[1000];
  double scint_photon_energies_transm[1000];

  double h_Planck = 6.626E-34;      // Joule second
  double electronVolt = 1.601E-19;  // Joule
  double c_vacuum = 3.0E8;         // m/s
  double integral=0;

  cout << "Reading in emission properties into arrays..." << endl;
 // cout << "wavelength\tphoton energy\trelative intensity" << endl;
  for (unsigned int i = 0; i < wl.size(); ++i)
  {
    scint_photon_energies[i] = 1/(wl.at(i)*1E-9)*h_Planck*c_vacuum*1/electronVolt*eV; // eV reduces the value by a factor of 10^-6
    scint_normalized_intensity[i] = intensity.at(i)/integrated_intensity;
    integral+=scint_normalized_intensity[i];
 //   cout << wl.at(i) << "\t\t" << scint_photon_energies[i] << "\t" << scint_normalized_intensity[i] << endl;
  }
  cout << "Reading in transmission properties into arrays... "<< endl;
 // cout "wavelength\tphoton energy\tabsortion length" << endl;
  for (unsigned int i = 0; i < wl_transm.size(); ++i)
  {
    scint_photon_energies_transm[i] = 1/(wl_transm.at(i)*1E-9)*h_Planck*c_vacuum*1/electronVolt*eV; // eV reduces the value by a factor of 10^-6
    if(useBGO) scint_abs_length[i] = -(3.24*cm)/log(transm.at(i));
    else scint_abs_length[i] = transm.at(i);
 //   cout << wl_transm.at(i) << "\t\t" << scint_photon_energies_transm[i] << "\t" << scint_abs_length[i] << endl;
  }

  const G4int numDataPoints_BGO = wl.size();
  cout << "---\nTotal number of data points:\t" << numDataPoints_BGO << "\tintegral (should be 1): " << integral << endl;


  G4MaterialPropertiesTable* SG_MPT= new G4MaterialPropertiesTable();
  const G4int scint_n = 6;
  G4double scint_ph_energies[scint_n] = {1.240*eV, 2.029*eV, 2.5*eV, 3.084*eV, 3.542*eV, 3.974*eV};
  G4double scint_rind[scint_n];
  if(useBGO){ 
    scint_rind[0] = 2.056;
    scint_rind[1] = 2.103;
    scint_rind[2] = 2.144;
    scint_rind[3] = 2.22;
    scint_rind[4] =  2.30;
    scint_rind[5] = 2.410;
  //  scint_rind[scint_n] = {2.056, 2.103, 2.144, 2.22, 2.30, 2.410}; //https://refractiveindex.info/?shelf=main&book=Bi4Ge3O12&page=Williams
  }
//  G4double scint_abs_length[scint_n] = {100*cm,100*cm,100*cm,100*cm,100*cm,100*cm};
  else {
   scint_rind[0] = 2.000;
   scint_rind[1] = 2.050;
   scint_rind[2] = 2.110;
   scint_rind[3] = 2.240;
   scint_rind[4] = 2.370;
   scint_rind[5] = 2.42;      // data taken from picture in the simulation folder!
  }
  G4double scint_norm_intens[scint_n] = {0.1,0.3,0.2,0.2,0.1,0.1};
  SG_MPT->AddProperty("RINDEX", scint_ph_energies, scint_rind, scint_n);
  SG_MPT->AddProperty("ABSLENGTH", scint_photon_energies_transm, scint_abs_length, wl_transm.size());

  SG_MPT->AddProperty("FASTCOMPONENT", scint_photon_energies, scint_normalized_intensity, numDataPoints_BGO);
//  SG_MPT->AddProperty("FASTCOMPONENT", scint_ph_energies, scint_norm_intens, scint_n);  
  SG_MPT->AddConstProperty("SCINTILLATIONYIELD", RC->ScintYield/MeV);
  SG_MPT->AddConstProperty("FASTTIMECONSTANT", 300*ns);
  //SG_MPT->AddConstProperty("SLOWTIMECONSTANT", 301*ns);
  SG_MPT->AddConstProperty("YIELDRATIO", 0.1);
  SG_MPT->AddConstProperty("RESOLUTIONSCALE",1.0);

  SG_mat->SetMaterialPropertiesTable(SG_MPT);









//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PLACEMENT OF PMMA AND SIPMs FOR CHERENKOV DETECTOR
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //
  // PMMA
  //

 G4Box* solidPMMA= 
    new G4Box("PMMA", 0.5*PMMA_size_x, 0.5*PMMA_size_y, 0.5*PMMA_thick);

  G4LogicalVolume* logicPMMA = 
    new G4LogicalVolume(solidPMMA, PMMA_mat,"PMMA");
   //new G4LogicalVolume(solidPMMA, BGO,"PMMA");

  G4VPhysicalVolume* physPMMA = 
    new G4PVPlacement(0, G4ThreeVector(0,0,0), logicPMMA, "PMMA", logicEnv, false, 0, checkOverlaps);

 G4OpticalSurface*OpSurface = new G4OpticalSurface("PMMA_surface");
 G4LogicalBorderSurface*Surface = new G4LogicalBorderSurface("PMMA_surface",physPMMA,physEnv,OpSurface);
 OpSurface -> SetFinish(groundbackpainted);



 G4LogicalVolume* logicVolumes_ref[100][100];
 char name_ref[1000];


 for(int k=0; k< p; ++k){
	for(int l=0; l<w; ++l){
		 for(int i=0; i < n;i++){
 			for(int j=0; j< v;j++){
			sprintf(name_ref,"%s_%d_%d","ref",k*n+i,l*v+j);

			G4ThreeVector ref_i_j_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),	
			 +PMMA_size_y/2-SiPMchannel_size_xy/2-space_y-j*(SiPMchannel_size_xy+space)-l*(((v-1)*space+v*SiPMchannel_size_xy)+spacem), 
			(PMMA_thick/2+ref_thick/2)
			);

  			G4Box* solidref_i_j = 
   			new G4Box(name_ref, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

  			logicVolumes_ref[i][j] = 
   			new G4LogicalVolume(solidref_i_j, ref_mat, name_ref);
  			new G4PVPlacement(0,ref_i_j_pos, logicVolumes_ref[i][j],name_ref, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }

 
 //SiPM_stop



 G4LogicalVolume* logicVolumes[100][100];
 char name[1000];


 for(int k=0; k< p; ++k){
	for(int l=0; l<w; ++l){
		 for(int i=0; i < n;i++){
 			for(int j=0; j< v;j++){
			sprintf(name,"%d_%d", k*n+i,l*v+j);

			G4ThreeVector SiPM_i_j_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),
			 +PMMA_size_y/2-SiPMchannel_size_xy/2-space_y-j*(SiPMchannel_size_xy+space)-l*(((v-1)*space+v*SiPMchannel_size_xy)+spacem), 
			PMMA_thick/2+ref_thick/2+ref_thick
			);

  			G4Box* solid_i_j = 
   			new G4Box(name, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

  			logicVolumes[i][j] = 
   			new G4LogicalVolume(solid_i_j, SiPM_mat, name);
  			new G4PVPlacement(0,SiPM_i_j_pos, logicVolumes[i][j],name, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }





 //Boundary SiPM


 //SiPM_ref, right

 G4LogicalVolume* logicVolumes_refright[100][100];
 char name_refright[1000];

 G4RotationMatrix rotSiPM_right  = G4RotationMatrix();
    rotSiPM_right.rotateY(90*deg); 
    rotSiPM_right.rotateZ(0*deg);
    rotSiPM_right.rotateX(0*deg);

 for(int k = 0; k < z1; ++k){					
 	for(int l = 0; l< w; ++l){
		 for(int i = 0; i < z2; i++){			
			for(int j = 0; j < v; j++){
			sprintf(name_refright,"%s_%d_%d","refright",n*p+k*z2+i,l*v+j);   

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-PMMA_size_x/2-ref_thick/2,
			-space_y+PMMA_size_y/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(spacem+(v-1)*space+v*SiPMchannel_size_xy),
			 PMMA_thick/2-SiPMchannel_size_xy/2-i*(SiPMchannel_size_xy+space)-k*((z2-1)*space+z2*SiPMchannel_size_xy+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_right,SiPM_sides_pos);

			G4Box* solidrefright_i_j = 
			new G4Box(name_refright, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_refright[i][j] = 
			new G4LogicalVolume(solidrefright_i_j, ref_mat, name_refright);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_refright[i][j],name_refright, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }



 //SiPM_stop, right

 G4LogicalVolume* logicVolumes_right[100][100];
 char name_right[1000];

 for(int k = 0; k < z1; ++k){					
 	for(int l = 0; l< w; ++l){
		 for(int i = 0; i < z2; i++){			
			for(int j = 0; j < v; j++){
			sprintf(name_right,"%s_%d_%d","right",n*p+k*z2+i,l*v+j);   

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-PMMA_size_x/2-ref_thick/2-ref_thick,
			-space_y+PMMA_size_y/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(spacem+(v-1)*space+v*SiPMchannel_size_xy),
			 PMMA_thick/2-SiPMchannel_size_xy/2-i*(SiPMchannel_size_xy+space)-k*((z2-1)*space+z2*SiPMchannel_size_xy+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_right,SiPM_sides_pos);

			G4Box* solidright_i_j = 
			new G4Box(name_right, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_right[i][j] = 
			new G4LogicalVolume(solidright_i_j, SiPM_mat, name_right);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_right[i][j],name_right, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }



 //SiPM_ref, left

 G4LogicalVolume* logicVolumes_refleft[100][100];
 char name_refleft[1000];

 G4RotationMatrix rotSiPM_left  = G4RotationMatrix();
    rotSiPM_left.rotateY(270*deg); 
    rotSiPM_left.rotateZ(0*deg);
    rotSiPM_left.rotateX(0*deg);

 for(int k = 0; k < z1; ++k){
	for(int l = 0; l < w; ++l){
		 for(int i = 0; i < z2; i++){
			for(int j = 0; j < v; j++){
			sprintf(name_refleft,"%s_%d_%d","refleft",-k*z2-i-1,l*v+j);			

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			PMMA_size_x/2+ref_thick/2,
			-space_y+PMMA_size_y/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(spacem+(v-1)*space+v*SiPMchannel_size_xy),
			PMMA_thick/2-SiPMchannel_size_xy/2-i*(SiPMchannel_size_xy+space)-k*((z2-1)*space+z2*SiPMchannel_size_xy+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_left,SiPM_sides_pos);

			G4Box* solidrefleft_i_j = 
			new G4Box(name_refleft, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_refleft[i][j] = 
			new G4LogicalVolume(solidrefleft_i_j, ref_mat, name_refleft);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_refleft[i][j],name_refleft, logicEnv, false, 0, checkOverlaps);
			}
		 }		
	}
 }

 //SiPM_stop, left

 G4LogicalVolume* logicVolumes_left[100][100];
 char name_left[1000];

 for(int k = 0; k < z1; ++k){
	for(int l = 0; l < w; ++l){
		 for(int i = 0; i < z2; i++){
			for(int j = 0; j < v; j++){
			sprintf(name_left,"%s_%d_%d","left",-k*z2-i-1,l*v+j);			

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			PMMA_size_x/2+ref_thick/2+ref_thick,
			-space_y+PMMA_size_y/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(spacem+(v-1)*space+v*SiPMchannel_size_xy),
			PMMA_thick/2-SiPMchannel_size_xy/2-i*(SiPMchannel_size_xy+space)-k*((z2-1)*space+z2*SiPMchannel_size_xy+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_left,SiPM_sides_pos);

			G4Box* solidleft_i_j = 
			new G4Box(name_left, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_left[i][j] = 
			new G4LogicalVolume(solidleft_i_j, SiPM_mat, name_left);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_left[i][j],name_left, logicEnv, false, 0, checkOverlaps);
			}
		 }		
	}
 }




 //SiPM_ref, top

 G4LogicalVolume* logicVolumes_refup[100][100];
 char name_refup[1000];

 G4RotationMatrix rotSiPM_up  = G4RotationMatrix();
    rotSiPM_up.rotateY(0*deg); 
    rotSiPM_up.rotateZ(0*deg);
    rotSiPM_up.rotateX(90*deg);

 for(int k = 0; k < p; ++k){
	for(int l = 0; l < z1; ++l){
		 for(int i = 0; i < n; i++){
			for(int j = 0; j < z2; j++){
			sprintf(name_refup,"%s_%d_%d","refup",k*n+i,-l*z2-j-1);		

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),
			+PMMA_size_y/2+ref_thick/2,
			PMMA_thick/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(((z2-1)*space+z2*SiPMchannel_size_xy)+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_up,SiPM_sides_pos);

			G4Box* solidrefup_i_j = 
			new G4Box(name_refup, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_refup[i][j] = 
			new G4LogicalVolume(solidrefup_i_j, ref_mat, name_refup);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_refup[i][j],name_refup, logicEnv, false, 0, checkOverlaps);
			}
 		}
	}
 }

 //SiPM_stop, top

 G4LogicalVolume* logicVolumes_up[100][100];
 char name_up[1000];

 for(int k = 0; k < p; ++k){
	for(int l = 0; l < z1; ++l){
		 for(int i = 0; i < n; i++){
			for(int j = 0; j < z2; j++){
			sprintf(name_up,"%s_%d_%d","up",k*n+i,-l*z2-j-1);

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),
			+PMMA_size_y/2+ref_thick/2+ref_thick,
			PMMA_thick/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(((z2-1)*space+z2*SiPMchannel_size_xy)+spacem)
			);

			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_up,SiPM_sides_pos);

			G4Box* solidup_i_j = 
			new G4Box(name_up, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_up[i][j] = 
			new G4LogicalVolume(solidup_i_j, SiPM_mat, name_up);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_up[i][j],name_up, logicEnv, false, 0, checkOverlaps);
			}
 		}
	}
 }


 //SiPM_ref, bottom

 G4LogicalVolume* logicVolumes_refdown[100][100];
 char name_refdown[1000];

 G4RotationMatrix rotSiPM_down  = G4RotationMatrix();
    rotSiPM_down.rotateY(0*deg); 
    rotSiPM_down.rotateZ(0*deg);
    rotSiPM_down.rotateX(90*deg);

 for(int k = 0; k < p; ++k){
	for(int l = 0; l < z1; ++l){
		 for(int i = 0; i < n; i++){
			for(int j = 0; j < z2; j++){
			sprintf(name_refdown,"%s_%d_%d","refdown",k*n+i,w*v+l*z2+j);			

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),
			-PMMA_size_y/2-ref_thick/2,
			PMMA_thick/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(((z2-1)*space+z2*SiPMchannel_size_xy)+spacem)
			);
	
			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_down,SiPM_sides_pos);

			G4Box* solidrefdown_i_j = 
			new G4Box(name_refdown, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_refdown[i][j] = 
			new G4LogicalVolume(solidrefdown_i_j, ref_mat, name_refdown);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_refdown[i][j],name_refdown, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }

 //SiPM_stop, bottom

 G4LogicalVolume* logicVolumes_down[100][100];
 char name_down[1000];


 for(int k = 0; k < p; ++k){
	for(int l = 0; l < z1; ++l){
		 for(int i = 0; i < n; i++){
			for(int j = 0; j < z2; j++){
			sprintf(name_down,"%s_%d_%d","down",k*n+i,w*v+l*z2+j);

			G4ThreeVector SiPM_sides_pos = 
			G4ThreeVector(
			-1*(-PMMA_size_x/2+SiPMchannel_size_xy/2+space_x+i*(SiPMchannel_size_xy+space)+k*(((n-1)*space+n*SiPMchannel_size_xy)+spacem)),
			-PMMA_size_y/2-ref_thick/2-ref_thick,
			PMMA_thick/2-SiPMchannel_size_xy/2-j*(SiPMchannel_size_xy+space)-l*(((z2-1)*space+z2*SiPMchannel_size_xy)+spacem)
			);
	
			G4Transform3D transformSiPM_sides = G4Transform3D(rotSiPM_down,SiPM_sides_pos);

			G4Box* soliddown_i_j = 
			new G4Box(name_down, 0.5*SiPMchannel_size_xy, 0.5*SiPMchannel_size_xy, 0.5*ref_thick);

			logicVolumes_down[i][j] = 
			new G4LogicalVolume(soliddown_i_j, SiPM_mat, name_down);
			new G4PVPlacement(transformSiPM_sides, logicVolumes_down[i][j],name_down, logicEnv, false, 0, checkOverlaps);
			}
		 }
	}
 }


/////////////////////////////////////////////////////////////////////
// assign optical properties to the materials of the SiPM channels
/////////////////////////////////////////////////////////////////////
const int NUMENTRIES = 701;
G4double photonEnergy[NUMENTRIES];
G4double abslength[NUMENTRIES];

ifstream i_abslength;
i_abslength.open("Photon_energy_and_absorption_length.csv");
double d1, d2;
int ct=0;
while(!i_abslength.eof()){
i_abslength >> d1 >> d2;
photonEnergy[ct]=d1*eV;
abslength[ct] = d2*mm;
ct++;
}

 G4double phenergy[]={1.3776*eV, 6.2280*eV};
 G4double refractiveIndexPMMA[] ={1.49,1.50};
 G4double refractiveIndexAir[] = {1.00,1.00};
 G4double refractiveIndexSiPM[] = {1.55,1.55};
 //G4double abslengthSiPM[] = {35*mm,35*mm};			/////
 const G4int nEntries = sizeof (phenergy) / sizeof (G4double);
 G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX", phenergy, refractiveIndexPMMA, nEntries) 
  ->SetSpline(true);

  myMPT1->AddProperty("ABSLENGTH", photonEnergy, abslength, NUMENTRIES) 
  ->SetSpline(true);

 G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();
  myMPT3->AddProperty("RINDEX", phenergy, refractiveIndexSiPM, nEntries) 
  ->SetSpline(true);

 //myMPT3->AddProperty("ABSLENGTH", phenergy, abslengthSiPM, nEntries) 
 //->SetSpline(true);

 G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", phenergy, refractiveIndexAir, nEntries) 
  ->SetSpline(true);

/////////////////////////////////////////////////////////////////////
// comment out these two lines to disable reflection and transmission at PMMA surfaces.
// all optical photon tracks are then terminated at the border
/////////////////////////////////////////////////////////////////////


 PMMA_mat->SetMaterialPropertiesTable(myMPT1); 
 env_mat->SetMaterialPropertiesTable(myMPT2);
 world_mat->SetMaterialPropertiesTable(myMPT2);
 ref_mat->SetMaterialPropertiesTable(myMPT3);

 //SiPM_mat->SetMaterialPropertiesTable(myMPT3);

  o_volumes.close();

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
