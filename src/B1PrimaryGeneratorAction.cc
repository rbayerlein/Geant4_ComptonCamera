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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"
#include "RunControl.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4VPrimaryGenerator.hh" 
#include "G4ParticleGun.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),																																				
  fParticleGun(0), 
  fEnvelopeBox(0),
  fPMMABox(0)
{
  // fParticleGun  = new G4GeneralParticleSource();
    fParticleGun = new G4ParticleGun(1);
  
  RunControl* RC = RunControl::GetInstance();
  // default particle kinematic-TEST
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle																		
    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);

 fParticleGun->SetParticleEnergy(RC->ParticleEnergy);
// fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-6.5*cm));

  // default particle kinematic
/*  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleEnergy(1.5*MeV);*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get PMMA volume
  // from G4LogicalVolumeStore.
 
  RunControl* RC = RunControl::GetInstance();

  if(RC->coverPMMASolidAngle){

    G4double PMMA_SizeXY = 0.0;

    if (!fPMMABox)
    {
      G4LogicalVolume* PMMA_LV = G4LogicalVolumeStore::GetInstance()->GetVolume("PMMA");
      if ( PMMA_LV ) fPMMABox = dynamic_cast<G4Box*>(PMMA_LV->GetSolid());
    }

    if ( fPMMABox ) 
    {
      PMMA_SizeXY = fPMMABox->GetXHalfLength()*2.;
    }  
    else  {
      G4ExceptionDescription msg;
      msg << "PMMA volume of box shape not found.\n"; 
      msg << "Perhaps you have changed geometry.\n";
      msg << "The gun will be place at the center. Maybe...";
      G4Exception("B1PrimaryGeneratorAction::GeneratePrimaries()",
       "MyCode0002",JustWarning,msg);
    }

// choose point on line in x-direction

      G4double lengthOfLine = RC->lengthOfLine;
      G4double pointOnLine = G4UniformRand()*lengthOfLine; 

// choose start position
     //ausgedehnte quelle
      
      double pi=CLHEP::pi;
      G4double phi = G4UniformRand()*2*pi;
      G4double theta = acos(G4UniformRand()*2-1);
      G4double ra = RC->r_source;
      G4ThreeVector source_shift = RC->source_shift;

     G4double PMMA_thick=RC->PMMA_thick;
     G4double dist_PMMAap = RC->dist_PMMAap;      // from aperture to center of PMMA
     G4double dist_apsource = RC->dist_apsource;

        G4double z = ra* sin(theta)*cos(phi);
      G4double y = ra*sin(theta)*sin(phi);
      G4double x = ra*cos(theta);
      double startposition_x = - lengthOfLine/2 + x + pointOnLine;
      double startposition_y = y;
      double startposition_z = z- dist_apsource- dist_PMMAap;
      G4ThreeVector startposition =  G4ThreeVector(startposition_x, startposition_y, startposition_z) + source_shift;
      fParticleGun->SetParticlePosition(startposition);

// choose momentum direction
      G4double rn_x = G4UniformRand();
      G4double rn_y = G4UniformRand();

      G4double point_on_PMMA_x = - PMMA_SizeXY/2 + rn_x * PMMA_SizeXY;
      G4double point_on_PMMA_y = - PMMA_SizeXY/2 + rn_y * PMMA_SizeXY;

      G4ThreeVector point_on_PMMA = G4ThreeVector(point_on_PMMA_x, point_on_PMMA_y, -PMMA_thick/2);

      G4ThreeVector startmomentum = point_on_PMMA - startposition;
      fParticleGun->SetParticleMomentumDirection(startmomentum);

  }else{

   G4double PMMA_thick=RC->PMMA_thick;
   G4double dist_PMMAap = RC->dist_PMMAap;      // from aperture to center of PMMA
   G4double dist_apsource = RC->dist_apsource;
   float alpha = RC->alpha;

   //ausgedehnte quelle
      double pi=CLHEP::pi;
      G4double phi = G4UniformRand()*2*pi;
      G4double theta = acos(G4UniformRand()*2-1);
      G4double ra = RC->r_source;
      G4ThreeVector source_shift = RC->source_shift;

      G4double z = ra* sin(theta)*cos(phi);
      G4double y = ra*sin(theta)*sin(phi);
      G4double x = ra*cos(theta);
      if(y>0){y=-y;}

    // G4ThreeVector startposition = G4ThreeVector((x),(y),(z-(dist_PMMAap+dist_apsource+PMMA_thick/2)));
    // G4ThreeVector startpositionalpha = G4ThreeVector(startposition.x()*cos(alpha*pi/180)+startposition.z()*sin(alpha*pi/180),startposition.y(),-startposition.x()*sin(alpha*pi/180)+startposition.z()*cos(alpha*pi/180));
    double startpositionalpha_x = x+sin(alpha*pi/180)*(dist_PMMAap - PMMA_thick/2 + dist_apsource);
    double startpositionalpha_y = y;
    double startpositionalpha_z = z-cos(alpha*pi/180)*(dist_PMMAap- PMMA_thick/2 +dist_apsource)-PMMA_thick/2;

    G4ThreeVector startpositionalpha = G4ThreeVector(startpositionalpha_x,startpositionalpha_y,startpositionalpha_z) + source_shift;
    fParticleGun->SetParticlePosition(startpositionalpha);

   float d = RC->d_ap;    // diameter of the aperture 
   float r;
   float x_i=0;
   float y_i=0;

   bool startmomentum =false;
   
   while(startmomentum==false){
   	  float r1 = G4UniformRand();
  	  float r2 = G4UniformRand();
         	  x_i= d*r1-0.5*d;				
  	  y_i= d*r2-0.5*d;
  	  r= sqrt(x_i*x_i+y_i*y_i);
  	  if(r<=0.5*d){startmomentum = true;}
   }

     // vector from PMMA to the aperture for alpha = 0 deg 
        double z_i = -dist_PMMAap+PMMA_thick/2;
        G4ThreeVector aperture_0deg = G4ThreeVector(x_i,y_i,z_i);

    // rotate vector
        double aperture_x = cos(-alpha/180*pi)*x_i + sin(-alpha/180*pi)*z_i;
        double aperture_y = y_i;
        double aperture_z = -sin(-alpha/180*pi)*x_i + cos(-alpha/180*pi)*z_i;
    // add pmma thickness to z component:
        aperture_z = aperture_z -PMMA_thick/2;

     // G4ThreeVector apreture = G4ThreeVector((x_i+0*mm),(y_i+0*mm), -(dist_PMMAap))+source_shift;
  /*   double aperture_x = x_i+sin(alpha*pi/180)*dist_PMMAap;
     double aperture_y = y_i;
     double aperture_z = -PMMA_thick/2-cos(alpha/180*pi)*dist_PMMAap;*/
     G4ThreeVector apreture = G4ThreeVector(aperture_x, aperture_y, aperture_z)+source_shift;

     G4ThreeVector startmom = apreture - startpositionalpha ;
     // RC->startmomalpha = G4ThreeVector(startmom.x()*cos(alpha*pi/180)+startmom.z()*sin(alpha*pi/180),startmom.y(),-startmom.x()*sin(alpha*pi/180)+startmom.z()*cos(alpha*pi/180));
     fParticleGun->SetParticleMomentumDirection(startmom);
  }
// G4cout << "countedPhotons" << "\t" << RC->primaryphotonscounter <<  endl;
 
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));   
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

