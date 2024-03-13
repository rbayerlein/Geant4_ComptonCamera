/// \file PMTDetector.cc
/// \brief Implementation of the Scattered Gamma Detector class

#include "PMTDetector.hh"
#include "Analysis.hh"
#include "G4Track.hh"
#include "B1DetectorConstruction.hh"
#include "HistoManager.hh"
#include "RunControl.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4VSensitiveDetector.hh"


using namespace std;

PMTDetector::PMTDetector(string name)
: G4VSensitiveDetector(name)
{
 collectionName.insert(name); // 'collectionName' is protected member of SDollectionName.insert(SDname); // 'collectionName' is protected member of SD
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PMTDetector::~PMTDetector(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PMTDetector::Initialize(G4HCofThisEvent* kHCEvent) /* create new Hit Collection at start of the event */
{
   /* create a collection to store info about the hits */
   gHitCollection = new hitCollection(SensitiveDetectorName, collectionName[0]);
 
   /* To insert the collection, we need to get an unique ID for it */
   if(mCollectionID<0) mCollectionID = GetCollectionID(0); // <<-- this is to get an ID for collectionName[0]
 
   /* add the info into the collection */
   kHCEvent->AddHitsCollection(mCollectionID, gHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PMTDetector::ProcessHits(G4Step* kStep, G4TouchableHistory* ){



	string particleName = kStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
	if(particleName == "gamma"){
	//	cout << "gamma in PMT" << endl;
	}
	if(particleName == "opticalphoton"){
		G4String volname = kStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
	//	cout << "optical photon in PMT " << volname << endl;
	}


   	/* create a hit and populate it with the information */
	PMTHit* hit = new PMTHit(mCollectionID);
  	hit->setName        (kStep->GetTrack()->GetParticleDefinition()->GetParticleName());
   	hit->setEdep        (kStep->GetTotalEnergyDeposit());
   	hit->setGlobalTime  (kStep->GetPreStepPoint()->GetLocalTime());
    hit->setTrackID     (kStep->GetTrack()->GetTrackID());
    hit->setVolume      (kStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName());
    hit->setEnergy      (kStep->GetTrack()->GetKineticEnergy());
    // track ID 
    // energy / wavelength
    // volume
 
  	/* store this hit in the collection, basically a push_back of the HC vector */
   	gHitCollection -> insert(hit);
   	// cout << "hit inserted: # " << mCollectionID << " particle " << hit->getName() << endl;

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PMTDetector::EndOfEvent(G4HCofThisEvent*)/* Action taken at the end of the event, but before EventAction */
{
    // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  Analysis* analysis = Analysis::GetInstance();

   /* no point continuing if there were no hits */
   int n = gHitCollection -> entries();
   if(!n) return;

   if(n == 1 && (*gHitCollection)[0]->getName() == "gamma") return; //skip writing this event if only a gamma entered the channel

   int savedID = -1;
   int photonCounter = 0;
   std::vector<G4String> names;
   G4int NtupleSColumn = 1;  
   G4int NtupleSColumn_max = 3;   

   int photonsPerVolume[3]={0};

   for (int i = 0; i < n; ++i)
   {
    if((*gHitCollection)[i]->getName() == "opticalphoton"){
      if(savedID == (*gHitCollection)[i]->getTrackID() ) continue;
      savedID = (*gHitCollection)[i]->getTrackID();
      photonCounter++;
      // This saves the wavelengths of all detected photons per event. Would only uncomment for small number of events.
      // analysis->PMT_scint_wavelengths.push_back(1.2398E-6*eV*m/(*gHitCollection)[i]->getEnergy()*1E6);
      G4String currentName = (*gHitCollection)[i]->getVolume();
      bool found = false;

      int volPos=0;
      for (unsigned int j = 0; j < names.size(); ++j) // check if current volume has already been added. 
      {
        if(currentName == names.at(j)) {
          found = true; 
          volPos=j;
          break;}
      }
      if(found) photonsPerVolume[volPos]++;

      if(!found && NtupleSColumn <= NtupleSColumn_max) { // add current volume if it is a new volume
//        cout << currentName << endl;
        names.push_back(currentName);
        analysisManager->FillNtupleSColumn(NtupleSColumn, (*gHitCollection)[i]->getVolume());
        NtupleSColumn++;
        photonsPerVolume[names.size()-1]++;
      }
    }
/*     if((*gHitCollection)[i]->getName() == "gamma") {
      cout  << (*gHitCollection)[i]->getName() << "\t" 
            << (*gHitCollection)[i]->getTrackID() << "\t"
            << (*gHitCollection)[i]->getVolume() << "\t"
            << (*gHitCollection)[i]->getEnergy() << endl;
    }*/
   }
//   cout << "Scintillator: Number of detected photons: " << photonCounter << endl;
   analysisManager->FillNtupleIColumn(0, 0, photonCounter);
   analysisManager->FillNtupleIColumn(0, 5, photonsPerVolume[0]);
   analysisManager->FillNtupleIColumn(0, 6, photonsPerVolume[1]);
   analysisManager->FillNtupleIColumn(0, 7, photonsPerVolume[2]);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
