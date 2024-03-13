/// \file SGDetector.cc
/// \brief Implementation of the Scattered Gamma Detector class

#include "SGDetector.hh"
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

SGDetector::SGDetector(string name)
: G4VSensitiveDetector(name)
{
 collectionName.insert(name); // 'collectionName' is protected member of SDollectionName.insert(SDname); // 'collectionName' is protected member of SD

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SGDetector::~SGDetector(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SGDetector::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist){

	Analysis* ana = Analysis::GetInstance();

	string particleName = aStep->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
	if(particleName == "gamma"){
	//	G4String volname = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

	//	cout << "Particle in scoring volume " << volname << ". Particle type: " << particleName << endl;
	//  G4ThreeVector SG_HitPosition = aStep->GetPreStepPoint()->GetPosition();
	//	cout << "at position (" << SG_HitPosition.x() << "," << SG_HitPosition.y() << "," << SG_HitPosition.z() << ")" << endl;


	}
	if(particleName == "e-"){
		if(ana->GetComptonCounter() < 1) return true; 	//only process events in the SG detector if there have been Compton events in the Chkv detector in this event.
		const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
		int trackID = aStep->GetTrack()->GetTrackID();
		//cout << "----\nElectron found: Track ID: " << trackID << ". Number of secondaries produced: " << secondaries->size() << endl;
		for (unsigned int i = 0; i < secondaries->size(); ++i)
		{	
		//	cout << secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() << endl;
			if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "opticalphoton"){
				G4double phenergy = secondaries->at(i)->GetKineticEnergy();
				G4double phlambda = 1.2398E-6*eV*m/phenergy*1E6;// wavelength [nm]
			//	cout << "photon found with wl: " << phlambda << endl;
				
				ana->FillScintVector(phlambda);
			}
		}

	}

	
	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
