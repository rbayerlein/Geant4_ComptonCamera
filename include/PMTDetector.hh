#ifndef PMTDetector_h
#define PMTDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"
#include "G4HCofThisEvent.hh"
#include "G4HCtable.hh"
#include "HistoManager.hh"
#include "PMTHit.hh"

#include <string>
#include <iostream>
#include <fstream>

class G4Step;

class PMTDetector : public G4VSensitiveDetector
{
public:
	PMTDetector(std::string name);

	virtual ~PMTDetector();

	virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

	/* (optional) methods */
	void Initialize(G4HCofThisEvent* kHCEvent);
	void EndOfEvent(G4HCofThisEvent* kHCEvent);


private:
	/* Pointer to collection of Hits */
	hitCollection* gHitCollection;
	G4int mCollectionID = -1;
	

//	G4AnalysisManager* analysisManager;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif