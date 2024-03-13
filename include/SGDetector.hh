#ifndef SGDetector_h
#define SGDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"
#include "HistoManager.hh"

#include <string>
#include <iostream>
#include <fstream>

class G4Step;

class SGDetector : public G4VSensitiveDetector
{
public:
	SGDetector(std::string name);

	virtual ~SGDetector();

	virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

//	G4AnalysisManager* analysisManager;
	
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

