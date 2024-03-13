/// \file Analysis.hh
/// \brief Definition of the Analysis class

#ifndef Analysis_h
#define Analysis_h 1

#include "globals.hh"
#include "G4Track.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"

#include <string>
#include <iostream>
#include <fstream>

/// Analysis class
/// 
using namespace std;


class B1DetectorConstruction;
class Analysis
{
  public:
	Analysis(B1DetectorConstruction*);
	static Analysis* GetInstance(){// Returns SINGLETON of Analysis
								// i.e. it only exists once
		if (Analysis::fInstance == NULL) {
			Analysis::fInstance = new Analysis();
			}
		return fInstance;
	}
	virtual ~Analysis();
	
	void PrepareNewRun(const G4Run* run);
	void EndOfRunAction(const G4Run* run);

	void PrepareNewEvent(const G4Event* event);
	void EndOfEventAction(const G4Event* event);
	
	void PreUserTrackingAction(const G4Track* track);
	void PostUserTrackingAction(const G4Track* track);
	
	void SteppingAction(const G4Step* step);

	void FillScintVector(double in){SG_Scint_Wavelengths.push_back(in);}

	int GetComptonCounter(){return fComptonCounter;};

	vector<double> wl_PDE;	// wavelength from the PDE file
	vector<double> PDE;	// corresponing PDE
	vector<double> PMT_scint_wavelengths;
	std::vector<double> SG_Scint_Wavelengths;

  private:
	Analysis(); // Private Constructor
	G4VPhysicalVolume*	D1;  
	G4VPhysicalVolume*	D2; 
	G4VPhysicalVolume*	D3; 
	G4VPhysicalVolume*	D4; 
	G4VPhysicalVolume*	D5; 
	G4VPhysicalVolume*	D6; 
	G4VPhysicalVolume*	D7; 
	G4VPhysicalVolume*	D8; 
	double in_energy;
	static Analysis * fInstance;
	int lastDist;
	int lastTrackID;
	int LastTrackID;
	int fCerenkovCounter;
	int fComptonCounter = 0;
	int fPhotoCounter = 0;
	bool savedElectronStartPosition;
	int fEventNumber;
	int detpn;
	//	G4LogicalVolume* fScoringVolume;
		G4LogicalVolume* D4ScoringVolume;
		G4LogicalVolume* D5ScoringVolume;

	G4double photonenergy=1.5;
	G4double eminusenergy=200;

	G4double electronEnergy;
	G4ThreeVector electronMomentumDirection;
	G4ThreeVector electronVertex;

	int nofchannels;
	int nofphotons;

	int nofchannelsarray;
	int nofchannelstop;
	int nofchannelsbottom;
	int nofchannelsleft;
	int nofchannelsright;
	int nofchannelstest;

	int nofchannelsarraynoise;
	int nofchannelstopnoise;
	int nofchannelsbottomnoise;
	int nofchannelsleftnoise;
	int nofchannelsrightnoise;
	int nofchannelstestnoise;

 
	int nofphotonsarray;
	int nofphotonstop;
	int nofphotonsbottom;
	int nofphotonsleft;
	int nofphotonsright;
	int nofphotonstest;

	int nofphotonsarraynoise;
	int nofphotonstopnoise;
	int nofphotonsbottomnoise;
	int nofphotonsleftnoise;
	int nofphotonsrightnoise;
	int nofphotonstestnoise;

	int nofphotons_sig_and_noise;
	int nofchannels_sig_and_noise;

	int n;
	int v;
	int p;

	int channelHit[256][256];  
	int channelHittest[256][256]; 
	int channelHitarray[256][256];
	int channelHitleft[256][256]; 
	int channelHitright[256][256]; 	
	int channelHittop[256][256]; 	
	int channelHitbottom[256][256];

 	int channelHitnoise[256][256];  
	int channelHittestnoise[256][256]; 
	int channelHitarraynoise[256][256];
	int channelHitleftnoise[256][256]; 
	int channelHitrightnoise[256][256]; 	
	int channelHittopnoise[256][256]; 	
	int channelHitbottomnoise[256][256];

	int channelNOP[256][256];
	int channelNOPtest[256][256];
	int channelNOParray[256][256];
	int channelNOPleft[256][256];
	int channelNOPright[256][256];
	int channelNOPtop[256][256];
	int channelNOPbottom[256][256];

	int channelNOPnoise[256][256];
	int channelNOPtestnoise[256][256];
	int channelNOParraynoise[256][256];
	int channelNOPleftnoise[256][256];
	int channelNOPrightnoise[256][256];
	int channelNOPtopnoise[256][256];
	int channelNOPbottomnoise[256][256];



string text_file;
ofstream dat_out;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
