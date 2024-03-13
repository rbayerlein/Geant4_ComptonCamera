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
//
// $Id: HistoManager.cc 72239 2013-07-12 08:41:30Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "RunControl.hh"
#include <PMTDetector.hh>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/////////////////////////////////////////////
// here you can change the file name for the root output file
// old files will be overwritten, so make sure to either change the name or copy your old files
// to a safe location
// it would be nice if the file name would include a timestamp to prevent this,
// but I don't know how to do that :(
////////////////////////////////////////////
time_t now = time(0);
char* dt = ctime(&now);

HistoManager::HistoManager(string fileName)
  : fFileName(fileName)

{
	/*RunControl* RC = RunControl::GetInstance();
	fFileName=RC->fileName;*/
  	Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  cout<< "HistoManager:\tSet \n";
  analysisManager->SetFileName(fFileName); 	cout << "file name " << fFileName << "\n";
  analysisManager->SetVerboseLevel(1);		cout << "verbose level 1\n";
  analysisManager->SetActivation(true);   	cout << "activation true\n";	//enable inactivation of histograms
 



///////////////////////////////////////////
// Define 1D histogram titles
// If you want to add a new histogram, increment kMaxHisto by one, add another id to the string and define the title.
// Don't forget to use a comma to separate them from previous parts of the string
///////////////////////////////////////////

  const G4int kMaxHisto = 39;
  const G4String id[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30","31","32","33","34", "35", "36", "37", "38"};
  const G4String title[] = 
                { "comptoncounter",  //"scattered primary particle: energy spectrum",        //0
                  "number of compton events caused by a single photon", // "Compton electron angle w.r.t. z-axis",  //1
                  "electron emission z-coordinate", //2
                  "leftover energy up to 1mm", //3
				  "leftover energy up to 2mm", //4
				  "leftover energy up to 3mm", //5
				  "leftover energy up to 4mm", //6
				  "leftover energy up to 5mm", //7
				  "leftover energy up to 6mm", //8
				  "photon step length (not working right now)", //9 
				  "cherenkov photons per electron",// 10
				  "photon wavelength", //11
				  "initial Compton electron energy", //12 
				  "number of detected photons", //13
				  "number of active channels including DCR and CT", //14
				  "number of detected photons, >1 channels", //15
				  "number of detected photons, >2 channels", //16
				  "number of detected photons, >3 channels", //17
				  "number of detected photons, >4 channels", //18
				  "number of detected photons, >5 channels", //19
                  "num. of det. photons including side arrays WITH DC, >=4 channels", //20
				  "num. of det. photons including side arrays WITH DC, >=6 channels", //21
				  "num. of det. photons including side arrays WITH DC, >=7 channels", //22
				  "num. of det. photons including side arrays WITH DC, >=9 channels", //23
				  "number of detected photons, >14 channels", //24
				  "number of detected photons, >19 channels", //25
				  "num. of det. photons on center array w/o DC, using rnocc", //26
				  "num. of det. photons including side arrays w/o DC using rnocc", //27
				  "TrackID, killed photons", //28
				  "TrackID, detected photons", //29
				  "number of counted photons", //30
				  "gammas reaching PMMA", //31
				  "num. of det. photons including side arrays WITH DC using rnocc", //32
				  "leftover energy after compton", //33	
				  "e- energy", //34		
				  "First Photon Scattering Angle for Accpted Events", // 35
				  "total energy involved in compton process", // 36
				  "scintillator: created photon wavelengths", // 37
				  "Number of created photons in BGO (scint+chkv)"	// 38
};  

  // Default values (to be reset via /analysis/h1/set command)               
	G4int nbins = 10;
	G4double vmin = 0.;
	G4double vmax = 10.;

// get array dimensions for correct display in histograms
	RunControl* RC = RunControl::GetInstance();
	int nbins_x = (RC->arraySideLength_x)*(RC->number_of_arrays_x); 		//Anzahl der Kanaele x
	int nbins_y = (RC->arraySideLength_y)*(RC->number_of_arrays_y);		//Anzahl der Kanaele y
	int nbins_z = (RC->arraySideLength_z)*(RC->number_of_arrays_z);		//Anzahl der Kanaele z

	G4double vmin_x = -0.5*(nbins_x*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_x-1)*(RC->number_of_arrays_x)*(RC->space) + (RC->number_of_arrays_x-1)*RC->spacem);
	G4double vmax_x = +0.5*(nbins_x*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_x-1)*(RC->number_of_arrays_x)*(RC->space) + (RC->number_of_arrays_x-1)*RC->spacem);

	G4double vmin_y = -0.5*(nbins_y*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_y-1)*(RC->number_of_arrays_y)*(RC->space) + (RC->number_of_arrays_y-1)*RC->spacem);
	G4double vmax_y = +0.5*(nbins_y*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_y-1)*(RC->number_of_arrays_y)*(RC->space) + (RC->number_of_arrays_y-1)*RC->spacem);

	// G4double vmin_z = -(nbins_z*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_z-1)*(RC->number_of_arrays_z)*(RC->space) + (RC->number_of_arrays_z-1)*RC->spacem);
	// G4double vmax_z = +(nbins_z*(RC->SiPMchannel_size_xy) + (RC->arraySideLength_z-1)*(RC->number_of_arrays_z)*(RC->space) + (RC->number_of_arrays_z-1)*RC->spacem);

/*		 int n = 40;
		 int v = 40;
 
		 G4int nbins_x = n;
		 G4int nbins_y = v;
		 G4double vmin_x = -n/2.;
		 G4double vmin_y = -v/2.;
		 G4double vmax_x = n/2.;
		 G4double vmax_y = v/2.;*/

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }


////////////////////////////////////////////////
// I didn't use a loop to create 2D histograms. Just copy and paste the last line and choose your titles
// the id will automatically be increased by one, you dont have to set it specifically
// default values are replaced in the main function
////////////////////////////////////////////////


// Create 2D histograms
analysisManager->CreateH2("photonHitMap", "Photon Hit Map", nbins_x, vmin_x, vmax_x, nbins_y, vmin_y, vmax_y);    //0

analysisManager->CreateH2("HitMap_Usable_Only", "Photon Hit Map Usable Events Only", nbins_x, vmin_x, vmax_x, nbins_y, vmin_y, vmax_y);   //1

analysisManager->CreateH2("HitMap_w_SideArrays", "Photon Hit Map with Side Arrays", nbins_x+2*nbins_z, -0.5*(nbins_x+2*nbins_z), +0.5*(nbins_x+2*nbins_z), nbins_y+2*nbins_z, -0.5*(nbins_y+2*nbins_z), +0.5*(nbins_y+2*nbins_z));  //2 
analysisManager->CreateH2("HitMap_w_SideArrays_incl_DC", "Photon Hit Map with Side Arrays including Dark Count", nbins_x+2*nbins_z, -0.5*(nbins_x+2*nbins_z), +0.5*(nbins_x+2*nbins_z), nbins_y+2*nbins_z, -0.5*(nbins_y+2*nbins_z), +0.5*(nbins_y+2*nbins_z));  //3
analysisManager->CreateH2("Photons_vs_Com_Angle", "Detected Photons vs Compton Scattering Angle", 90, 0, 180, 100, 0, 100);  //4
analysisManager->CreateH2("Photons_vs_El_Energy", "Detected Photons vs Electron Energy", 150, 0, 1.500, 100, 0, 100);  //5

// Former used H2 histograms:
// analysisManager->CreateH2("1scatter", "wavelength vs. step length", nbins, vmin, vmax, nbins, vmin, vmax);  //1
// analysisManager->CreateH2("2new HitMap", "Photon Hit Map >4 channels", nbins_x, vmin_x, vmax_x, nbins_y, vmin_y, vmax_y);   //2
// analysisManager->CreateH2("4new HitMap", "Photon Hit Map >4 photons", nbins_x, vmin_x, vmax_x, nbins_y, vmin_y, vmax_y);   //4
// analysisManager->CreateH2("5new HitMap", "Photon Hit Map >14 channels", nbins_x, vmin_x, vmax_x, nbins_y, vmin_y, vmax_y);   //5
// analysisManager->CreateH2("6photohits","Photon Hit Map 0.0mm < z < 0.5mm", nbins, vmin, vmax, nbins, vmin, vmax);  //6  
// analysisManager->CreateH2("7photohits","Photon Hit Map 0.5mm < z < 1.0mm", nbins, vmin, vmax, nbins, vmin, vmax);   //7
// analysisManager->CreateH2("8photohits","Photon Hit Map 1.0mm < z < 1.5mm", nbins, vmin, vmax, nbins, vmin, vmax);   //8
// analysisManager->CreateH2("9photohits","Photon Hit Map 1.5mm < z < 2.0mm", nbins, vmin, vmax, nbins, vmin, vmax);   //9
// analysisManager->CreateH2("10photohits","Photon Hit Map 2.0mm < z < 2.5mm", nbins, vmin, vmax, nbins, vmin, vmax);   //10
// analysisManager->CreateH2("11photohits","Photon Hit Map 2.5mm < z < 3.0mm", nbins, vmin, vmax, nbins, vmin, vmax);   //11
// analysisManager->CreateH2("12photohits","Photon Hit Map 3.0mm < z < 3.5mm", nbins, vmin, vmax, nbins, vmin, vmax);   //12
// analysisManager->CreateH2("13photohits","Photon Hit Map 0.0mm < z < 5.0mm", nbins, vmin, vmax, nbins, vmin, vmax);    //13
// analysisManager->CreateH2("14photohits","Photon Hit Map 5.0mm < z < 8.0mm", nbins, vmin, vmax, nbins, vmin, vmax);   //14

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
