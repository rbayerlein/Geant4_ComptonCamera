
/// \file Analysis.cc
/// \brief Implementation of the Analysis class

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
#include "G4SDManager.hh"
//#include <TTree.h>


using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Analysis* Analysis::fInstance = 0;

Analysis::Analysis()
{
////////////////////////////
// reset default/initial values
////////////////////////////


	lastDist = -100;
	lastTrackID = -1;
	LastTrackID = -1;
	fCerenkovCounter = 0;
	fComptonCounter = 0;	// counts all occurencies of a compton scattering in ALL detector volumes.
	savedElectronStartPosition=false;
  	fEventNumber = -1;
	detpn = 0;
	nofchannels = 0;
	nofchannelstest = 0;
	nofchannelsarray = 0;
	nofchannelstop = 0;
	nofchannelsbottom = 0;
	nofchannelsleft = 0;
	nofchannelsright = 0;
	nofphotons = 0;
	nofphotonstest = 0;
	nofphotonsarray = 0;
	nofphotonstop = 0;
	nofphotonsbottom = 0;
	nofphotonsleft = 0;
	nofphotonsright = 0;

	nofchannelstestnoise = 0;
	nofchannelsarraynoise = 0;
	nofchannelstopnoise = 0;
	nofchannelsbottomnoise = 0;
	nofchannelsleftnoise = 0;
	nofchannelsrightnoise = 0;
	nofphotonstestnoise = 0;
	nofphotonsarraynoise = 0;
	nofphotonstopnoise = 0;
	nofphotonsbottomnoise = 0;
	nofphotonsleftnoise = 0;
	nofphotonsrightnoise = 0;



	RunControl* RC = RunControl::GetInstance();
	
	 n = RC->arraySideLength_x*RC->number_of_arrays_x; 		
	 v = RC->arraySideLength_y*RC->number_of_arrays_y;
	 p = RC->arraySideLength_z*RC->number_of_arrays_z;
	
///////////////////////////
// reset all counters/flags for all the individual channels
///////////////////////////

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHit[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHitarray[i][j] = 0;
		channelHitleft[i][j] = 0;
		channelHitright[i][j] = 0;
		channelHittop[i][j] = 0;
		channelHitbottom[i][j] = 0;
		channelHittest[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHitarraynoise[i][j] = 0;
		channelHitleftnoise[i][j] = 0;
		channelHitrightnoise[i][j] = 0;
		channelHittopnoise[i][j] = 0;
		channelHitbottomnoise[i][j] = 0;
		channelHittestnoise[i][j] = 0;
	}
 }


 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOP[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOParray[i][j] = 0;
		channelNOPleft[i][j] = 0;
		channelNOPright[i][j] = 0;
		channelNOPtop[i][j] = 0;
		channelNOPbottom[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOParraynoise[i][j] = 0;
		channelNOPleftnoise[i][j] = 0;
		channelNOPrightnoise[i][j] = 0;
		channelNOPtopnoise[i][j] = 0;
		channelNOPbottomnoise[i][j] = 0;
	}
 }


/*stringstream ss_data_out;
ss_data_out << RC->RAWfileName << "_data.txt";
text_file = ss_data_out.str();

dat_out.open(text_file.c_str());
dat_out.close();*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Analysis::~Analysis()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





////////////////////////////////////////////////////////////// PREPARE NEW RUN //////////////////////////////////////////////////////////////////////////////////////////////


void Analysis::PrepareNewRun(const G4Run* aRun){

	G4int runN = aRun->GetRunID();
  	G4cout << "Ana:\tRun : " << runN << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
	analysisManager->OpenFile();    
  }else cout << "analysisManager IS NOT ACTIVE! " << endl;

 // define Ntuple for PMTs
  analysisManager->CreateNtuple("PMTHits", "PMTHits");
  analysisManager->CreateNtupleIColumn("PMT_NumberOfDetectedPhotons");
//vector<G4int> vv;
  analysisManager->CreateNtupleSColumn("PMT_Volume_Name");
  analysisManager->CreateNtupleSColumn("PMT_Volume_Name_2");
  analysisManager->CreateNtupleSColumn("PMT_Volume_Name_3");
  analysisManager->CreateNtupleDColumn("PMT_PhotonWavelength", PMT_scint_wavelengths);	// PMT_scint_wavelengths is a vector to be filled in each event
  analysisManager->CreateNtupleIColumn("Photons_In_First_Volume");
  analysisManager->CreateNtupleIColumn("Photons_In_Second_Volume");
  analysisManager->CreateNtupleIColumn("Photons_In_Third_Volume");
 /* analysisManager->CreateNtupleDColumn("Egap");
  analysisManager->CreateNtupleDColumn("Labs");
  analysisManager->CreateNtupleDColumn("Lgap");*/
  analysisManager->FinishNtuple();
//  analysisManager->SetNtupleMerging(true);

// define Ntuple for SiPMs
  analysisManager->CreateNtuple("SiPMHits", "SiPMHits");
  analysisManager->CreateNtupleIColumn("SiPM_NumberOfDetectedPhotons");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronEnergy");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronMomentumDirectionX");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronMomentumDirectionY");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronMomentumDirectionZ");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronVertexX");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronVertexY");
  analysisManager->CreateNtupleDColumn("PrimaryComptonElectronVertexZ");
  analysisManager->CreateNtupleDColumn("ScatteredGammaEnergy");
  analysisManager->FinishNtuple();

/*  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector * SG_Detector = SDManager->FindSensitiveDetector ("GammaDetector", true); // (string name, bool warning)
  SG_Detector->PrepareDetector();*/

//////////////////////////
// read PDE table
//////////////////////////
	wl_PDE.clear();
	PDE.clear();
	ifstream in;
	in.open("Relative_PDE.csv");
	double dummy1, dummy2;
	if(!in) cout << "Relative_PDE.csv: file does not exist" << endl;
	else{
		while(in >> dummy1 >> dummy2){
			wl_PDE.push_back(dummy1);
			PDE.push_back(dummy2);
		}	
///////////////////////////////
// uncomment this to print PDE table
///////////////////////////////

		//cout << "file content:" << endl;
		//for(unsigned int i =0; i<wl_PDE.size(); i++){cout << wl_PDE.at(i) << "\t" << PDE.at(i) << endl;}
	}
	in.close();
}


/////////////////////////////////////////////////////////////////// END OF RUN ACTION /////////////////////////////////////////////////////////////////////////////////////


void Analysis::EndOfRunAction(const G4Run* aRun){
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
	G4int runN = aRun->GetRunID();
  	cout << "Ana:\tEnd of run number: " << runN << G4endl;

////////////////////////////////
//	 save histograms, close file
////////////////////////////////

  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }

}


/////////////////////////////////////////////////////////////////// PREPARE NEW EVENT //////////////////////////////////////////////////////////////////////////////////////

void Analysis::PrepareNewEvent(const G4Event* event){	
/////////////////////////
// plot initial energy distribution of primary particle
/////////////////////////


  in_energy = event -> GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();


/* dat_out.open(text_file.c_str(),ios_base::app | ios_base::out);
 if(!dat_out){cout<< "Ausgabefehler dings" << endl; return;}
 dat_out << event->GetEventID() << ":\t PrepareNewEvent\t primary energy: " ;
 dat_out << in_energy <<endl;
 dat_out.close();*/

 //cout.flush();
 cout << "\r" << "\t" << "eventNumber "<< event->GetEventID() ;

 SG_Scint_Wavelengths.clear();
}


//////////////////////////////////////////////////////////////////// END OF EVENT ACTION ////////////////////////////////////////////////////////////////////////////////////



void Analysis::EndOfEventAction(const G4Event* event){

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	lastTrackID = -1;

//////////////////////////////
// fill histogram h1-10: number of cherenkov photons produced
//////////////////////////////

 analysisManager->FillH1(10,fCerenkovCounter);

//////////////////////////////
// reset cherenkov counter for next particle
//////////////////////////////

 fCerenkovCounter = 0;
 



//now I try to implement the DarkCount

for(int i=0; i<n; i++){
	for(int j=0; j<v; j++){
		G4double rndc = G4UniformRand();
		if(rndc<2.499E-4){
			channelHitarraynoise[i][j]=1;channelNOParraynoise[i][j]++;
			G4double rnct = G4UniformRand();
			if(rnct <0.15){channelHitarraynoise[i][j]=1;channelNOParraynoise[i][j]++;}
		}
	}						
}					

for(int i=0; i<p; i++){
	for(int j=0; j<v; j++){
		G4double rndc = G4UniformRand();
		if(rndc<2.499E-4){
			channelHitrightnoise[i][j]=1;channelNOPrightnoise[i][j]++;
			G4double rnct = G4UniformRand();
			if(rnct <0.15){channelHitrightnoise[i][j]=1;channelNOPrightnoise[i][j]++;}
		}
	}						
}

for(int i=0; i<p; i++){
	for(int j=0; j<v; j++){
		G4double rndc = G4UniformRand();
		if(rndc<2.499E-4){
			channelHitleftnoise[i][j]=1;channelNOPleftnoise[i][j]++;
			G4double rnct = G4UniformRand();
			if(rnct <0.15){channelHitleftnoise[i][j]=1;channelNOPleftnoise[i][j]++;}
		}
	}						
}


for(int i=0; i<n; i++){
	for(int j=0; j<p; j++){
		G4double rndc = G4UniformRand();
		if(rndc<2.499E-4){
			channelHittopnoise[i][j]=1;channelNOPtopnoise[i][j]++;
			G4double rnct = G4UniformRand();
			if(rnct <0.15){channelHittopnoise[i][j]=1;channelNOPtopnoise[i][j]++;}
		}
	}						
}

for(int i=0; i<n; i++){
	for(int j=0; j<p; j++){
		G4double rndc = G4UniformRand();
		if(rndc<2.499E-4){
			channelHitbottomnoise[i][j]=1;channelNOPbottomnoise[i][j]++;
			G4double rnct = G4UniformRand();
			if(rnct <0.15){channelHitbottomnoise[i][j]=1;channelNOPbottomnoise[i][j]++;}
		}
	}						
}

//////////////////////////////
// add accumulative channel and photon numbers
//////////////////////////////


for(int i = 0; i < n; i++){
	for(int j = 0; j < v; j++){
		nofchannels+=channelHit[i][j];
	}
}

for(int i = 0; i < n; i++){
	for(int j = 0; j < v; j++){
		nofchannelsarray+=channelHitarray[i][j];
		nofchannelsleft+=channelHitleft[i][j];
		nofchannelsright+=channelHitright[i][j];
		nofchannelstop+=channelHittop[i][j];
		nofchannelsbottom+=channelHitbottom[i][j];
	}
}

 nofchannelstest=nofchannelsarray+nofchannelsleft+nofchannelsright+nofchannelstop+nofchannelsbottom;

for(int i = 0; i < n; i++){
	for(int j = 0; j < v; j++){
		nofchannelsarraynoise+=channelHitarraynoise[i][j];
		nofchannelsleftnoise+=channelHitleftnoise[i][j];
		nofchannelsrightnoise+=channelHitrightnoise[i][j];
		nofchannelstopnoise+=channelHittopnoise[i][j];
		nofchannelsbottomnoise+=channelHitbottomnoise[i][j];
	}
}

 nofchannelstestnoise=nofchannelsarraynoise+nofchannelsleftnoise+nofchannelsrightnoise+nofchannelstopnoise+nofchannelsbottomnoise;

nofchannels_sig_and_noise=nofchannelstestnoise+nofchannelstest;

for(int i = 0; i < n; i++){
	for(int j =0;j<v;j++){
		nofphotons += channelNOP[i][j];
	}
}

for(int i = 0; i < n; i++){
	for(int j = 0; j < v; j++){
		nofphotonsarray += channelNOParray[i][j];
		nofphotonsleft += channelNOPleft[i][j];
		nofphotonsright += channelNOPright[i][j];
		nofphotonstop += channelNOPtop[i][j];
		nofphotonsbottom += channelNOPbottom[i][j];
	}
}

 nofphotonstest = nofphotonstop+nofphotonsbottom+nofphotonsright+nofphotonsleft+nofphotonsarray;


for(int i = 0; i < n; i++){
	for(int j = 0; j < v; j++){
		nofphotonsarraynoise += channelNOParraynoise[i][j];
		nofphotonsleftnoise += channelNOPleftnoise[i][j];
		nofphotonsrightnoise += channelNOPrightnoise[i][j];
		nofphotonstopnoise += channelNOPtopnoise[i][j];
		nofphotonsbottomnoise += channelNOPbottomnoise[i][j];
	}
}

 nofphotonstestnoise = nofphotonstopnoise+nofphotonsbottomnoise+nofphotonsrightnoise+nofphotonsleftnoise+nofphotonsarraynoise;

 nofphotons_sig_and_noise = nofphotonstestnoise + nofphotonstest;

//G4cout << "nofphotons/channels" << "\t" << nofphotons << "\t" << nofchannels << endl;

////////////////////////
// here you can set the conditions under which the total number of detected photons is filled into histogram h1-26
////////////////////////

 RunControl* RC = RunControl::GetInstance(); 
 int primaryphotonscounter = RC->primaryphotonscounter; 
 analysisManager->FillH1(30,primaryphotonscounter);

if(nofchannelsarray >= RC->rnocc){
	//if(nofphotonsarray >= RC->rnocc){
 analysisManager->FillH1(26,nofphotonsarray);
	//}
}

if(nofchannelstest >= RC->rnocc){
	//if(nofphotonstest >= RC->rnocc){
 analysisManager->FillH1(27,nofphotonstest);
	//}
}

if(nofchannels_sig_and_noise >= RC->rnocc){
	//if(nofphotons_sig_and_noise >= RC->rnocc){
 analysisManager->FillH1(32,nofphotons_sig_and_noise);
	//}
}

////////////////////////
// here you can set the conditions under which the hit map h2-3 is filled. there is a for-loop for each detector cell.
// you can copy-and-paste them if you want to examine multiple conditions, you have to change the histogram id. ids 2 to 5 are reserved for different conditions
// make sure you adjust the histogram title in the histo-manager
////////////////////////

// the first two only for primary photons
//if(fPhotoCounter == 0){
//	if(fComptonCounter == 1){


		if(nofchannelstest > 2){
			if(nofphotonstest > 4){
				analysisManager->FillH1(13,nofphotonsarray);
				for(int i=0; i<n;i++){
					for(int j=0; j<v;j++){
						for(int l=0; l<channelNOParray[i][j];l++){analysisManager->FillH2(1,-n/2.+1/2+i,v/2.-1-j);}
					}							
				}
			}
		}
//	}
//}
					

//if(fPhotoCounter == 0){
//	if(fComptonCounter == 1){


		if(nofchannelstest > 2){
			if(nofphotonstest > 4){
				for(int i=0; i<n;i++){
					for(int j=0; j<v;j++){
						//for(int k=0; k< z;++k){
							for(int l=0; l<channelNOParray[i][j];l++){analysisManager->FillH2(2,-n/2.+1/2+i,v/2.-1-j);}
							for(int l=0; l<channelNOPleft[i][j];l++){analysisManager->FillH2(2,-n/2.-1+1/2-i,v/2.-1-j);}
							for(int l=0; l<channelNOPright[i][j];l++){analysisManager->FillH2(2,n/2.-1/2+i,v/2.-1-j);}
							for(int l=0; l<channelNOPtop[i][j];l++){analysisManager->FillH2(2,-n/2.+1/2+i,v/2.-1/2+j);}
							for(int l=0; l<channelNOPbottom[i][j];l++){analysisManager->FillH2(2,-n/2.+1/2+i,-v/2.-1+1/2-j);}
						//}
					}							
				}
			}
		}
//	}
//}

//if(fPhotoCounter == 0){
//	if(fComptonCounter == 1){


		if(nofchannels_sig_and_noise > 2){
			if(nofphotons_sig_and_noise > 4){
				for(int i=0; i<n;i++){
					for(int j=0; j<v;j++){
							for(int l=0; l<channelNOParraynoise[i][j]+channelNOParray[i][j];l++){analysisManager->FillH2(3,-n/2.+1/2+i,v/2.-1-j);}
							for(int l=0; l<channelNOPleftnoise[i][j]+channelNOPleft[i][j];l++){analysisManager->FillH2(3,-n/2.-1+1/2-i,v/2.-1-j);}
							for(int l=0; l<channelNOPrightnoise[i][j]+channelNOPright[i][j];l++){analysisManager->FillH2(3,n/2.-1/2+i,v/2.-1-j);}
							for(int l=0; l<channelNOPtopnoise[i][j]+channelNOPtop[i][j];l++){analysisManager->FillH2(3,-n/2.+1/2+i,v/2.-1/2+j);}
							for(int l=0; l<channelNOPbottomnoise[i][j]+channelNOPbottom[i][j];l++){analysisManager->FillH2(3,-n/2.+1/2+i,-v/2.-1+1/2-j);}
					}							
				}
			}
		}
//	}
//}

//here the new histogramm
		
if(fComptonCounter>0){
//if(RC->theta_compt>-50){
	if(nofchannelstest >= RC->rnocc){
			analysisManager->FillH2(4,RC->theta_compt,nofphotonstest);
			analysisManager->FillH2(5,electronEnergy, nofphotonstest);
			analysisManager->FillH1(35,RC->theta_compt);
			analysisManager->FillNtupleIColumn(1, 0, nofchannelstest);
			analysisManager->FillNtupleDColumn(1, 1, electronEnergy);
			double em_x = electronMomentumDirection.x();
			double em_y = electronMomentumDirection.y();
			double em_z = electronMomentumDirection.z();
			analysisManager->FillNtupleDColumn(1,2, em_x);
			analysisManager->FillNtupleDColumn(1,3, em_y);
			analysisManager->FillNtupleDColumn(1,4, em_z);
			double ev_x = electronVertex.x();
			double ev_y = electronVertex.y();
			double ev_z = electronVertex.z();
			analysisManager->FillNtupleDColumn(1,5, ev_x);
			analysisManager->FillNtupleDColumn(1,6, ev_y);
			analysisManager->FillNtupleDColumn(1,7, ev_z);
			analysisManager->FillNtupleDColumn(1,8, photonenergy);
			
		//	if(PMT_scint_wavelengths.size()>0){
				analysisManager->AddNtupleRow(0);
				analysisManager->AddNtupleRow(1);
		//	}
	}PMT_scint_wavelengths.clear();
//}
}

if(fComptonCounter>0){analysisManager->FillH1(1,fComptonCounter);}


//////////////////////////////////
// uncomment this section to reset photo/compton counter after each event
// you can use this to distinguish primary photon events by electron production mechanism
/////////////////////////////////

// fPhotoCounter = 0;
 fComptonCounter=0;
 savedElectronStartPosition=false;

/////////////////////////////////
// uncomment this section to write a photon number output file in .txt format
/////////////////////////////////

/* dat_out.open(text_file.c_str(),ios_base::app | ios_base::out);
 if(!dat_out){cout<< "Ausgabefehler dings" << endl; return;}
 dat_out << event->GetEventID() << ":\tend of event Action\n===========" << endl;
 dat_out.close();*/


/////////////////////////////
// here, detected photon number distributions are filled for different minimum channel number conditions
// most of these are obsolete, some are commented out
// standard condition: nofchannels > 4
/////////////////////////////


    if(nofchannels_sig_and_noise>0){
  	analysisManager->FillH1(14,nofchannels_sig_and_noise);
  	}
	if(nofchannelstest>1){
	analysisManager->FillH1(15,nofphotonstest);
	}
	if(nofchannelstest>2){
	analysisManager->FillH1(16,nofphotonstest);
	}
	if(nofchannelstest>3){
	analysisManager->FillH1(17,nofphotonstest);
	}
	if(nofchannelstest>4){
	analysisManager->FillH1(18,nofphotonstest);
	}
	if(nofchannelstest>5){
	analysisManager->FillH1(19,nofphotonstest);
	}

	if(nofchannels_sig_and_noise >=4){
	analysisManager->FillH1(20,nofphotons_sig_and_noise);
	}
	if(nofchannels_sig_and_noise>=6){
	analysisManager->FillH1(21,nofphotons_sig_and_noise);
	}
	if(nofchannels_sig_and_noise>=7){
	analysisManager->FillH1(22,nofphotons_sig_and_noise);
	}
	if(nofchannels_sig_and_noise>=9){
	analysisManager->FillH1(23,nofphotons_sig_and_noise);
	}
	if(nofchannelstest>14){
	analysisManager->FillH1(24,nofphotonstest);
	}
	if(nofchannelstest>19){
	analysisManager->FillH1(25,nofphotonstest);
	}


///////////////////////////////
// reset all counters and flags
///////////////////////////////


 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHit[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHitarray[i][j] = 0;
		channelHittop[i][j] = 0;
		channelHitbottom[i][j] = 0;
		channelHitleft[i][j] = 0;
		channelHitright[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelHitarraynoise[i][j] = 0;
		channelHittopnoise[i][j] = 0;
		channelHitbottomnoise[i][j] = 0;
		channelHitleftnoise[i][j] = 0;
		channelHitrightnoise[i][j] = 0;
	}
 }


 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOP[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOParray[i][j] = 0;
		channelNOPleft[i][j] = 0;
		channelNOPright[i][j] = 0;
		channelNOPtop[i][j] = 0;
		channelNOPbottom[i][j] = 0;
	}
 }

 for(int i=0; i < n;i++){
 	for(int j=0; j< v;j++){
		channelNOParraynoise[i][j] = 0;
		channelNOPleftnoise[i][j] = 0;
		channelNOPrightnoise[i][j] = 0;
		channelNOPtopnoise[i][j] = 0;
		channelNOPbottomnoise[i][j] = 0;
	}
 }

nofchannels = 0;
nofphotons = 0;

nofchannelstest = 0;
nofphotonstest = 0;

nofchannelsarray = 0;
nofphotonsarray = 0;

nofchannelstop = 0;
nofphotonstop = 0;

nofchannelsbottom = 0;
nofphotonsbottom = 0;

nofchannelsleft = 0;
nofphotonsleft = 0;

nofchannelsright = 0;
nofphotonsright = 0;


nofchannelstestnoise = 0;
nofphotonstestnoise = 0;

nofchannelsarraynoise = 0;
nofphotonsarraynoise = 0;

nofchannelstopnoise = 0;
nofphotonstopnoise = 0;

nofchannelsbottomnoise = 0;
nofphotonsbottomnoise = 0;

nofchannelsleftnoise = 0;
nofphotonsleftnoise = 0;

nofchannelsrightnoise = 0;
nofphotonsrightnoise = 0;

	for (int i = 0; i < SG_Scint_Wavelengths.size(); ++i)
	{
		//cout << SG_Scint_Wavelengths.at(i) << endl;
		analysisManager->FillH1(37, SG_Scint_Wavelengths.at(i));
	}
	analysisManager->FillH1(38, SG_Scint_Wavelengths.size());

/*	if(PMT_scint_wavelengths.size()>0){
		analysisManager->AddNtupleRow(0);
		analysisManager->AddNtupleRow(1);
	}
	PMT_scint_wavelengths.clear();*/
}	
	
//////////////////////////////////////////////////////////////////////////////////////// Pre User Tracking Action ///////////////////////////////////////////////////////////////

void Analysis::PreUserTrackingAction(const G4Track* /*track*/){

}

//////////////////////////////////////////////////////////////////////////////////////// Post User tracking Action ////////////////////////////////////////////////////////////////

void Analysis::PostUserTrackingAction(const G4Track* /*track*/){

}

//////////////////////////////////////////////////////////////////////////////////////// Stepping Action //////////////////////////////////////////////////////////////////////////

void Analysis::SteppingAction(const G4Step* step){
 RunControl* RC = RunControl::GetInstance();

	const G4StepPoint* PSP = step->GetPreStepPoint();
	double z = PSP->GetPosition().z();
	// double x = PSP->GetPosition().x();
	// double y = PSP->GetPosition().y();

	if(z > RC->PMMA_thick/2+RC->SiPM_thick+RC->ref_thick+0.0001) return;	// skip this step if it happens behind PMMA as it has nothing to do with Chkv then.
	// if(x > RC->PMMA_thick/2+0.0001) return;
	// if(y > RC->PMMA_thick/2+RC->SiPM_thick+RC->ref_thick+0.0001) return;

///////////////////////
// instantiate analysis manager
///////////////////////

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();	

///////////////////////
// get particle name
///////////////////////


  G4Track* track = step->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();


/* dat_out.open(text_file.c_str(),ios_base::app | ios_base::out);
 if(!dat_out){cout<< "Ausgabefehler dings in stepping action" << endl; return;}
 if(ParticleName == "gamma"){
 	dat_out << "TrackID: " << step->GetTrack()->GetTrackID() << "\t(gamma)" << endl;
 }

*/

///////////////////////////
// get volume name
///////////////////////////

//G4VPhysicalVolume* lvolume 
//    = step->GetPreStepPoint()->GetTouchableHandle()
//      ->GetVolume();
 G4String volname2 = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetName();					
// G4cout << "volume name:\t" << volname2 << "\tparticle name:\t" << ParticleName << endl;


///////////////////////////
// Save gamma z position and count if inside PMMA
///////////////////////////

 G4double PMMA_thick = RC->PMMA_thick;
 int gammahit = 0;
 G4String volname3 = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetName();		
 if(ParticleName == "gamma"){
 	const G4StepPoint* gammaprep = step->GetPreStepPoint();
 	G4ThreeVector gammapos = gammaprep->GetPosition(); 
 	G4double gammaposz = gammapos.z();
 	// dat_out << gammaposz << endl;
	if(gammaposz>-(PMMA_thick/2+0.000001)*mm && gammaposz< -(PMMA_thick/2-0.000001)*mm && volname3 == "PMMA"){
		gammahit++;
	//	dat_out << "counted in histo 31 as gammahit" << endl;
	}
 }
 analysisManager->FillH1(31,gammahit);


char name[1000];
char ref_name[1000];
char name_refright[1000];
char name_refleft[1000];
char name_refup[1000];
char name_refdown[1000];






if(ParticleName == "opticalphoton"){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < v;j++){
			sprintf(name,"%d_%d",i,j);
			if(volname2 == name){
				G4double ID = step->GetTrack()->GetTrackID();							//to fill H1-28
				track->SetTrackStatus(fStopAndKill);	
			//	G4cout << "volume match\t" << volname2.c_str() << "\t" << name << endl;
				analysisManager->FillH1(28,ID);									
			}
		}
 	}


/////////////////////////////////////////////////
// here, the photon detection efficiency (PDE) is implemented. the wavelength dependent PDE is extracted from file
// and then compared to a random number
// photon wavelength is filled into h1-11
/////////////////////////////////////////////////




	G4double rn = G4UniformRand();
	const G4StepPoint* phendPoint = step->GetPreStepPoint();
	G4double length = step->GetStepLength();
	G4ThreeVector phposition = phendPoint->GetPosition();
	G4double phenergy = phendPoint->GetKineticEnergy();
	G4double phlambda = 1.2398E-6*eV*m/phenergy;// wavelength [mm]
	//G4cout << phlambda << endl;
	//analysisManager->FillH1(11,phlambda);
	G4double xcoord = phposition.x();
	G4double ycoord = phposition.y();
	G4double zcoord2 = phposition.z();
	int wlposition = wl_PDE.size()-1;
	for (unsigned int i = 0 ; i < wl_PDE.size(); ++i){
		if (wl_PDE.at(i) > phlambda*1E6){
			wlposition = i;
				break;
		}	
	}

			
	if (ParticleName == "opticalphoton"){
		const G4StepPoint* phot_prepoint = step->GetPreStepPoint();
		G4ThreeVector phot_position = phot_prepoint->GetPosition();
	    	G4double phot_prez = phot_position.z();
			//if(phot_prez < 4.501*mm){
		if(rn < PDE.at(wlposition)*0.582){	//0.519 @ V=3V
											//0.552 @ V=3.5V
											//0.582 @ V=4V
											//0.614 @ V=4.5V
											//0.631 @ V=5V
			// analysisManager->FillH1(11,zcoord2);
			analysisManager->FillH2(0,xcoord,ycoord);
			G4double rn2 = G4UniformRand();


			const G4StepPoint* phot_postpoint = step->GetPostStepPoint();
			G4ThreeVector phot_post = phot_postpoint->GetPosition();
			G4double phot_postx = phot_post.x();
			G4double phot_posty = phot_post.y();
			G4double phot_postz = phot_post.z();

			G4ThreeVector momentum = phot_prepoint->GetMomentumDirection();    	//test impuls des teilchens auf der oberflaeche
			G4double momentumz = momentum.z();					//test definiere z-komponente
			G4double momentumx = momentum.x();
			G4double momentumy = momentum.y();

			G4String volname = step->GetPreStepPoint()->GetTouchableHandle()
					      ->GetVolume()->GetName();
			

/////////////////////////////////////
// Here h2-6 to h2-14 are filled with the photon x and y coordinates at the detector plane sorted by production z coordinate
// These hit maps don't take dead space between detectors into account and also rejection conditions cannot be applied
// I haven't thought of a nice way to implement this yet and maybe I never will.
// You can adjust the prez-intervals in the respective if-clauses. Make sure to adjust the title in the HistoManager if you change anything.
/////////////////////////////////////
/*
			if(3.999*mm < phot_postz && phot_postz < 4.001*mm){
				if(-4.0*mm < phot_prez && phot_prez < -3.5*mm){analysisManager->FillH2(6,phot_postx,phot_posty);}
				if(-3.5*mm < phot_prez && phot_prez < -3.0*mm){analysisManager->FillH2(7,phot_postx,phot_posty);}
				if(-3.0*mm < phot_prez && phot_prez < -2.5*mm){analysisManager->FillH2(8,phot_postx,phot_posty);}
				if(-2.5*mm < phot_prez && phot_prez < -2.0*mm){analysisManager->FillH2(9,phot_postx,phot_posty);}
				if(-2.0*mm < phot_prez && phot_prez < -1.5*mm){analysisManager->FillH2(10,phot_postx,phot_posty);}
				if(-1.5*mm < phot_prez && phot_prez < -1.0*mm){analysisManager->FillH2(11,phot_postx,phot_posty);}
				if(-1.0*mm < phot_prez && phot_prez < -0.5*mm){analysisManager->FillH2(12,phot_postx,phot_posty);}
				if(-4.0*mm < phot_prez && phot_prez < 1.0*mm){analysisManager->FillH2(13,phot_postx,phot_posty);}
				if(1.0*mm < phot_prez && phot_prez < 4.0*mm){analysisManager->FillH2(14,phot_postx,phot_posty);}
			}*/

//////////////////////////////////////
// final photon volume is called.
// if a detector volume is hit the respective flag is set and the photon counter is incremented by 1
//////////////////////////////////////



			for(int i=0; i<n; i++){
				for(int j=0; j<v; j++){
					sprintf(ref_name,"%s_%d_%d","ref",i,j);
					if(volname == ref_name && momentumz > 0){				//momentumz>0 to count only photons going into the SiPM (not the reflecting 																	ones)
						//if(step->GetTrack()->GetTrackID()!= LastTrackID){		//(**)
						channelHit[i][j]=1;channelNOP[i][j]++;
						//LastTrackID=step->GetTrack()->GetTrackID();}     		//(**)
						}
					}
			}

			for(int i=0; i<n; i++){
				for(int j=0; j<v; j++){
					sprintf(ref_name,"%s_%d_%d","ref",i,j);
					sprintf(name_refright,"%s_%d_%d","refright",n+i,j);
					sprintf(name_refleft,"%s_%d_%d","refleft",-i-1,j);
					sprintf(name_refup,"%s_%d_%d","refup",i,-j-1);
					sprintf(name_refdown,"%s_%d_%d","refdown",i,v+j);

					if(volname == ref_name && momentumz > 0){
						if(step->GetTrack()->GetTrackID()!= LastTrackID){			//uncomment this line..
						channelHitarray[i][j]=1;channelNOParray[i][j]++; 
						LastTrackID=step->GetTrack()->GetTrackID();}				//..and this line to make sure that each photon is only counted once
					}										//but you have to comment the two lines above (**)

					if(volname == name_refright && momentumx < 0){						
						if(step->GetTrack()->GetTrackID()!= LastTrackID){
						channelHitright[i][j]=1;channelNOPright[i][j]++; 
						LastTrackID=step->GetTrack()->GetTrackID();}
					}

					if(volname == name_refleft && momentumx > 0){
						if(step->GetTrack()->GetTrackID()!= LastTrackID){
						channelHitleft[i][j]=1;channelNOPleft[i][j]++; 
						LastTrackID=step->GetTrack()->GetTrackID();}
					}

					if(volname == name_refup && momentumy > 0){
						if(step->GetTrack()->GetTrackID()!= LastTrackID){
						channelHittop[i][j]=1;channelNOPtop[i][j]++; 
						LastTrackID=step->GetTrack()->GetTrackID();}
					}

					if(volname == name_refdown && momentumy < 0){
						if(step->GetTrack()->GetTrackID()!= LastTrackID){
						channelHitbottom[i][j]=1;channelNOPbottom[i][j]++; 
						LastTrackID=step->GetTrack()->GetTrackID();}
					}
				}
			}
					

//////////////////////////////////////////
// Crosstalk: there is a 15% probability that a photon is counted twice
// this is implemented by simply repeating the above steps in 15% of the cases here:
//////////////////////////////////////////


			if(rn2 < 0.15){
				analysisManager->FillH2(0,xcoord,ycoord);

/*				if(3.999*mm < phot_postz && phot_postz < 4.001*mm){
					if(-4.0*mm < phot_prez && phot_prez < -3.5*mm){analysisManager->FillH2(6,phot_postx,phot_posty);}
					if(-3.5*mm < phot_prez && phot_prez < -3.0*mm){analysisManager->FillH2(7,phot_postx,phot_posty);}
					if(-3.0*mm < phot_prez && phot_prez < -2.5*mm){analysisManager->FillH2(8,phot_postx,phot_posty);}
					if(-2.5*mm < phot_prez && phot_prez < -2.0*mm){analysisManager->FillH2(9,phot_postx,phot_posty);}
					if(-2.0*mm < phot_prez && phot_prez < -1.5*mm){analysisManager->FillH2(10,phot_postx,phot_posty);}
					if(-1.5*mm < phot_prez && phot_prez < -1.0*mm){analysisManager->FillH2(11,phot_postx,phot_posty);}
					if(-1.0*mm < phot_prez && phot_prez < -0.5*mm){analysisManager->FillH2(12,phot_postx,phot_posty);}
					if(-4.0*mm < phot_prez && phot_prez < 1.0*mm){analysisManager->FillH2(13,phot_postx,phot_posty);}
					if(1.0*mm < phot_prez && phot_prez < 4.0*mm){analysisManager->FillH2(14,phot_postx,phot_posty);}
				}*/



				for(int i=0; i<n; i++){
					for(int j=0; j<v; j++){
						sprintf(ref_name,"%s_%d_%d","ref",i,j);
						if(volname == ref_name && momentumz > 0){
							//if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHit[i][j]=1;channelNOP[i][j]++;
							//LastTrackID=step->GetTrack()->GetTrackID();}     			//test: momentumz>0 heisst photon wird nicht reflektiert
							}
						}
				}

				for(int i=0; i<n; i++){
					for(int j=0; j<v; j++){
						sprintf(ref_name,"%s_%d_%d","ref",i,j);
						sprintf(name_refright,"%s_%d_%d","refright",n+i,j);
						sprintf(name_refleft,"%s_%d_%d","refleft",-i-1,j);
						sprintf(name_refup,"%s_%d_%d","refup",i,-j-1);
						sprintf(name_refdown,"%s_%d_%d","refdown",i,v+j);

						if(volname == ref_name && momentumz > 0){
							if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHitarray[i][j]=1;channelNOParray[i][j]++; 
							LastTrackID=step->GetTrack()->GetTrackID();}
						}

						if(volname == name_refright && momentumx < 0){
							if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHitright[i][j]=1;channelNOPright[i][j]++; 
							LastTrackID=step->GetTrack()->GetTrackID();}
						}

						if(volname == name_refleft && momentumx > 0){
							if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHitleft[i][j]=1;channelNOPleft[i][j]++; 
							LastTrackID=step->GetTrack()->GetTrackID();}
						}

						if(volname == name_refup && momentumy > 0){
							if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHittop[i][j]=1;channelNOPtop[i][j]++; 
							LastTrackID=step->GetTrack()->GetTrackID();}
						}

						if(volname == name_refdown && momentumy < 0){
							if(step->GetTrack()->GetTrackID()!= LastTrackID){
							channelHitbottom[i][j]=1;channelNOPbottom[i][j]++; 
							LastTrackID=step->GetTrack()->GetTrackID();}
						}
					}
				}


			}// if(rn2 < 0.15) CHECK OCT

		}// if (rn < PDE.at(wlposition)*0.582) CHECK PDE
				
	}// if(ParticleName == "opticalphoton") CHECK OPICAL PHOTON


////////////////////////////////////////////////////////
// uncomment this if you want to fill histograms with photon step lengths and photon energy vs. step length
// this was done to check if photon absorption works correctly
// to get reasonable results, the PMMA volume has to be increased significantly
// everything looks good, you usually won't need this anymore
/////////////////////////////////////////////////////// 



	//		analysisManager->FillH2(1,1.2398E-6*eV*m/phenergy,length);


}	//test, sonst wieder entfernen aber vermutlich Klammer zur if(opticalphoton)			



//////////////////////////////////////////////
// uncomment this if you want to fill the photon angle w.r.t. the z-axis into h1-1
// to get reasonable results, you should disable reflection at PMMA boundaries (remove the assignment of a refractive index to
// the envelope material in the detector construction
// don't forget to reverse that step when you're done!
//////////////////////////////////////////////



//if (ParticleName == "opticalphoton"){
//	const G4StepPoint* PhotoPoint = step->GetPostStepPoint();
//	G4ThreeVector photodir = PhotoPoint->GetMomentumDirection();
//   	G4double costheta = photodir.z();
//	G4double theta = acos (costheta) * 180.0 / 3.14159;
//    	analysisManager->FillH1(1,theta);
//}




      
///////////////////////////////
// here, the cherenkov photon counter and detected photon counter are is reset after all secondaries are tracked and a new event is ready to start
// this could probably also be done in the "prepare new event" section, but this works, so whatever
///////////////////////////////


G4int eventNumber = G4RunManager::GetRunManager()->
                                              GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
     fEventNumber = eventNumber;
     fCerenkovCounter = 0;
     detpn = 0;
	

  }

//////////////////////////////////////////
// this section was used to plot the emission angle of photon induced (Compton-)electrons
// to only plot the initial momentum direction, "kill" the particle track after the first step
// take care to reverse that step when you're done investigating the emission angle
//////////////////////////////////////////

	if(ParticleName == "e-"){
		const G4StepPoint* prepoint = step->GetPreStepPoint();

		G4ThreeVector photodir = prepoint->GetMomentumDirection();
	  	G4double costheta = photodir.z();
		G4double theta = acos (costheta) * 180.0 / 3.14159;			//das hier auch fuer PGA, um fehler mit e- zu finden
	    //	analysisManager->FillH1(1,theta);

		G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
		// if(volume->GetName() == "PMMA"){
		if (step->GetTrack()->GetCreatorProcess()->GetProcessName() == "compt"){
			G4ThreeVector compt_pre = prepoint->GetPosition();
			G4double compt_prez = compt_pre.z();
			if(!savedElectronStartPosition){
				analysisManager->FillH1(2,compt_prez);
				savedElectronStartPosition=true;
			}
		}
		////////////////////////////
// uncomment this to stop electron tracks after first step to only fill initial energy and momentum direction
////////////////////////////
		//track->SetTrackStatus(fStopAndKill);
	}

///////////////////////////////////
// here, secondary particles are are counted: Cherenkov photons, and, in case of photons as primary particles,
// Compton and photo electrons
///////////////////////////////////


  const std::vector<const G4Track*>* secondaries =
                                            step->GetSecondaryInCurrentStep();



if (secondaries->size()>0) {
	G4double totalComptonEnergy=0;
    for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {

        	if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
           		if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov") fCerenkovCounter++; 
			}	
        	
			if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4Electron::ElectronDefinition()){             
				if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "compt"){
		 			fComptonCounter++;
		 			// dat_out << "COMPTON SCATTER - track ID " << secondaries->at(i)->GetTrackID() << endl;
		 			// dat_out << "electron energy: " << secondaries->at(i)->GetKineticEnergy() << "\tgamma energy: " << step->GetTrack()->GetKineticEnergy() << 
		 			// "sum: " << secondaries->at(i)->GetKineticEnergy()+step->GetTrack()->GetKineticEnergy() << endl;
		 			totalComptonEnergy+=secondaries->at(i)->GetKineticEnergy();
		 			if(fComptonCounter==1){
				    	totalComptonEnergy+=step->GetTrack()->GetKineticEnergy();
						analysisManager->FillH1(36,totalComptonEnergy);
						analysisManager->FillH1(33,step->GetTrack()->GetKineticEnergy());

						electronEnergy = secondaries->at(i)->GetKineticEnergy();
						electronMomentumDirection = secondaries->at(i)->GetMomentumDirection();
						electronVertex = secondaries->at(i)->GetPosition();

						G4ThreeVector photodira = RC->startmomalpha;
						G4double ax = photodira.x();
						G4double ay = photodira.y();
						G4double az = photodira.z();
						const G4StepPoint* postpoint = step->GetPostStepPoint();
						G4ThreeVector photodirb = postpoint->GetMomentumDirection();
						G4double bx = photodirb.x();
						G4double by = photodirb.y();
						G4double bz = photodirb.z();
						G4double costheta = ((ax*bx)+(ay*by)+(az*bz))/(sqrt((ax*ax)+(ay*ay)+(az*az))*sqrt((bx*bx)+(by*by)+(bz*bz)));
						if(costheta<-1){RC->theta_compt=180;}
						if(costheta>1){RC->theta_compt=0;}
						if(costheta>=-1&&costheta<=1){RC->theta_compt=acos(costheta)*180.0 / 3.14159;}
						// G4cout << "THETA" << RC->theta_compt<<endl;
						photonenergy = postpoint->GetKineticEnergy();
						// track->SetTrackStatus(fStopAndKill); //to only plot the first compton, make sure to uncomment again

					}
				}
			}

			if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4Electron::ElectronDefinition()){             
				if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "phot") fPhotoCounter++; 
			}
        }
    }

}

//test
if(fComptonCounter==1){
	if(ParticleName=="gamma"){
		if (secondaries->size()>0) {
			for(unsigned int i=0; i<secondaries->size(); ++i) {
				if (secondaries->at(i)->GetParentID()>0) {
					if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4Electron::ElectronDefinition()){
						if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "phot"){
							int k = 1; analysisManager->FillH1(34,k);
						} 
					}  
				}
			}
		}

	}

}
 // dat_out.close();
}//}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

