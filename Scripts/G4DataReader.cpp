#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdint.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <signal.h>
#include <dirent.h>
#include <sys/stat.h>
#include <stdio.h>

#include "TStyle.h"
#include "TFile.h"
#include "TAttText.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TF1.h"
#include "TGraph.h"
#include "TColor.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TMath.h"
#include "THStack.h"
#include "Rtypes.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TMarker.h"

using namespace std;

string inputFileName = "/home/bayerlein/Documents/Analysis/Analysis_and_Results/20200323_Simulation_Analysis/Scattered_Gamma_Detector/TEST_02.root";
string VolumesFileName = "VolumeList.txt";
int tPMT_NumPhotons, tSiPM_NumPhotons;
Char_t PMT_Volume_Name[16];
Char_t PMT_Volume_Name_2[16];
Char_t PMT_Volume_Name_3[16];

string volumeNames[1024];
double PMT_pos_x[1024];
double PMT_pos_y[1024];
double PMT_pos_z[1024];
int num_PMT_Volumes;

double pCee, pCevX, pCevY, pCevZ; //primary Compton electron energy in Chkv detector

int GetContinuousVolumeNumber(char *c);
void fillVectors(string input);
void G4DataReader(){
	cout << "defining output file and graphs" << endl;
	TH1F *h_cos_scat_angle= new TH1F("cos_scat_angle", "cos_scat_angle", 40, -1, 1);
	TH1F *h_scat_angle= new TH1F("scat_angle", "scat_angle", 90, 0, 180);

	cout << "starting to read from file " << VolumesFileName << endl;
	fillVectors(VolumesFileName);

	cout << "starting to read from file " << inputFileName << endl;
	TFile *f = new TFile(inputFileName.c_str(), "READ");

// read from PMT tree
	TTree *PMT_Tree = (TTree*)f->Get("PMTHits");
	if (PMT_Tree==NULL){cout << "No PMT tree found." << endl; return;}

/*	PMT_Tree->SetBranchAddress("Channel",&tChannel);
	PMT_Tree->SetBranchAddress("ToT", &tEnergy);
	PMT_Tree->SetBranchAddress("Time", &tTime);
	PMT_Tree->SetBranchAddress("numChannels", &tnumCh);*/
	PMT_Tree->SetBranchAddress("PMT_NumberOfDetectedPhotons", &tPMT_NumPhotons);
	PMT_Tree->SetBranchAddress("PMT_Volume_Name", &PMT_Volume_Name);
	PMT_Tree->SetBranchAddress("PMT_Volume_Name_2", &PMT_Volume_Name_2);
	PMT_Tree->SetBranchAddress("PMT_Volume_Name_3", &PMT_Volume_Name_3);


// read from SiPM tree
	TTree *SiPM_Tree = (TTree*)f->Get("SiPMHits");
	if (SiPM_Tree==NULL){cout << "No SiPM tree found." << endl; return;}
	SiPM_Tree->SetBranchAddress("SiPM_NumberOfDetectedPhotons", &tSiPM_NumPhotons);
	SiPM_Tree->SetBranchAddress("PrimaryComptonElectronEnergy", &pCee);
	SiPM_Tree->SetBranchAddress("PrimaryComptonElectronVertexX", &pCevX);
	SiPM_Tree->SetBranchAddress("PrimaryComptonElectronVertexY", &pCevY);
	SiPM_Tree->SetBranchAddress("PrimaryComptonElectronVertexZ", &pCevZ);

	int n_entries=PMT_Tree->GetEntries();
	cout << n_entries << " entries in the PMT_Tree and "<< SiPM_Tree->GetEntries()  << " in the SiPM tree." << endl;
	int it=0;
	while(it < n_entries){
		SiPM_Tree->GetEntry(it);
		PMT_Tree->GetEntry(it);
		char test[10]="P";
		if(PMT_Volume_Name[0] != test[0]) {cout << "event " << it << ": no photon hit in PMT" << endl; it++; continue;}
		int CVN = GetContinuousVolumeNumber(PMT_Volume_Name);
		double px = PMT_pos_x[CVN];
		double py = PMT_pos_y[CVN];
		double pz = PMT_pos_z[CVN];
	//	cout << "Position:\n" 	<< px << "\t" << py << "\t"	<< pz << endl;
		double e_gamma_vector[3];
		e_gamma_vector[0] = px - pCevX;
		e_gamma_vector[1] = py - pCevY;
		e_gamma_vector[2] = pz - pCevZ;
		double cos_scat_angle = e_gamma_vector[2]/sqrt(pow(e_gamma_vector[0],2) + pow(e_gamma_vector[1],2) + pow(e_gamma_vector[2],2));
		cout << "cos_scat_angle " << TMath::ACos(cos_scat_angle)*180/TMath::Pi() << endl;
		h_cos_scat_angle->Fill(cos_scat_angle);
		h_scat_angle->Fill(TMath::ACos(cos_scat_angle)*180/TMath::Pi() );
		it++;
	};

	h_scat_angle->GetXaxis()->SetTitle("Scattering Angle [deg]");
	h_scat_angle->GetYaxis()->SetTitle("Counts per 2 deg");
	h_scat_angle->SetTitle(" ");
	h_cos_scat_angle->Draw();
	h_scat_angle->Draw();

	cout << "done." << endl;
}

void fillVectors(string V){
	ifstream in;
	in.open(V.c_str());

	string dummy;
	double dummy_x, dummy_y, dummy_z;
	int it=0;
	while(!in.eof()){
		in >> dummy >> dummy_x >> dummy_y >> dummy_z;
		volumeNames[it] = dummy;
		PMT_pos_x[it] = dummy_x;
		PMT_pos_y[it] = dummy_y;
		PMT_pos_z[it] = dummy_z;
		it++;
	}
	cout << "PMT volume list (shows the center of the corresponding scintillator crystal): " << endl;
	for (int i = 0; i < it; ++i)
	{
		cout << volumeNames[i] << "\t" << PMT_pos_x[i] << "\t" << PMT_pos_y[i] << "\t" << PMT_pos_z[i]<< endl;
	}
	in.close();
	num_PMT_Volumes=it;
}

int GetContinuousVolumeNumber(char *c){
	//cout << "---\n" << c << endl;
	int contNum=-1;
	for (int i = 0; i < num_PMT_Volumes; ++i)
	{
		if(c == volumeNames[i]) {
		//	cout << "found volume number: " << i << endl;
			contNum = i; 
		}
	}
	return contNum;
	
}