/// \file RunControl.hh
/// \brief Definition of the RunControl class

#ifndef RunControl_h
#define RunControl_h 1

#include "globals.hh"
#include <string>
#include <iostream>
#include "G4Run.hh"
#include "Analysis.hh"

/// RunControl Class
using namespace std;
class HistoManager;
class RunControl
{
	public:
		static RunControl* GetInstance(){
			if(RunControl::fInstance == NULL) {
				RunControl::fInstance = new RunControl;
			}
			return fInstance;
		}
		
		virtual ~RunControl();
		
		int NUM_Events;
		bool use_VIS;
		G4long GetSeed();
		string fileName;
		string RAWfileName;
		bool HandleParameters(int numArgs, char** argv);
		int arraySideLength_x;
		int arraySideLength_y;
		int arraySideLength_z;
		int number_of_arrays_x;
		int number_of_arrays_y;
		int number_of_arrays_z;

		G4ThreeVector startmomalpha;
		G4double theta_compt;
		G4ThreeVector source_shift;

		 G4double space; 
		 G4double spacem;
		 G4double SiPMchannel_size_xy; 
		 G4double space_x; 
		 G4double space_y; 
		 G4double PMMA_thick;
		 G4double SiPM_thick;
		 G4double ref_thick;
		 float alpha;
		 int primaryphotonscounter;
		 G4double r_source;
		 float d_ap;
		 G4double dist_PMMAap;
		 G4double dist_apsource;
		 G4double PMMA_cyl_outerRadius;

		 G4int rnocc; // required number of coincident channels

		 int SG_SideLength_x; 
		 int SG_SideLength_y; 
		 int SG_SideLength_z; 
		 double SG_ChannelSize_xy;
		 double SG_ChannelSize_z;
		 double SG_ChannelGap_xy;

		 double SG_Distance_behind_PMMA;
		 double SG_DetectorPosition_z;

		 G4double ScintYield;

		 G4double PMTs_ChannelSize_z;
		 G4double ParticleEnergy;

		 bool coverPMMASolidAngle;
		 G4double lengthOfLine;


	private:
		RunControl();
		
		static RunControl* fInstance;
		void create_RAWfileName(char* f);
		bool parametersReceived;
		HistoManager* hm;
		string info_file;
		ofstream info_out;
};
#endif
