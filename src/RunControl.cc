/// \file RunControl.cc
/// \brief Implementation of the RunControl class

#include <sys/stat.h>   // stat
#include <stdbool.h>    // bool type

#include "RunControl.hh"
#include "HistoManager.hh"
#include "G4SystemOfUnits.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunControl * RunControl::fInstance =0;

bool file_exists(char *f);
void create_RAWfileName(char *f);
void PrintHelp();

RunControl::RunControl()
:hm(0)
{
	parametersReceived = false;
	use_VIS=false;
	NUM_Events=1;
	fileName = "/home/rbayerlein/Documents/Analyses/FromSiegen/Scattered_Gamma_Detector/TEST.root";
	RAWfileName = "/home/rbayerlein/Documents/Analyses/FromSiegen/Scattered_Gamma_Detector/TEST";
//	hm = new HistoManager(fileName);
	arraySideLength_x = 4;			//remember to change the size of the histogamms in the HistoManager
	arraySideLength_y = 4;
	arraySideLength_z = 4;
	number_of_arrays_x = 2;
	number_of_arrays_y = 2;
	number_of_arrays_z = 2;			// set 0 for simulating without sidearrays

	startmomalpha = G4ThreeVector(1,1,1);
	theta_compt = -50.;
	source_shift = G4ThreeVector(0*mm, 0*mm, 0);		// shift in x-direction perpendicular to beam line
	
	space = 0.01*cm; 			//space between SiPM-channels 
	spacem = 0.05*cm;

	SiPMchannel_size_xy = 0.3*cm; 		//xy-size SiPM-channels (0.3cm)
	space_x = 0.*cm;			//overlapp PMMA above SiPM x
	space_y = 0.*cm;			//overlapp PMMA above SiPM y

	PMMA_thick = 24.7*mm;			//old 8mm, 8*8*8 cube 24.7mm
	SiPM_thick = 0.5*cm;
	ref_thick = 0.1*mm;      
	
	alpha =0;							//angle of incidence in deg
	d_ap = 6.25;  				            //diameter of the aperture in mm, 10 in real, 0 to simulate without momentum distribution
	r_source = 0;						//radius of the extended source, 0 to simulate with point source, 1 to have an half sphere, rember to chance the stat position in the PGA

	dist_PMMAap = 7.5+PMMA_thick/(2*mm);	//distance between PMMA and apreture in mm (middle of the PMMA). In set-up: 36 + 8mm/2. 
										//Removed this term (19/10/2020): -(PMMA_thick-8)
	dist_apsource = 7.5;			        //distance between appreture and the source (middle of the source). In set-up: 25
	primaryphotonscounter = 10000;
	PMMA_cyl_outerRadius = 1.5/2.*cm;

	rnocc = 5;

// parameters for SG Detector

	SG_SideLength_x = 4;
	SG_SideLength_y = 4;
	SG_SideLength_z = 4;

	SG_ChannelSize_xy = 15*mm;
	SG_ChannelSize_z = 15*mm;

	SG_ChannelGap_xy = 2.0*mm;

	SG_Distance_behind_PMMA=0.5*SG_SideLength_x*(SG_ChannelGap_xy+SG_ChannelSize_xy);
	SG_DetectorPosition_z = SG_Distance_behind_PMMA;

	ScintYield = 6300; 	// given as number of photons per MeV; for BGO: 10000, for LaBr3: 63000

	PMTs_ChannelSize_z = 0.05*cm;

	ParticleEnergy=1500*keV; // in keV

	coverPMMASolidAngle = true;		// determines if the incident particles are sent towards the PMMA surface to cover the full solid 
											// angle of its front face. Otherwise, a point source or extended source with apterture is used
	lengthOfLine = 0*mm;					// in case you want to have an extended source distribution along a line in x-direction, use this. 
											// only works, if coverPMMASolidAngle is set to TRUE

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	
RunControl::~RunControl(){
	delete hm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4long RunControl::GetSeed(){
	time_t t = time(0);   // get time now
	struct tm * now = localtime( & t );
	G4long seed = ((now->tm_hour)*(now->tm_min)*(now->tm_sec))*now->tm_mday;
	G4cout << "RC:\tSeed for Random Machine: " << seed << G4endl;
	return seed;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool RunControl::HandleParameters(int numArgs, char** argv){
	if(parametersReceived) {
		cout << "RC:\tCannot overwrite parameters" << endl;
		return true;
	}

	cout << "RC:\tNumber of Arguments passed: " << numArgs << endl;
	if(numArgs>1){
		string mainArgs[numArgs];
		for (int i = 1; i < numArgs; ++i)  // read in arguments into array of strings
		{
		  mainArgs[i]=argv[i];
		}

		for (int i = 1; i < numArgs; ++i)
		{
		  if(mainArgs[i]=="VIS") use_VIS=true;
		  else if(mainArgs[i] == "-h") {
		  	PrintHelp();
		  	return false;
		  }
		  else if(mainArgs[i]=="numE") {
		    if(i==numArgs-1) {
		      cout << "WARNING - event number missing." << endl;
		      return false;
		    }else if(mainArgs[i+1].find_first_not_of("0123456789 ") != std::string::npos){
		      cout << "no valid number entered. default is being used: NUM_Events = "<< NUM_Events << endl;
		      return false;
		    }else {
		      NUM_Events = std::atoi(mainArgs[i+1].c_str());
		      cout<< "number of events entered:\t"<< NUM_Events << endl;
		    }
		  }
		  else if(mainArgs[i] == "file"){
		  	if(i==numArgs-1){
		  		cout << "File name missing. " << endl;
		  		return false;
		  	}else{
		  		// here check if file exists and if it does, add something to it in the end. 
		  		char tempF[2048];
		  		sprintf(tempF, mainArgs[i+1].c_str());
		  		cout << "tempF " << tempF << endl;
		  		create_RAWfileName(tempF);

		  		char tempF2[2048];
		  		sprintf(tempF2, "%s%s", RAWfileName.c_str(), ".root");
		  		cout << "tempF2 " << tempF2 << endl;

		  		fileName = tempF2;

		  		int increment = 2;
		  		while(file_exists(tempF2)){
		  			cout << "File exists. Will add _" << increment << " to the name." << endl;
		  			stringstream ss_fileName, ss_RAWfileName;
		  			if(increment == 2){
			  			ss_fileName << RAWfileName << "_" << increment << ".root" ;
			  			ss_RAWfileName << RAWfileName  << "_" << increment;
			  		}else{
			  			ss_fileName << RAWfileName.substr(0,RAWfileName.length()-2) << "_" << increment << ".root" ;
			  			ss_RAWfileName << RAWfileName.substr(0,RAWfileName.length()-2)   << "_" << increment;
					}
		  			fileName = ss_fileName.str();
		  			RAWfileName = ss_RAWfileName.str();
		  			increment++;
		  			sprintf(tempF2, "%s%s", RAWfileName.c_str(), ".root");

		  		}


		  	}
		  }
		  else if(mainArgs[i]=="size_z") {
		    if(i==numArgs-1) {
		      cout << "WARNING - size_z missing." << endl;
		      return false;
		      PrintHelp();
		    }else if(mainArgs[i+1].find_first_not_of("0123456789. ") != std::string::npos){
		      cout << "no valid number entered for size_z: "<< mainArgs[i+1] << endl;
		      return false;
		      PrintHelp();
		    }else {
		      SG_ChannelSize_z = std::atoi(mainArgs[i+1].c_str())*mm;
		      cout<< "size_z entered:\t"<< SG_ChannelSize_z << endl;
		  	}
		  }
/*		  else if(mainArgs[i]=="quartz_z") {
		    if(i==numArgs-1) {
		      cout << "WARNING - quartz thickness missing." << endl;
		      return false;
		      PrintHelp();
		    }else if(mainArgs[i+1].find_first_not_of("0123456789. ") != std::string::npos){
		      cout << "no valid number entered for Quartz_Thickness: "<< mainArgs[i+1] << endl;
		      return false;
		      PrintHelp();
		    }else {
		      Quartz_Thickness = std::atoi(mainArgs[i+1].c_str())*mm;
		      cout<< "Quartz_Thickness entered:\t"<< Quartz_Thickness << endl;
		  	}
		  }*/
		  else if(mainArgs[i]=="ScintYield") {
		    if(i==numArgs-1) {
		      cout << "WARNING - ScintYield missing." << endl;
		      return false;
		      PrintHelp();
		    }else if(mainArgs[i+1].find_first_not_of("0123456789. ") != std::string::npos){
		      cout << "no valid number entered for ScintYield: "<< mainArgs[i+1] << endl;
		      return false;
		      PrintHelp();
		    }else {
		      ScintYield = std::atoi(mainArgs[i+1].c_str())*mm;
		      cout<< "ScintYield entered:\t"<< ScintYield << endl;
		  	}
		  }
		  else if(mainArgs[i]=="size_xy") {
		    if(i==numArgs-1) {
		      cout << "WARNING - size_xy missing." << endl;
		      return false;
		      PrintHelp();
		    }else if(mainArgs[i+1].find_first_not_of("0123456789. ") != std::string::npos){
		      cout << "no valid number entered for size_xy: "<< mainArgs[i+1] << endl;
		      return false;
		      PrintHelp();
		    }else {
		      SG_ChannelSize_xy = std::atoi(mainArgs[i+1].c_str())*mm;
		      cout<< "size_xy entered:\t"<< SG_ChannelSize_xy << endl;
		  	}
		  }
		  else if(mainArgs[i]=="energy") {
		    if(i==numArgs-1) {
		      cout << "WARNING - energy missing." << endl;
		      return false;
		      PrintHelp();
		    }else if(mainArgs[i+1].find_first_not_of("0123456789. ") != std::string::npos){
		      cout << "no valid number entered for energy: "<< mainArgs[i+1] << endl;
		      return false;
		      PrintHelp();
		    }else {
		      ParticleEnergy = std::atoi(mainArgs[i+1].c_str())*keV;
		      cout<< "energy entered:\t"<< ParticleEnergy << endl;
		  	}
		  }
		  
		  else cout << "...taking next entry..." << endl;
		}// end of for loop
	} // end of "if(numArgs > 1)"

	parametersReceived=true;

	stringstream ss_info_out;
	ss_info_out << RAWfileName << "_info.txt";
	info_file = ss_info_out.str();
	info_out.open(info_file.c_str());
	info_out << "Output data:\n" << RAWfileName << "_*.*" << endl;
	info_out << "\n####################" << endl;
	info_out << "Geometry Information SiPM arrays:" << endl;
	info_out << "arraySideLength_x\t" << arraySideLength_x << endl;
	info_out << "arraySideLength_y\t" << arraySideLength_y << endl;
	info_out << "arraySideLength_z\t" << arraySideLength_z << endl;
	info_out << "number_of_arrays_x\t" << number_of_arrays_x << endl;
	info_out << "number_of_arrays_y\t" << number_of_arrays_y << endl;
	info_out << "number_of_arrays_z\t" << number_of_arrays_z << endl;
	info_out << "SiPMchannel_size_xy\t" << SiPMchannel_size_xy << endl;
	info_out << "Space between channels:\t" << space << endl;

	info_out << "\n####################" << endl;
	info_out << "Geometry Information PMT arrays:" << endl;
	info_out << "SG_Distance_behind_PMMA\t" << SG_Distance_behind_PMMA << endl;
	info_out << "SG_SideLength_x\t" << SG_SideLength_x << endl;
	info_out << "SG_SideLength_y\t" << SG_SideLength_y << endl;
	info_out << "SG_SideLength_z\t" << SG_SideLength_z << endl;
	info_out << "SG_ChannelSize_xy\t" << SG_ChannelSize_xy << endl;
	info_out << "SG_ChannelSize_z\t" << SG_ChannelSize_z << endl;
	info_out << "SG_ChannelGap_xy\t" << SG_ChannelGap_xy << endl;
	
	info_out << "\n####################" << endl;
	info_out << "Pre-Analysis info:" << endl;
	if(coverPMMASolidAngle) {
		info_out << "Solid angle between source and PMMA covered:\t" << coverPMMASolidAngle << endl;
		info_out << "Extended line - length\t" << lengthOfLine << endl;
	}
	info_out << "NUM_Events:\t" << NUM_Events << endl;
	info_out << "Energy:\t" << ParticleEnergy << endl;
	info_out << "rnocc\t " << rnocc << endl;
	info_out << "ScintYield\t" << ScintYield << endl;
	info_out << "Angle of incidence\t" << alpha << endl;
	info_out << "Distance from source to PMMA\t" << dist_PMMAap +dist_apsource - PMMA_thick/2 << endl;
	info_out << "Aperture diameter\t" << d_ap << endl;
	info_out << "Distance from source to aperture\t" << dist_apsource << endl;
	info_out.close();


	cout << "RC:\tParameters successfully received and handled." << endl;

	hm = new HistoManager(fileName);
	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


bool file_exists(char *f){
	struct stat buffer;
	return (stat (f,&buffer) == 0);
}

void RunControl::create_RAWfileName(char *f){
	string str=f;
	cout << "str " << str << endl;
	size_t found = str.find_last_of(".");
	if(found != string::npos){
		if(found+1==str.length()){
			RAWfileName = str.substr(0,found-4);
		}
	}else{
		RAWfileName = str;
	}
	cout << "RAWfileName: " << RAWfileName << endl;
}

void PrintHelp(){
	cout << "###\nSupport message.\nUse the following parameter names and values:" << endl;
	cout 	<< "parameter\ttype[unit]\n"
			<< "---\t\t---" << "\n"
			<< "VIS\t\tnone" << "\n"
			<< "energy\t\tint[keV]" << "\n"
			<< "numE\t\tint" << "\n"
			<< "file\t\tstring" << "\n"
			<< "size_z\t\tint[mm]" << "\n"
			<< "size_xy\t\tint[mm]" << "\n"
			<< "ScintYield\tint[#/MeV]" << "\n";

	cout << "Print this message any time using -h as argument." << endl;
}
