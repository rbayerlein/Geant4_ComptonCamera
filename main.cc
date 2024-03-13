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
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"
#include "RunControl.hh"
#include "HistoManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

/////////////////////////////////////////
// use op novice physics list if you want to enable cherenkov light emission
// QBBC is a standard physics list, no cherenkov light emission
/////////////////////////////////////////

#include <cstddef>
#include <sys/stat.h>
#include <time.h>

//#include "QBBC.hh"
#include "OpNovicePhysicsList.hh"


//////////////////////////////////////////
// uncomment #undef G4VIS_USE ~AND~ #undef G4UI_USE if you want to simulate without graphical user interface
// the number of events to be simulated is specified further down
//////////////////////////////////////////


//#undef G4VIS_USE
//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif
//#undef G4UI_USE 
//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
//#endif

#include "Randomize.hh"

G4int NUM_Events=1;
G4bool use_VIS=false; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

/*  if(argc>1){
    string mainArgs[argc];
    for (int i = 1; i < argc; ++i)  // read in arguments into array of strings
    {
      mainArgs[i]=argv[i];
    }

    for (int i = 1; i < argc; ++i)
    {
      if(mainArgs[i]=="VIS") use_VIS=true;

      else if(mainArgs[i]=="numE") {
        if(i==argc-1) {
          cout << "WARNING - event number missing." << endl;
          return 0;
        }else if(mainArgs[i+1].find_first_not_of("0123456789 ") != std::string::npos){
          cout << "no valid number entered. default is used." << endl;
          return 0;
        }else {
          NUM_Events = std::atoi(mainArgs[i+1].c_str());
          cout<< "number of events entered:\t"<< NUM_Events << endl;
        }
      }
      else cout << mainArgs[i] << " - entry ignored." << endl;
    }

  }*/
  // Get an instance of RunControl:
  RunControl *RC = RunControl::GetInstance();

  // Pass parameters to RunControl:
  if(!RC->HandleParameters(argc, argv)) return 0; // check if correctly received. 
  NUM_Events=RC->NUM_Events;
  use_VIS=RC->use_VIS;
  cout << "Saving results in\n" << RC->fileName << endl;

  cout << "run control ready" << endl;

  // Choose the Random engine 
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(RC->GetSeed());

  // start histo manager
  // HistoManager* hm = new HistoManager(RC->fileName);

  // Construct the default run manager
  //
//#ifdef G4MULTITHREADED
//  G4MTRunManager* runManager = new G4MTRunManager;
//#else
  G4RunManager* runManager = new G4RunManager;
//#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new B1DetectorConstruction());

/////////////////////////////////
// Physics list: choose again between QBBC and OpNovice
/////////////////////////////////


 // G4VModularPhysicsList* physicsList = new QBBC;
 runManager-> SetUserInitialization(new OpNovicePhysicsList());


  //physicsList->SetVerboseLevel(0);
  //runManager->SetUserInitialization(physicsList);
    
  // User action initialization
  runManager->SetUserInitialization(new B1ActionInitialization());
  
  // Initialize visualization
  
//#ifdef G4VIS_USE
G4VisManager* visManager;

if(use_VIS){
visManager = new G4VisExecutive("all");   // alternatives: "quiet", "startup", "errors", "warnings", "confirmations", "parameters"

  	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
}
//#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

//////////////////////////////////////////////////////////////
// here you can specify a set of commands that are to be executed in each simulation run
// set verbosity, initialize kernel, set parameters for histogram output
//////////////////////////////////////////////////////////////





  UImanager->ApplyCommand("/control/verbose 2");  // 0=silent; 1=valid commands are shown; 2= comment lines are shown additionally
  UImanager->ApplyCommand("/run/verbose 2");      // 0=silent; 1=main topics; 2=main topics and run summary
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/run/initialize");     // initialize G4 kernel - this is crucial for the simulation to run!


///////////////////////////////////////////////
// set number of bins, xmin, xmax, unit for histograms
// not all of these histograms are currently filled, check HistoManager for details
///////////////////////////////////////////////


////////
// 1D histograms:
////////
//  UImanager->ApplyCommand("/analysis/h1/set 0  10 0 10 none");
  UImanager->ApplyCommand("/analysis/h1/set 1  10 0 10 none");
  UImanager->ApplyCommand("/analysis/h1/set 2  200 -50  50 mm");
  // UImanager->ApplyCommand("/analysis/h1/set 3  1000 0.  1600. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 4  1000 0.  1600. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 5  1000 0.  1600. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 6  1000 0.  1200. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 7  1000 0.  1200. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 8  1000 0.  1200. keV");
  // UImanager->ApplyCommand("/analysis/h1/set 9  1000 0. 25. mm"); 
  // UImanager->ApplyCommand("/analysis/h1/set 10  1000 0 700");
  // UImanager->ApplyCommand("/analysis/h1/set 11  1000 200.  900. nm");
  // UImanager->ApplyCommand("/analysis/h1/set 12  1000 0.  2500 keV");
  // UImanager->ApplyCommand("/analysis/h1/set 13  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 14  300 0. 300 none");
  // UImanager->ApplyCommand("/analysis/h1/set 15  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 16  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 17  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 18  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 19  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 20  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 21  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 22  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 23  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 24  300 0. 300. none");
  // UImanager->ApplyCommand("/analysis/h1/set 25  300 0. 300. none");
  UImanager->ApplyCommand("/analysis/h1/set 26  300 0. 300 none");
  UImanager->ApplyCommand("/analysis/h1/set 27  300 0. 300 none"); //testhistogramm
  // UImanager->ApplyCommand("/analysis/h1/set 28  300 0. 300 none"); //testhistogramm
  // UImanager->ApplyCommand("/analysis/h1/set 29  300 0. 300 none"); //testhistogramm
  // UImanager->ApplyCommand("/analysis/h1/set 30  5 -2 2 none"); //primaryphotonscounter
  UImanager->ApplyCommand("/analysis/h1/set 31  1 1 2 none"); //gammas
  UImanager->ApplyCommand("/analysis/h1/set 32  300 0. 300 none"); //darkcount
  UImanager->ApplyCommand("/analysis/h1/set 33  2000 0. 2000 keV"); 
  // UImanager->ApplyCommand("/analysis/h1/set 34  2000 0. 2000 keV"); 
  UImanager->ApplyCommand("/analysis/h1/set 35  180 0. 180 none"); 
  UImanager->ApplyCommand("/analysis/h1/set 36  2000 0. 2000 keV"); 
  UImanager->ApplyCommand("/analysis/h1/set 37 450 100 1000 none");
  UImanager->ApplyCommand("/analysis/h1/set 38 400 0 20000 none");

////////
// 2D histograms
////////

 //  UImanager->ApplyCommand("/analysis/h2/set 0  8 -13.  13. mm none linear 8 -13. 13 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 1  1000 0. 1000. nm none linear 1000 0. 25. mm none linear");
 // UImanager->ApplyCommand("/analysis/h2/set 2  8 -4. 4. none none linear 8 -4. 4. none none linear");
 // UImanager->ApplyCommand("/analysis/h2/set 3  8 -4. 4. none none linear 8 -4. 4. none none linear");
 // UImanager->ApplyCommand("/analysis/h2/set 4  8 -4. 4. none none linear 8 -4. 4. none none linear");
 // UImanager->ApplyCommand("/analysis/h2/set 5  8 -4. 4. none none linear 8 -4. 4. none none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 6  8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 7  8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 8  8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 9  8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 10 8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 11 8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 12 8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 13 8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");
 //  UImanager->ApplyCommand("/analysis/h2/set 14 8 -12.5 12.5 mm none linear 8 -12.5 12.5 mm none linear");

/////////
// set x axis title [title must be one word]
/////////
  UImanager->ApplyCommand("/analysis/h1/setXaxis 11 wavelength[nm]");
//  UImanager->ApplyCommand("/analysis/h1/setXaxis 1 electron_angle_[degrees]");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 2 production_depth_[mm]");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 10 number_of_produced_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 12 initial_Compton_electron_energy_[keV]");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 13 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 14 number_of_active_channels");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 15 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 16 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 17 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 18 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 19 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 26 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 27 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 32 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 35 photon_scattering_angle");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 37 Wavelength_[nm]");
  UImanager->ApplyCommand("/analysis/h1/setXaxis 38 Number_of_emitted_scintillation_photons");

  UImanager->ApplyCommand("/analysis/h2/setXaxis 4 Compton_scattering_angle");
  UImanager->ApplyCommand("/analysis/h2/setYaxis 4 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h2/setZaxis 4 counts");
  UImanager->ApplyCommand("/analysis/h2/setXaxis 5 electron_energy");
  UImanager->ApplyCommand("/analysis/h2/setYaxis 5 number_of_detected_photons");
  UImanager->ApplyCommand("/analysis/h2/setZaxis 5 counts");

//////////////////////////////////////////////////////////////
// here, a macro is chosen depending on wether or not visualisation is enabled.
// vis.mac : visualisation enabled
// nogui.mac : visualisation disabled
// IMPORTANT: any changes to the particle source have to be made in each of these macros individually
//////////////////////////////////////////////////////////////



if (use_VIS)
{

  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
                                              // GUI AND VIS
  UImanager->ApplyCommand("/control/execute vis.mac"); 

/*  if (ui->IsGUI()) {
  UImanager->ApplyCommand("/control/execute gui.mac");
  }*/

/*#else                                                         // GUI NO VIS
  UImanager->ApplyCommand("/control/execute init.mac");*/
  cout << "MAIN:\tUI session start." << endl;
  ui->SessionStart();

  cout << "MAIN:\tUI session end." << endl;
  
  delete ui; 
}else{

            
                                            //NO GUI AND NO VIS
  G4cout << "No G4UI_USE" << G4endl;
 // UImanager->ApplyCommand("/control/execute nogui.mac");

////////////////////////////////////////////////////////////////
// THIS IS WHERE YOU DEFINE THE NUMBER OF EVENTS IN nogui.mac-MODE!!!
////////////////////////////////////////////////////////////////


  char command[128];
  sprintf(command, "%s%d", "/run/beamOn ", NUM_Events );
  UImanager->ApplyCommand(command);
}
/*
// alternative: Define a macro that contains run and beamOn commands as follows:
  UImanager->ApplyCommand("/control/execute MAKRO.mac"); 
// one needs to add the macro's name to CMakeLists under set(SIMULATION_SCRIPTS ... )  
*/

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  if(use_VIS){
    cout << "deleting VIS Manager" << endl;
    delete visManager;
    cout << "done" << endl;
  }
//  delete hm;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
