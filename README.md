# Geant4_Simulation

G4 simulation of the creation and propagation of Cherenkov Light from energetic electrons in PMMA or other optically transparent materials with the aim to detect them on an array of silicon photomultipliers. 

Build simulation by executing build script:
./build_simulation.sh

---

execute simulation from inside the Build folder:
./simulation [parameters]

---

Support message printed if run execution command contains wrong parameters
Use the following parameter names and values:
parameter	type[unit]
---		---
VIS		none		visualization on
energy		int[keV]	energy of incident photon
numE		int		number of events to simulate
file		string		output file name
size_z		int[mm]	size of PMMA material
size_xy	int[mm]	size of PMMA material
ScintYield	int[#/MeV]	yield of scintillator material [number of photons per deposited energy]
---		---
Print this message any time using "./simulation -h"

---

parameter VIS allows for displaying the designed set up and zoom in and rotate. It also alows to simulate small amounts of events and to visualize the particle tracks.

