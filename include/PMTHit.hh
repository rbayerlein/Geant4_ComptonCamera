#ifndef PMTHit_h
#define PMTHit_h 1
 
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
 
class PMTHit : public G4VHit{
 
   public:
      /* the hit must be given an ID */
      PMTHit(const G4int kHitIDnumber);
      ~PMTHit();
      void Print();
      void Draw();
 
   public:
      /* SET methods */
      void setName       (const G4String STR) {mName = STR;};
      void setEdep       (const G4double E)   {mEdep = E;};
      void setGlobalTime (const G4double T)   {mGlobalTime = T;};
      void setVolume     (const G4String V)   {mVolume = V;};
      void setEnergy     (const G4double En)  {mEnergy = En;};
      void setTrackID    (const G4int ID)     {mTrackID = ID;};
 
      /* GET methods */
      G4String getName()        const {return mName;};
      G4double getEdep()        const {return mEdep;};
      G4double getGlobalTime()  const {return mGlobalTime;};
      G4String getVolume()      const {return mVolume;};
      G4double getEnergy()      const {return mEnergy;};
      G4int    getTrackID()     const {return mTrackID;};
 
   private:
      G4int     gHitID;
      G4String  mName = "unknown";
      G4double  mEdep = -1.0;
      G4double  mGlobalTime = -1.0;
      G4String  mVolume = "unknownVol";
      G4double  mEnergy = -1.0;
      G4int     mTrackID = -1;
};
 
/* define the hit collection type we used in the SD header */
typedef G4THitsCollection<PMTHit> hitCollection;
 
#endif