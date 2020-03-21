#ifndef GammaSD_h
#define GammaSD_h 1

#include "G4VSensitiveDetector.hh"
#include "Gamma_Hit.hh"

class GammaSD : public G4VSensitiveDetector {
  
public:

  GammaSD(G4String name);
  ~GammaSD();

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*);
  void Initialize(G4HCofThisEvent*);
  void EndOfEvent(G4HCofThisEvent*);
  
private:
  
  Gamma_Hit_Collection* HC;

};

#endif
