#ifndef Gamma_Hit_h
#define Gamma_Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class Gamma_Hit : public G4VHit {
  
public:

  Gamma_Hit();
  ~Gamma_Hit();

  inline void* operator new(size_t);
  inline void operator delete(void *aHit);

  inline void SetEdep(G4double de) {edep = de;}
  inline void SetPos(G4ThreeVector vec) {pos = vec;}
  
  inline void SetDet(G4int d) {det = d;}
  inline void SetSeg(G4int sg) {seg = sg;}
  
  void SetDetSeg(G4int id);
  
  G4double GetEdep() const {return edep;}
  G4ThreeVector GetPos() const {return pos;}
 
  G4int GetSegment() const {return seg;}
  G4int GetDetector() const {return det;}
  
private:
  
  G4double edep;
  G4ThreeVector pos;

  G4int det, seg;
  
};

typedef G4THitsCollection<Gamma_Hit> Gamma_Hit_Collection;

extern G4Allocator<Gamma_Hit>* Gamma_Hit_Allocator;

inline void* Gamma_Hit::operator new(size_t) {
  
  if(!Gamma_Hit_Allocator) {
    Gamma_Hit_Allocator = new G4Allocator<Gamma_Hit>;
  }
  
  return (void*)Gamma_Hit_Allocator->MallocSingle();
}

inline void Gamma_Hit::operator delete(void* aHit) {
  Gamma_Hit_Allocator->FreeSingle((Gamma_Hit*)aHit);
}

#endif
