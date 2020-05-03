#include "Gamma_Hit.hh"

G4Allocator<Gamma_Hit>* Gamma_Hit_Allocator = 0;

Gamma_Hit::Gamma_Hit() : fep(false), pfep(false) {}

Gamma_Hit::~Gamma_Hit() {}

void Gamma_Hit::SetDetSeg(G4int id) {

  seg = id%100;

  det = (id-seg)/100;
  
  return;
  
}
