#ifndef Gamma_Decay_h
#define Gamma_Decay_h 1

#include "G4GeneralPhaseSpaceDecay.hh"

class Gamma_Decay : public G4GeneralPhaseSpaceDecay {

public:

  Gamma_Decay(const G4ParticleDefinition* Parent, const G4ParticleDefinition* daughter, G4double BR)
    : G4GeneralPhaseSpaceDecay(Parent->GetParticleName(),BR,2, daughter->GetParticleName(),"gamma") {

    SetParent(Parent);
    SetParentMass(Parent->GetPDGMass());
    
  }
  
  ~Gamma_Decay() {}
  
};

#endif
