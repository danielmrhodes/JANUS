#ifndef Gamma_Decay_h
#define Gamma_Decay_h 1

//#include "G4GeneralPhaseSpaceDecay.hh"
#include "G4VDecayChannel.hh"

//class Gamma_Decay : public G4GeneralPhaseSpaceDecay {
class Gamma_Decay : public G4VDecayChannel {

public:

  Gamma_Decay(const G4ParticleDefinition* Parent, const G4ParticleDefinition* daughter, G4double BR);
  ~Gamma_Decay();

  G4DecayProducts* DecayIt(G4double);
  static G4double Pmx(G4double e, G4double p1, G4double p2);

protected:

  G4DecayProducts* TwoBodyDecayIt();
  
private:

  G4double parentmass;
  const G4double* theDaughterMasses;
  
};

#endif
