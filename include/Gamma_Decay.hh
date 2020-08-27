#ifndef Gamma_Decay_h
#define Gamma_Decay_h 1

#include "G4VDecayChannel.hh"
#include "Polarized_Particle.hh"
#include "G4PolarizationTransition.hh"
#include "G4NuclearPolarization.hh"

class Gamma_Decay : public G4VDecayChannel {

public:

  Gamma_Decay(Polarized_Particle* Parent, Polarized_Particle* daughter, G4double BR);
  ~Gamma_Decay();

  G4DecayProducts* DecayIt(G4double);
  static G4double Pmx(G4double e, G4double p1, G4double p2);
  
private:

  Gamma_Decay(G4ParticleDefinition* Parent, G4ParticleDefinition* daughter, G4double BR);
  
  G4double parentmass;
  const G4double* theDaughterMasses;

  G4PolarizationTransition* trans;

  Polarized_Particle* pParent;
  Polarized_Particle* pDaughter;
  
};

#endif
