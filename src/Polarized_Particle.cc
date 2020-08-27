#include "Polarized_Particle.hh"

Polarized_Particle::Polarized_Particle(G4ParticleDefinition* def, G4int Z, G4int A, G4int J,
				       G4double ex) : part(def), spin(J) {

  polar = new G4NuclearPolarization(Z,A,ex);
  
}

Polarized_Particle::~Polarized_Particle() {}

