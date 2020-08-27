#ifndef Polarized_Particle_h
#define Polarized_Particle_h 1

#include "G4ParticleDefinition.hh"
#include "G4NuclearPolarization.hh"
#include "G4PolarizationTransition.hh"

class Polarized_Particle {

public:

  Polarized_Particle(G4ParticleDefinition* def, G4int Z, G4int A, G4int J, G4double ex);
  ~Polarized_Particle();

  G4ParticleDefinition* GetDefinition() const {return part;}
  G4NuclearPolarization* GetNuclearPolarization() const {return polar;}
  G4int GetSpin() const {return spin;}
  
  std::vector< std::vector<G4complex> >& GetPolarization() const {return polar->GetPolarization();}
  
  void SetPolarization(std::vector< std::vector<G4complex> >& p) {polar->SetPolarization(p);}
  
private:

  G4ParticleDefinition* part;
  G4int spin;
  G4NuclearPolarization* polar;

};

#endif
