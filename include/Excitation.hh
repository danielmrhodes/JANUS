#ifndef Excitation_h
#define Excitation_h 1

#include "G4ParticleDefinition.hh"
#include "Polarized_Particle.hh"
#include "Excitation_Messenger.hh"
#include "G4DataInterpolation.hh"

class Excitation_Messenger;
class Excitation {

public:
   
  Excitation();
  ~Excitation();

  void BuildLevelSchemes(int projZ, int projA, int recZ, int recA);
  void BuildProbabilities();

  G4int ChooseProjectileState(const G4double thetaCM);
  G4int ChooseRecoilState(const G4double thetaCM);

  G4ParticleDefinition* GetProjectileState(int index) {return pLevels.at(index)->GetDefinition();}
  G4ParticleDefinition* GetRecoilState(int index) {return rLevels.at(index)->GetDefinition();}

  std::vector<Polarized_Particle*> GetProjectileLevels() {return pLevels;}
  std::vector<Polarized_Particle*> GetRecoilLevels() {return rLevels;}

  G4double GetProjectileExcitation(int index);
  G4double GetRecoilExcitation(int index);

  void SetProjLSFile(G4String name) {pFN = name;}
  void SetProjPrbFile(G4String name) {pPF = name;}

  void SetSimpleProj() {pSimple = true;}
  void SetSimpleProjEn(G4double En) {pSimpleEn = En;}
  void SetSimpleProjLt(G4double Lt) {pSimpleLt = Lt;}
  void FixProjState(G4int index) {pSelected = index;} 

  void SetRecPrbFile(G4String name) {rPF = name;}
  void SetRecLSFile(G4String name) {rFN = name;}

  void SetSimpleRec() {rSimple = true;}
  void SetSimpleRecEn(G4double En) {rSimpleEn = En;}
  void SetSimpleRecLt(G4double Lt) {rSimpleLt = Lt;}
  void FixRecState(G4int index) {rSelected = index;} 
  
private:

  inline void BuildProjectileLS(int Z, int A);
  inline void BuildRecoilLS(int Z, int A);

  inline void BuildProjectileSplines();
  inline void BuildRecoilSplines();
  
  Excitation_Messenger* messenger;

  G4String pFN; //Projectile LS file name
  G4String rFN; //Recoil LS file name

  G4String pPF; //Projectile probabilities file name
  G4String rPF; //Recoil probabilities file name

  std::vector<Polarized_Particle*> pLevels; //projectile states
  std::vector<G4DataInterpolation*> pSplines; //projectile probability splines

  G4bool pSimple;
  G4double pSimpleEn;
  G4double pSimpleLt;
  G4int pSelected;

  std::vector<Polarized_Particle*> rLevels; //recoil states
  std::vector<G4DataInterpolation*> rSplines; //recoil probability splines
  
  G4bool rSimple;
  G4double rSimpleEn;
  G4double rSimpleLt;
  G4int rSelected;

  //G4bool popGS;

};

#endif
