#ifndef Excitation_h
#define Excitation_h 1

#include "G4ParticleDefinition.hh"
#include "Polarized_Particle.hh"
#include "Excitation_Messenger.hh"
#include "Polarization.hh"
#include "G4DataInterpolation.hh"

class Excitation_Messenger;
class Excitation {

public:
   
  Excitation();
  ~Excitation();

  void BuildLevelSchemes(G4int projZ, G4int projA, G4int recZ, G4int recA);
  void BuildProbabilities();
  void BuildStatisticalTensors(G4int projZ, G4int projA, G4double beamEn, G4int recZ, G4int recA);

  G4int ChooseProjectileState(const G4double thetaCM);
  G4int ChooseRecoilState(const G4double thetaCM);

  G4ParticleDefinition* GetProjectileDefinition(G4int index) {return pLevels.at(index)->GetDefinition();}
  G4ParticleDefinition* GetRecoilDefinition(G4int index) {return rLevels.at(index)->GetDefinition();}

  std::vector<Polarized_Particle*> GetProjectileLevels() {return pLevels;}
  std::vector<Polarized_Particle*> GetRecoilLevels() {return rLevels;}

  G4double GetProjectileExcitation(G4int index);
  G4double GetRecoilExcitation(G4int index);

  void SetProjLSFile(G4String name) {pFN = name;}
  void SetProjPrbFile(G4String name) {pPF = name;}

  void FixProjState(G4int index) {pSelected = index;}
  void OnlyConsiderProjState(G4int index) {pConsidered = index;}
  void SetProjGSS(G4double gss) {pGSS = gss;}

  void SetRecPrbFile(G4String name) {rPF = name;}
  void SetRecLSFile(G4String name) {rFN = name;}

  void FixRecState(G4int index) {rSelected = index;}
  void OnlyConsiderRecState(G4int index) {rConsidered = index;}
  void SetRecGSS(G4double gss) {rGSS = gss;}

  void SetProjTensorFile(G4String fn) {polar->SetProjectileFile(fn);}
  void SetRecTensorFile(G4String fn) {polar->SetRecoilFile(fn);}

  void SetProjCalcGk(G4bool calc) {polar->SetProjCalcGk(calc);}
  void SetRecCalcGk(G4bool calc) {polar->SetRecCalcGk(calc);}

  void Polarize(G4int pIndex, G4int rIndex, G4double th, G4double ph);
  void Unpolarize();

  void SetProjAverageJ(G4double avj) {polar->SetProjAverageJ(avj);}
  void SetProjGamma(G4double gam) {polar->SetProjGamma(gam);}
  void SetProjLambdaStar(G4double lam) {polar->SetProjLambdaStar(lam);}
  void SetProjTauC(G4double tc) {polar->SetProjTauC(tc);}
  void SetProjGFac(G4double gf) {polar->SetProjGFac(gf);}
  void SetProjFieldCoef(G4double coef) {polar->SetProjFieldCoef(coef);}
  void SetProjFieldExp(G4double ex) {polar->SetProjFieldExp(ex);}

  void SetRecAverageJ(G4double avj) {polar->SetRecAverageJ(avj);}
  void SetRecGamma(G4double gam) {polar->SetRecGamma(gam);}
  void SetRecLambdaStar(G4double lam) {polar->SetRecLambdaStar(lam);}
  void SetRecTauC(G4double tc) {polar->SetRecTauC(tc);}
  void SetRecGFac(G4double gf) {polar->SetRecGFac(gf);}
  void SetRecFieldCoef(G4double coef) {polar->SetRecFieldCoef(coef);}
  void SetRecFieldExp(G4double ex) {polar->SetRecFieldExp(ex);}
  
private:

  inline void BuildProjectileLS(G4int Z, G4int A);
  inline void BuildRecoilLS(G4int Z, G4int A);

  inline void BuildProjectileSplines();
  inline void BuildRecoilSplines();
  
  Excitation_Messenger* messenger;
  Polarization* polar;

  G4String pFN; //Projectile LS file name
  G4String rFN; //Recoil LS file name

  G4String pPF; //Projectile probabilities file name
  G4String rPF; //Recoil probabilities file name
  
  std::vector<Polarized_Particle*> pLevels; //projectile states
  std::vector<G4DataInterpolation*> pSplines; //projectile probability splines
  
  G4int pSelected;
  G4int pConsidered;
  G4double pGSS;

  std::vector<Polarized_Particle*> rLevels; //recoil states
  std::vector<G4DataInterpolation*> rSplines; //recoil probability splines
  
  G4int rSelected;
  G4int rConsidered;
  G4double rGSS;

};

#endif
