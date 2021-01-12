#ifndef Polarization_h
#define Polarization_h 1

#include "Polarized_Particle.hh"
#include "G4DataInterpolation.hh"
#include <map>

class Polarization {

public:
  
  struct TensorComp {

    int fK;
    int fKappa;  
    std::vector<double> fVal;
    
  };
  
  struct State {

    int fStateIndex;
    int fMaxK;
  
    std::vector<TensorComp> fTensor;

    bool HasK(int k) {

      for(unsigned int i=0;i<fTensor.size();i++) {
	if(fTensor.at(i).fK == k) {
	  return true;
	}
      }

      return false;
    }

    bool HasKKappa(int k, int kappa) {

      for(unsigned int i=0;i<fTensor.size();i++) {
	if(fTensor.at(i).fK == k && fTensor.at(i).fKappa == kappa) {
	  return true;
	}
      }

      return false;
    }

    int Index(int k, int kappa) {

      for(unsigned int i=0;i<fTensor.size();i++) {
	if(fTensor.at(i).fK == k && fTensor.at(i).fKappa == kappa) {
	  return i;
	}
      }

      return -1;
    
    }
  
  };

  Polarization();
  ~Polarization();

  void BuildStatisticalTensors(G4int projZ, G4int projA, G4double projEn, G4int recZ, G4int recA,
			       std::vector<Polarized_Particle*> projLevels,
			       std::vector<Polarized_Particle*> recLevels);

  std::vector< std::vector<G4complex> > GetProjectilePolarization(const G4int state, const G4double th,
								  const G4double ph);
  std::vector< std::vector<G4complex> > GetRecoilPolarization(const G4int state, const G4double th,
							      const G4double ph);

  void SetProjectileFile(G4String fn) {pFN = fn;}
  void SetRecoilFile(G4String fn) {rFN = fn;}

  void SetProjCalcGk(G4bool calc) {pCalcGk = calc;}
  void SetRecCalcGk(G4bool calc) {rCalcGk = calc;}

  void SetAverageJ(G4double avj) {Avji = avj;}
  void SetGamma(G4double gam) {Gam = gam;}
  void SetLambdaStar(G4double lam) {Xlamb = lam;}
  void SetTauC(G4double tc) {TimeC = tc;}
  void SetGFacMult(G4double mult) {GfacMult = mult;}
  void SetFieldCoef(G4double coef) {Field = coef;}
  void SetFieldExp(G4double ex) {Power = ex;}

private:

  void ReadTensorFile(G4String fn, std::vector<State>& states, std::vector<double>& thetas);
  void BuildProjectileTensors(G4int projZ, G4int projA, G4double projM, G4double projEn, G4double recM,
			      std::vector<Polarized_Particle*> projLevels);
  
  void BuildRecoilTensors(G4double projM, G4double projEn, G4int recZ, G4int recA, G4double recM,
			  std::vector<Polarized_Particle*> recLevels);

  //calculate Gk coefficients
  std::array<G4double,7> GKK(const G4int iz, const G4int ia, const G4double beta, const G4double spin,
			     const G4double time);
  G4double ATS(const G4int Ne);
  void XSTATIC(const G4int iz, const G4double beta, G4int& id, G4int& iu, G4double& qcen, G4double& dq,
	       G4double& xnor);

  inline void Print(const std::vector< std::vector<G4complex> > polar) const;

  std::vector< std::vector<G4complex> > unpolarized; //Unpolarized statistical tensor

  std::map<G4int,G4int> offsets;
  std::vector< std::vector<G4DataInterpolation*> > pTensors; //Projectile statistical tensors
  std::vector< std::vector<G4DataInterpolation*> > rTensors; //Recoil statistical tensors

  std::vector<G4int> pMaxK; //Max k for projectile states
  std::vector<G4int> rMaxK; //Max k for recoil states

  G4String pFN; //Projectile tensor file name
  G4String rFN; //Recoil tensor file name

  G4bool pCalcGk; //Flag to calculate depolarization coefficients for projectile
  G4bool rCalcGk; //Flag to calculate depolarization coefficients for recoil

  //model parameters
  G4double Avji; // Average atomic spin
  G4double Gam; // FWHM of frequency distribution
  G4double Xlamb; //Fluctuating state to static state transition rate
  G4double TimeC; //Mean time between random reorientations of fluctuating state 
  G4double GfacMult; //Scaling factor for nuclear gyromagnetic ratio (g = Z/A)
  G4double Field; //Hyperfine field coefficient
  G4double Power; //Hyperfine field exponent

};

#endif
