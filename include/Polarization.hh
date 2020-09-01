#ifndef Polarization_h
#define Polarization_h 1

#include "G4DataInterpolation.hh"

//class Polarization_Messenger;
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

  void BuildStatisticalTensors();

  Polarization();
  ~Polarization();

private:

  void ReadTensorFile(G4String fn, std::vector<State>& states, std::vector<double>& thetas);
  void BuildProjectileTensors();
  void BuildRecoilTensors();

  std::vector<std::vector<G4DataInterpolation*>> pTensors; //Projectile statistical tensors
  std::vector<std::vector<G4DataInterpolation*>> rTensors; //Recoil statistical tensors

  G4String pFN;
  G4String rFN;

};

#endif
