#include "Polarization.hh"
#include <fstream>

Polarization::Polarization() {

  pFN = "cd106.ten";
  rFN = "";
  
}

Polarization::~Polarization() {}

void Polarization::BuildStatisticalTensors() {

  BuildProjectileTensors();
  //BuildRecoilTensors();
  
  return;
}

void Polarization::ReadTensorFile(G4String fn, std::vector<State>& states, std::vector<double>& thetas) {

  std::ifstream file;
  file.open(fn.c_str(),std::ios::in);

  std::string line, word;
  while(std::getline(file,line)) {

    std::stringstream ls(line);
    ls >> word;
    ls >> word;
    ls >> word;

    double theta;
    std::stringstream ss(word);
    ss >> theta;

    thetas.push_back(theta);
    G4cout << "Theta: " << theta << G4endl;
    
    std::getline(file,line);
    std::getline(file,line);
    while(std::getline(file,line)) {

      if(line.empty()) {
	break;
      }
	
      std::stringstream ls1(line);

      ls1 >> word; //index
      unsigned int index;
      std::stringstream ss1(word);
      ss1 >> index;
	  
      if(states.size() < index) {
	    
	states.emplace_back();
	states.back().fStateIndex = (int)index;
	states.back().fMaxK = 0;
	   
	states.back().fTensor.emplace_back();
	states.back().fTensor.back().fK = 0;
	states.back().fTensor.back().fKappa = 0;
	   
      }

      ls1 >> word; //k
      int k;
      std::stringstream ss2(word);
      ss2 >> k;

      if(!(states.at(index-1).HasK(k))) { //new k value

	states.at(index-1).fMaxK = k;
    
	states.at(index-1).fTensor.emplace_back();
	states.at(index-1).fTensor.back().fK = k;
	states.at(index-1).fTensor.back().fKappa = 0;
	
      }

      ls1 >> word; //kappa
      int kappa;
      std::stringstream ss3(word);
      ss3 >> kappa;

      if(!(states.at(index-1).HasKKappa(k,kappa))) { //new kappa value (for current k)
    
	states.at(index-1).fTensor.emplace_back();
	states.at(index-1).fTensor.back().fK = k;
	states.at(index-1).fTensor.back().fKappa = kappa;
	
      }

      ls1 >> word; //value
      double val;
      std::stringstream ss4(word);
      ss4 >> val;

      int iindex = states.at(index-1).Index(k,kappa);
      states.at(index-1).fTensor.at(iindex).fVal.push_back(val);

      G4cout << " " << index << "\t" << k << "\t" << kappa << "\t" << val << G4endl;
	
    } //End second while(getline()) (1 ThetaCM Block)

    //std::cout << std::endl;
    
  } //End first while(getline()) (All ThetaCM Blocks)
  
  return;
}

void Polarization::BuildProjectileTensors() {

  std::vector<State> states;
  std::vector<double> thetas;
  ReadTensorFile(pFN,states,thetas);

  G4cout << "thetas.size(): " << thetas.size() << ", states.size(): " << states.size() << G4endl;
  
  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<states.size();i++) {

    pTensors.emplace_back();
    State st = states.at(i);

    for(unsigned int j=0;j<st.fTensor.size();j++) {

      pTensors.at(i).push_back(new G4DataInterpolation(&thetas[0],&(st.fTensor.at(j).fVal)[0],size,0,0));
      
    }
    
  }

  return;
  
}
