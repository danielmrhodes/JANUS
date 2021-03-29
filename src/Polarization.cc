#include "Polarization.hh"
#include "Reaction.hh"
#include <fstream>

#include "G4Clebsch.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

Polarization::Polarization() {

  offsets = {{0,0},{2,1},{4,4},{6,9}};
  
  pFN = "";
  rFN = "";

  pCalcGk = true;
  rCalcGk = true;

  G4complex val;
  val.real(1.0);
  val.imag(0.0);
  
  unpolarized.emplace_back();
  unpolarized.at(0).push_back(val);

  //deorientation effect model parameters for the projectile
  //default values
  Avji_P = 3.0;
  Gam_P = 0.02;
  Xlamb_P = 0.0345;
  TimeC_P = 3.5;
  Gfac_P = -1.0; //default is Z/A, which is assigned later
  Field_P = 6.0*std::pow(10.0,-6.0);
  Power_P = 0.6;

  //deorientation effect model parameters for the recoil
  //default values
  Avji_R = Avji_P;
  Gam_R = Gam_P;
  Xlamb_R = Xlamb_P;
  TimeC_R = TimeC_P;
  Gfac_R = -1.0; //default is Z/A, which is assigned later
  Field_R = Field_P;
  Power_R = Power_P;
  
}

Polarization::~Polarization() {}

void Polarization::BuildStatisticalTensors(G4int pZ, G4int pA, G4double pEn, G4int rZ, G4int rA,
					   std::vector<Polarized_Particle*> pLevels,
					   std::vector<Polarized_Particle*> rLevels) {
  
  G4double pM = pLevels.at(0)->GetDefinition()->GetPDGMass();
  G4double rM = rLevels.at(0)->GetDefinition()->GetPDGMass();
  
  BuildProjectileTensors(pZ,pA,pM,pEn,rM,pLevels);
  BuildRecoilTensors(pM,pEn,rZ,rA,rM,rLevels);
  
  return;
}

void Polarization::ReadTensorFile(G4String fn, std::vector<State>& states, std::vector<double>& thetas) {

  std::ifstream file;
  file.open(fn.c_str(),std::ios::in);
  if(!file.is_open()) {
    G4cout << "\nCould not open tensor file " << fn << "!!" << G4endl;
    return;
  }
  
  std::string line, word;
  while(std::getline(file,line)) {

    std::stringstream ls(line);
    ls >> word;
    ls >> word;
    ls >> word;

    G4double theta;
    std::stringstream ss(word);
    ss >> theta;

    thetas.push_back(theta);
    
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
      G4int k;
      std::stringstream ss2(word);
      ss2 >> k;

      if(!(states.at(index-1).HasK(k))) { //new k value

	if(k > states.at(index-1).fMaxK) {
	  states.at(index-1).fMaxK = k;
	}
    
	states.at(index-1).fTensor.emplace_back();
	states.at(index-1).fTensor.back().fK = k;
	states.at(index-1).fTensor.back().fKappa = 0;
	
      }

      ls1 >> word; //kappa
      G4int kappa;
      std::stringstream ss3(word);
      ss3 >> kappa;

      if(!(states.at(index-1).HasKKappa(k,kappa))) { //new kappa value (for current k)
    
	states.at(index-1).fTensor.emplace_back();
	states.at(index-1).fTensor.back().fK = k;
	states.at(index-1).fTensor.back().fKappa = kappa;
	
      }
      
      ls1 >> word; //value
      G4double val;
      std::stringstream ss4(word);
      ss4 >> val;

      int iindex = states.at(index-1).Index(k,kappa);
      states.at(index-1).fTensor.at(iindex).fVal.push_back(val);
	
    } //End second while(getline()) (1 ThetaCM Block)
    
  } //End first while(getline()) (All ThetaCM Blocks)
  
  return;
}

void Polarization::BuildProjectileTensors(G4int pZ, G4int pA, G4double pM, G4double pEn, G4double rM,
					  std::vector<Polarized_Particle*> levels) {

  if(levels.size() == 1) {
    return;
  }

  if(pFN == "") {
    G4cout << "\nNo projectile polarization\n";
    return;
  }
  
  std::vector<State> states;
  std::vector<double> thetas;
  ReadTensorFile(pFN,states,thetas);

  if(states.size()) {
    
    G4cout << "\nBuilding projectile statistical tensors from " << pFN << "\n";
    if(pCalcGk) {
      G4cout << " Gk coefficients will be applied\n";
    }
    else {
      G4cout << " No Gk coefficients\n";
    }
    
  }

  for(unsigned int i=0;i<states.size();i++) {
    pMaxK.push_back(states.at(i).fMaxK);
  }
  
  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<states.size();i++) {

    pTensors.emplace_back();
    State st = states.at(i);

    G4int spin = levels.at(i+1)->GetSpin();
    G4double time = levels.at(i+1)->GetDefinition()->GetPDGLifeTime();

    //Calculate and store Gk coefficients for one level at all ThetaCM values
    std::vector< std::array<G4double,7> > Gks;
    if(pCalcGk) {
      for(unsigned int l=0;l<size;l++) {
      
	G4double beta = Reaction::Beta_LAB(thetas.at(l),pEn,pM,rM,0.0*MeV);
	Gks.push_back(GKK(pZ,pA,beta,spin,time/ps,false));
      
      }
    }
    
    for(unsigned int j=0;j<st.fTensor.size();j++) {

      TensorComp comp = st.fTensor.at(j);
      
      if(comp.fVal.size() != size) {
	G4cout << "Projectile state " << i << " statistical tensor component " <<  comp.fK
	       << " " << comp.fKappa << " (index " << j << ") has" << comp.fVal.size()
	       << " polarization splines points, while there are " << size
	       << " thetaCM spline points! This is likely a mistake and could break thins!\n"; 
      }

      //apply Gk coefficients
      if(pCalcGk) {
	for(unsigned int l=0;l<comp.fVal.size();l++) {
	  comp.fVal.at(l) *= Gks.at(l).at(comp.fK);
	}
      }

      pTensors.at(i).push_back(new G4DataInterpolation(&thetas[0],&(comp.fVal)[0],size,0,0));
      
    }
  }

  if(pTensors.size() == levels.size() - 1) {
    G4cout << "All " << pTensors.size() << " statistical tensors built for the projectile!\n";
  }
  else {
    G4cout << pTensors.size() << " statistical tensors were built for the projectile, which has "
	   << levels.size() - 1 << " excited states! Make sure this was intentional." << G4endl;;
  }

  return;
  
}

void Polarization::BuildRecoilTensors(G4double pM, G4double pEn, G4int rZ, G4int rA, G4double rM,
				      std::vector<Polarized_Particle*> levels) {

  if(levels.size() == 1) {
    return;
  }
  
  if(rFN == "") {
    G4cout << "\nNo recoil polarization\n";
    return;
  }
  
  std::vector<State> states;
  std::vector<double> thetas;
  ReadTensorFile(rFN,states,thetas);

  if(states.size()) {
    
    G4cout << "\nBuilding recoil statistical tensors from " << rFN << "\n";
    if(rCalcGk) {
      G4cout << " Gk coefficients will be applied\n";
    }
    else {
      G4cout << " No Gk coefficients\n";
    }
    
  }
  
  for(unsigned int i=0;i<states.size();i++) {
    rMaxK.push_back(states.at(i).fMaxK);
  }

  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<states.size();i++) {

    rTensors.emplace_back();
    State st = states.at(i);

    G4int spin = levels.at(i+1)->GetSpin();
    G4double time = levels.at(i+1)->GetDefinition()->GetPDGLifeTime();

    //Calculate and store Gk coefficients for one level at all ThetaCM values
    std::vector< std::array<G4double,7> > Gks;
    if(rCalcGk) {
      for(unsigned int l=0;l<size;l++) {
      
	G4double beta = Reaction::Recoil_Beta_LAB(thetas.at(l),pEn,pM,rM,0.0*MeV);
	Gks.push_back(GKK(rZ,rA,beta,spin,time/ps,true));
	//Gks.back().at(2) *= 0.5;
	//Gks.back().at(4) *= 0.5;
	//Gks.back().at(6) *= 0.5;
	
      }
    }
    
    for(unsigned int j=0;j<st.fTensor.size();j++) {

      TensorComp comp = st.fTensor.at(j);
      
      if(comp.fVal.size() != size) {
	G4cout << "Recoil state " << i << " statistical tensor component " <<  comp.fK
	       << " " << comp.fKappa << " (index " << j << ") has" << comp.fVal.size()
	       << " polarization splines points, while there are " << size
	       << " thetaCM spline points! This is likely a mistake and could break things.\n"; 
      }

      //apply Gk coefficients
      if(rCalcGk) {
	for(unsigned int l=0;l<comp.fVal.size();l++) {
	  comp.fVal.at(l) *= Gks.at(l).at(comp.fK);
	}
      }
      
      rTensors.at(i).push_back(new G4DataInterpolation(&thetas[0],&(comp.fVal)[0],size,0,0));
      
    }
  }
  
  if(rTensors.size() == levels.size() - 1) {
    G4cout << "All " << rTensors.size() << " statistical tensors built for the recoil!\n";
  }
  else {
    G4cout << rTensors.size() << " statistical tensors were built for the recoil, which has "
	   << levels.size() - 1 << " excited states! Make sure this was intentional." << G4endl;;
  }

  return;
  
}

std::vector< std::vector<G4complex> > Polarization::GetProjectilePolarization(const G4int state,
									      const G4double th,
									      const G4double ph) {

  if((unsigned int)state >= pTensors.size()) {
    return unpolarized;
  }

  G4complex imp;
  imp.real(0.0);
  
  std::vector<G4DataInterpolation*> tensor = pTensors.at(state-1);
  G4int maxK = pMaxK.at(state-1);

  std::vector< std::vector<G4complex> > polar;
  for(G4int k=0;k<maxK+1;k++) {

    polar.emplace_back();
    if(k%2) {
      continue;
    }

    G4int off = offsets[k];
    for(G4int kp=0;kp<k+1;kp++) {

      G4complex val;
      val.real(tensor.at(off + kp)->CubicSplineInterpolation(th));
      val.imag(0);

      imp.imag(kp*ph);
      polar.at(k).push_back(std::exp(-imp)*val);
      
    }
  }

  G4complex p00 = polar[0][0];
  for(G4int k=0;k<maxK+1;k++) {

    if(k%2) {
      continue;
    }
   
    for(G4int kp=0;kp<k+1;kp++) {
      polar[k][kp] /= p00;
    }
  }
  
  return polar;
  
}

std::vector< std::vector<G4complex> > Polarization::GetRecoilPolarization(const G4int state,
									  const G4double th,
									  const G4double ph) {

  if((unsigned int)state >= rTensors.size()) {
    return unpolarized;
  }

  G4complex imp;
  imp.real(0.0);
  
  std::vector<G4DataInterpolation*> tensor = rTensors.at(state-1);
  G4int maxK = rMaxK.at(state-1);

  std::vector< std::vector<G4complex> > polar;
  for(G4int k=0;k<maxK+1;k++) {

    polar.emplace_back();
    if(k%2) {
      continue;
    }

    G4int off = offsets[k];
    for(G4int kp=0;kp<k+1;kp++) {

      G4complex val;
      val.real(tensor.at(off + kp)->CubicSplineInterpolation(th));
      val.imag(0);

      imp.imag(kp*ph);
      polar.at(k).push_back(std::exp(-imp)*val);
      
    }
  }

  G4complex p00 = polar[0][0];
  for(G4int k=0;k<maxK+1;k++) {

    if(k%2) {
      continue;
    }
   
    for(G4int kp=0;kp<k+1;kp++) {
      polar[k][kp] /= p00;
    }
  }
  
  return polar;
  
}

void Polarization::Print(const std::vector< std::vector<G4complex> > polar) const {

  G4cout << " P = [ {";
  size_t kk = polar.size();
  for(size_t k=0; k<kk; ++k) {
    if(k>0) { G4cout << "       {"; }
    size_t kpmax = (polar[k]).size();
    for(size_t kappa=0; kappa<kpmax; ++kappa) {
      if(kappa > 0) { G4cout << "}  {"; }
      G4cout << polar[k][kappa].real() << " + " 
	  << polar[k][kappa].imag() << "*i";
    }
    if(k+1 < kk) { G4cout << "}" << G4endl; }
  }
  G4cout << "} ]" << G4endl;

  return;
}

std::array<G4double,7> Polarization::GKK(const G4int iz, const G4int ia, const G4double beta,
					 const G4double spin, const G4double time, const G4bool rec) {

  //model parameters
  /*
  const G4double Avji = 3.0; // Average atomic spin
  const G4double Gam = 0.02; // FWHM of frequency distribution
  const G4double Xlamb = 0.0345; //Fluctuating state to static state transition rate
  const G4double TimeC = 3.5; //Mean time between random reorientations of fluctuating state 
  const G4double Gfac = iz/G4double(ia); //Nuclear gyromagnetic factor
  const G4double Field = 6.0*std::pow(10.0,-6.0); //Hyperfine field coefficient (600 T)
  const G4double Power = 0.6; //Hyperfine field exponent
  */

  G4double Avji;
  G4double Gam;
  G4double Xlamb;
  G4double TimeC;
  G4double Gfac = iz/G4double(ia);
  G4double Field;
  G4double  Power;

  if(rec) {
    Avji = Avji_R;
    Gam = Gam_R;
    Xlamb = Xlamb_R;
    TimeC = TimeC_R;

    if(Gfac_R > 0.0) {
      Gfac = Gfac_R;
    }
    
    Field = Field_R;
    Power = Power_R;
  }
  else {
    Avji = Avji_P;
    Gam = Gam_P;
    Xlamb = Xlamb_P;
    TimeC = TimeC_P;

    if(Gfac_P > 0.0) {
      Gfac = Gfac_P;
    }
    
    Field = Field_P;
    Power = Power_P;
  }
  
  G4int  inq, ifq;
  G4double qcen, dq, xnor;
  XSTATIC(iz,beta,inq,ifq,qcen,dq,xnor);
  
  G4double aks[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; //alpha_k
  G4double sum[3] = {0.0,0.0,0.0}; //stores sum over 6j symbols

  for(G4int j=inq;j<ifq+1;j++) {

    G4int nz = iz - j;
    G4double xji = ATS(nz);
    G4double sm = spin;
    if(spin > xji) {
      sm = xji;
    }

    G4int ncoup = G4int(2.0*sm + 0.5) + 1;
    sum[0] = 0.0;
    sum[1] = 0.0;
    sum[2] = 0.0;

    G4double valmi = spin - xji;
    if(valmi < 0.0) {
      valmi *= -1;
    }

    for(G4int mi=0;mi<ncoup;mi++) {

      G4double f = valmi + mi;
      for(G4int k=0;k<3;k++) {
	
	G4double rk = 2.0*k + 2.0;
	G4int if2 = G4int(2.0*f + 0.0001);
	G4int irk2 = G4int(2.0*rk + 0.0001);
	G4int ispin2 = G4int(2.0*spin + 0.0001);
	G4int ixji2 = G4int(2.0*xji + 0.0001);

	sum[k] += std::pow((2.0*f + 1.0)*G4Clebsch::Wigner6J(if2,if2,irk2,ispin2,ispin2,ixji2),2.0)/(2.0*xji + 1.0);
      }
      
    }

    for(G4int k=0;k<3;k++) {
      G4int k1 = 2*k;
      aks[k1] += sum[k]*std::exp(-std::pow((qcen-j)/dq,2)/2.0)/xnor;
    }
    
  } //end loop over j (charge states)

  G4double xji = Avji;
  G4double sm = spin;
  if(spin > xji) {
    sm = xji;
  }

  G4int ncoup = G4int(2.0*sm + 0.5) + 1;
  sum[0] = 0.0;
  sum[1] = 0.0;
  sum[2] = 0.0;

  G4double valmi = spin - xji;
  if(valmi < 0.0) {
    valmi *= -1;
  }

  for(G4int mi=0;mi<ncoup;mi++) {

    G4double f = valmi + mi;
    for(G4int k=0;k<3;k++) {
	
      G4double rk = 2.0*k + 2.0;
      G4int if2 = G4int(2.0*f + 0.0001);
      G4int irk2 = G4int(2.0*rk + 0.0001);
      G4int ispin2 = G4int(2.0*spin + 0.0001);
      G4int ixji2 = G4int(2.0*xji + 0.0001);

      sum[k] += std::pow((2.0*f + 1.0)*G4Clebsch::Wigner6J(if2,if2,irk2,ispin2,ispin2,ixji2),2.0)/(2.0*xji + 1.0);
    }
      
  }

  for(G4int k=0;k<3;k++) {   
    G4int k1 = 2*k + 1;
    aks[k1] += sum[k];
  }

  G4double hmean = Field*iz*std::pow(beta,Power);
  G4double wsp = 4789.0*Gfac*hmean/Avji;
  wsp *= TimeC;
  wsp *= (wsp*Avji*(Avji+1.0)/3.0);

  std::array<G4double,7> Gk = {1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  for(G4int k=0;k<3;k++) {
    
    G4int k2 = 2*k + 2;
    G4int k1 = 2*k + 1;

    G4double wrt = wsp*k2*(k2+1);
    G4double w2 = wrt;

    wrt *= (-1.0/(1.0 - aks[k2-1]));

    G4double xlam = (1.0 - aks[k2-1])*(1.0 - std::exp(wrt))/TimeC;
    G4double up = (Gam*time*aks[k1-1] + 1.0)/(time*Gam + 1.0);
    up *= (Xlamb*time);
    up += 1.0;

    G4double down = time*(xlam+Xlamb) + 1.0;
    Gk.at(k2) = up/down;
    
    G4double alp = std::sqrt(9.0*xlam*xlam + 8.0*xlam*TimeC*(w2 - xlam*xlam)) - 3.0*xlam;
    alp /= (4.0*xlam*TimeC);
    
    G4double upc = xlam*time*(down - 2.0*alp*alp*time*TimeC);
    G4double dwc = (down + alp*time)*(down + 2.0*alp*time);
    G4double ccf = 1.0 + upc/dwc;
    Gk.at(k2) *= ccf;

  }
  
  return Gk;

}

void Polarization::XSTATIC(const G4int iz, const G4double beta, G4int& id, G4int& iu, G4double& qcen,
			   G4double& dq, G4double& xnor) {

  G4double h = 1.0/(1.0 + std::pow(std::pow(iz,0.45)*0.012008/beta,5.0/3.0));
  qcen = iz*std::pow(h,0.6);
  dq = std::sqrt(qcen*(1.0-h))/2.0;

  iu = (G4int)(qcen + 3.0*dq + 0.5);
  id = (G4int)(qcen - 3.0*dq - 0.5);

  if(iu > iz) {
    iu = iz;
  }
  if(id < 1) {
    id = 1;
  }

  xnor = 0.0;
  for(G4int i=id;i<iu+1;i++) {
    xnor += std::exp(-std::pow((qcen-i)/dq,2)/2.0);
  }

  return;
  
}

G4double Polarization::ATS(G4int Ne) {
  
  if ( Ne <= 0 || Ne > 96 ) {
    return 0.0;
  }
    
  G4int mi = Ne/2 + 1;
  if (Ne%2) {
    
    if ( mi==1 || mi==2 || mi==3 || mi==6 || 
	 mi==7 || mi==10 || mi==15 || mi==16 || 
	 mi==19 || mi==24 || mi==25 || mi==28 || 
	 mi==31 || mi==35 || mi==37 || mi==40 || 
	 mi==41 || mi==44 ) {
      return 0.5;
    }
    else if ( mi==4 || mi==5 || mi==8 || mi==9 || 
	      mi==11 || mi==17 || mi==18 || mi==20 ||
	      mi==26 || mi==27 || mi==36 || mi==42 ||
	      mi==43 || mi==45 ) {
      return 1.5;
    }
    else if ( mi==12 || mi==14 || mi==21 || mi==23 ||
	      mi==32 || mi==39 ) {
      return 2.5;
    }
    else if ( mi==13 || mi==22 || mi==38 ) {
      return 4.5;
    }
    else if ( mi==29 || mi==30 || mi==48 ) {
      return 3.5;
    }
    else if ( mi==33 ) {
      return 7.5;
    }
    else if ( mi==34 ) {
      return 6.5;
    }
    else if ( mi==46 || mi==47 ) {
      return 5.5;
    }
  }
	   
  mi -= 1;
  if ( mi==4 || mi==8 || mi==17 || mi==26 || 
       mi==28 || mi==30 || mi==32 || mi==42 || 
       mi==45 || mi==48 ) {
    return 2.0;
  }
  else if ( mi==10 || mi==36 ) {
    return 1.0;
  }
  else if ( mi==12 || mi==21 || mi==37 ) {
    return 3.0;
  }
  else if ( mi==13 || mi==22 || mi==29 || mi==31 || 
	    mi==34 || mi==38 || mi==47 ) {
    return 4.0;
  }
  else if ( mi==33 ) {
    return 8.0;
  }
  else if ( mi==46 ) {
    return 6.0;
  }

  return 0.0;
   
}
