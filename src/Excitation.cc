#include "Excitation.hh"
#include "Gamma_Decay.hh"

#include "G4ProcessManager.hh"
#include "G4IonTable.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

Excitation::Excitation() {

  messenger = new Excitation_Messenger(this);
  polar = new Polarization();
  
  pFN = "";
  pPF = "";

  pSimple = false;
  pSimpleEn = 0.0*MeV;;
  pSimpleLt = 0.0*ps;
  pSelected = -1;
  pConsidered = 0;
  pGSS = 0.0;
    
  rFN = "";
  rPF = "";

  rSimple = false;
  rSimpleEn = 0.0*MeV;;
  rSimpleLt = 0.0*ps;
  rSelected = -1;
  rConsidered = 0;
  rGSS = 0.0;
  
}

Excitation::~Excitation() {

  delete messenger;
  delete polar;
}

void Excitation::BuildLevelSchemes(G4int pZ, G4int pA, G4int rZ, G4int rA) {

  BuildProjectileLS(pZ,pA);
  BuildRecoilLS(rZ,rA);
  
  return;
}

void Excitation::BuildProbabilities() {

  BuildProjectileSplines();
  BuildRecoilSplines();
  
  return;
}

void Excitation::BuildStatisticalTensors(G4int pZ, G4int pA, G4double pEn, G4int rZ, G4int rA) {
  
  polar->BuildStatisticalTensors(pZ,pA,pEn,rZ,rA,pLevels,rLevels);
  
  return;
}

void Excitation::BuildProjectileLS(G4int Z, G4int A) {
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition* projGS = table->GetIon(Z,A,0.0*MeV);
  projGS->SetPDGLifeTime(-1.0);
  
  Polarized_Particle* polGS = new Polarized_Particle(projGS,Z,A,pGSS,0.0*MeV);
  pLevels.push_back(polGS);

  if(pSimple) {

    G4cout << "\nBuilding simple projectile level scheme!\n";
    
    G4ParticleDefinition* part = table->GetIon(Z,A,pSimpleEn);
    part->SetPDGLifeTime(pSimpleLt);
    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,2.0,pSimpleEn);
    part->GetDecayTable()->Insert(new Gamma_Decay(ppart,pLevels.at(0),1,2,0,0));

    pLevels.push_back(ppart);

    G4cout << " 1 " << pSimpleEn/keV << " 2 " << pSimpleLt/ps << " 1\n  0 1 2 0 0\nSuccess!" << G4endl;  

    return;
    
  }

  if(pFN == "") {
    G4cout << "\nNo projectile excitations." << G4endl;
    return;
  }
  
  std::ifstream file;
  file.open(pFN.c_str(),std::ios::in);
  
  if(!file.is_open()) {
    G4cout << "\nCould not open projectile level scheme file " << pFN << "! No levels will be built!"
	   << G4endl;
    return;
  }

  G4cout << "\nBuilding projectile level scheme from " << pFN << G4endl;
  
  std::string line, word;
  while(std::getline(file,line)) {
      
    std::stringstream linestream1(line);
    linestream1 >> word;

    G4int state_index;
    std::stringstream ss(word);
    ss >> state_index;

    linestream1 >> word;

    G4double energy;
    std::stringstream ss1(word);
    ss1 >> energy;
    energy *= keV;
    
    linestream1 >> word;
    
    G4double spin;
    std::stringstream ss2(word);
    ss2 >> spin;

    linestream1 >> word;

    G4double lifetime;
    std::stringstream ss3(word);
    ss3 >> lifetime;
    lifetime *= ps;

    linestream1 >> word;

    G4int nbr;
    std::stringstream ss4(word);
    ss4 >> nbr;
    
    if(nbr == 0) {
      G4cout << "Problem reading projectile level scheme file " << pFN
	     << "! No decay braches declared for state " << state_index << G4endl;
    }
	
    G4ParticleDefinition* part = table->GetIon(Z,A,energy);
    if(!pConsidered || state_index == pConsidered) {
      
      part->SetPDGLifeTime(lifetime);
      part->SetDecayTable(new G4DecayTable());
      part->GetProcessManager()->SetParticleType(part);
      part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);
      
    }
    else {
      part->SetPDGLifeTime(-1.0);
    }

    G4cout << " " << state_index << " " << energy/keV << " " << spin << " " << lifetime/ps << " " << nbr
	   << G4endl;
    
    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,spin,energy);
    for(int i=0;i<nbr;i++) {

      std::getline(file,line);
      std::stringstream linestream2(line);
      
      linestream2 >> word;

      G4int index;
      std::stringstream ss5(word);
      ss5 >> index;

      linestream2 >> word;
      
      G4double BR;
      std::stringstream ss6(word);
      ss6 >> BR;

      linestream2 >> word;

      G4int L0;
      std::stringstream ss7(word);
      ss7 >> L0;

      linestream2 >> word;

      G4int Lp;
      std::stringstream ss8(word);
      ss8 >> Lp;

      linestream2 >> word;

      G4double del;
      std::stringstream ss9(word);
      ss9 >> del;

      G4cout << "  " << index << " " << BR << " " << L0 << " " << Lp << " " << del << G4endl;
      
      if(!pConsidered || state_index == pConsidered) {
	part->GetDecayTable()->Insert(new Gamma_Decay(ppart,pLevels.at(index),BR,L0,Lp,del));
      }
      
    }

    if((unsigned int)state_index != pLevels.size()) {
      G4cout << "States are out of order in projectile level scheme file " << pFN << "!" << G4endl;
    }

    pLevels.push_back(ppart);
    
  }

  G4cout << pLevels.size()-1  << " excited states built for the projectile!" << G4endl;

  return;
  
}

void Excitation::BuildRecoilLS(G4int Z, G4int A) {
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition* recGS = table->GetIon(Z,A,0.0*MeV);
  recGS->SetPDGLifeTime(-1.0);
  
  Polarized_Particle* polGS = new Polarized_Particle(recGS,Z,A,rGSS,0.0*MeV);
  rLevels.push_back(polGS);

  if(rSimple) {

    G4cout << "\nBuilding simple recoil level scheme!\n";
    
    G4ParticleDefinition* part = table->GetIon(Z,A,rSimpleEn);
    part->SetPDGLifeTime(rSimpleLt);
    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,2.0,rSimpleEn);
    part->GetDecayTable()->Insert(new Gamma_Decay(ppart,rLevels.at(0),1,2,0,0));
    
    rLevels.push_back(ppart);

    G4cout << " 1 " << rSimpleEn/keV << " 2 " << rSimpleLt/ps << " 1\n  0 1 2 0 0\nSuccess!" << G4endl;  

    return;
    
  }

  if(rFN == "") {
    G4cout << "\nNo recoil excitations." << G4endl;
    return;
  }
  
  std::ifstream file;
  file.open(rFN.c_str(),std::ios::in);
  
  if(!file.is_open()) {
    G4cout << "\nCould not open recoil level scheme file " << rFN << "! No levels will be built!"
	   << G4endl;
    return;
  }

  G4cout << "\nBuilding recoil level scheme from " << rFN << G4endl;
  
  std::string line, word;
  while(std::getline(file,line)) {
      
    std::stringstream linestream1(line);
    linestream1 >> word;

    G4int state_index;
    std::stringstream ss(word);
    ss >> state_index;

    linestream1 >> word;

    G4double energy;
    std::stringstream ss1(word);
    ss1 >> energy;
    energy *= keV;
    
    linestream1 >> word;
    
    G4double spin;
    std::stringstream ss2(word);
    ss2 >> spin;

    linestream1 >> word;

    G4double lifetime;
    std::stringstream ss3(word);
    ss3 >> lifetime;
    lifetime *= ps;

    linestream1 >> word;

    G4int nbr;
    std::stringstream ss4(word);
    ss4 >> nbr;
    
    if(nbr == 0) {
      G4cout << "Problem reading recoil level scheme file " << rFN
	     << "! No decay braches declared for state " << state_index << G4endl;
    }
	
    G4ParticleDefinition* part = table->GetIon(Z,A,energy);
    if(!rConsidered || state_index == rConsidered) {
      
      part->SetPDGLifeTime(lifetime);
      part->SetDecayTable(new G4DecayTable());
      part->GetProcessManager()->SetParticleType(part);
      part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);
      
    }
    else {
      part->SetPDGLifeTime(-1.0);
    }

    G4cout << " " << state_index << " " << energy/keV << " " << spin << " " << lifetime/ps << " " << nbr
	   << G4endl;
    
    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,spin,energy);
    for(int i=0;i<nbr;i++) {

      std::getline(file,line);
      std::stringstream linestream2(line);
      
      linestream2 >> word;

      G4int index;
      std::stringstream ss5(word);
      ss5 >> index;

      linestream2 >> word;
      
      G4double BR;
      std::stringstream ss6(word);
      ss6 >> BR;

      linestream2 >> word;

      G4int L0;
      std::stringstream ss7(word);
      ss7 >> L0;

      linestream2 >> word;

      G4int Lp;
      std::stringstream ss8(word);
      ss8 >> Lp;

      linestream2 >> word;

      G4double del;
      std::stringstream ss9(word);
      ss9 >> del;

      G4cout << "  " << index << " " << BR << " " << L0 << " " << Lp << " " << del << G4endl;

      if(!rConsidered || state_index == rConsidered) {
	part->GetDecayTable()->Insert(new Gamma_Decay(ppart,rLevels.at(index),BR,L0,Lp,del));
      }
	
    }

    if((unsigned int)state_index != rLevels.size()) {
      G4cout << "States are out of order in recoil level scheme file " << rFN << "!" << G4endl;
    }

    rLevels.push_back(ppart);
    
  }

  G4cout << rLevels.size()-1  << " excited states built for the recoil!" << G4endl;

  return;
  
}

void Excitation::BuildProjectileSplines() {

  if(pSimple) {
    G4cout << "\nProjectile simple state will always be populated" << G4endl;
    return;
  }

  if(pSelected > -1) {
    G4cout << "\nProjectile state " << pSelected << " will always be populated" << G4endl;
    return;
  }

  if(pFN == "") {
    return;
  }
  
  std::ifstream file;
  file.open(pPF.c_str(),std::ios::in);

  if(!file.is_open()) {
    G4cout << "\nCould not open projectile excitation probability file " << pPF
	   << "! No probabilities will be built!" << G4endl;
    return;
  }

  G4cout << "\nBuilding projectile excitation probabilities from " << pPF << G4endl;

  std::vector<G4double> thetas;
  std::vector< std::vector<G4double> > probs;

  std::string line, word;

  int line_num = 0;
  while(std::getline(file,line)) {
    
    std::stringstream linestream1(line);
    
    int word_num = 0;
    while(linestream1 >> word) {

      G4double temp;
      std::stringstream ss(word);
      ss >> temp;

      if(!line_num && word_num) {
	std::vector<G4double> tmpvec;
	probs.push_back(tmpvec);
      }

      if(!word_num) {
	thetas.push_back(temp);
      }
      else {
	probs.at(word_num-1).push_back(temp);
      }
      
      word_num++;
    }
    line_num++;
  }
  
  if(pConsidered) {

    G4cout << " Renormalizing for projectile considered state " << pConsidered << "..." << G4endl;
    
    std::vector<G4double> tmp1 = probs.at(pConsidered);
    probs.clear();

    G4double max = 0;
    for(unsigned int i=0;i<tmp1.size();i++) {
      if(tmp1.at(i) > max) {
	max = tmp1.at(i);
      }
    }

    for(unsigned int i=0;i<tmp1.size();i++) {
      tmp1.at(i) /= max;
    }

    std::vector<G4double> tmp0;
    for(unsigned int i=0;i<tmp1.size();i++) {
      tmp0.push_back(1.0-tmp1.at(i));
    }

    probs.push_back(tmp0);
    probs.push_back(tmp1);
    
  }

  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<probs.size();i++) {

    if(probs.at(i).size() != size) {

      G4int state = i;
      if(pConsidered && i == 1) {
	state = pConsidered;
      }
      
      G4cout << "Projectile state " << state << " has " << probs.at(i).size()
	     << " probability splines points, while there are " << size
	     << " thetaCM spline points! This is likely a mistake and could break things." << G4endl; 
    }
    
    pSplines.push_back(new G4DataInterpolation(&thetas[0],&(probs.at(i))[0],size,0,0));
    
  }

  if(!pConsidered) {
    if(pSplines.size() == pLevels.size()) {
      G4cout << "All " << pSplines.size() << " projectile probability splines successfully built!"
	     << G4endl;
    }
    else {
      G4cout << "Created " << pSplines.size() << " projectile probability splines for " << pLevels.size()
	     << " levels! The simulation might not behave properly!" << G4endl;
    }
  }
  else {
    if(pSplines.size() == 2) {
      G4cout << "Both projectile probability splines successfully built!" << G4endl;
    }
    else {
      G4cout << "Created " << pSplines.size() << " projectile probability splines instead of 2!"; 
      G4cout << "The simulation might not behave properly!" << G4endl;
    }
  }

  return;
  
}

void Excitation::BuildRecoilSplines() {

  if(rSimple) {
    G4cout << "\nRecoil simple state will always be populated" << G4endl;
    return;
  }

  if(rSelected > -1) {
    G4cout << "\nRecoil state " << rSelected << " will always be populated" << G4endl;
    return;
  }

  if(rFN == "") {
    return;
  }

  std::ifstream file;
  file.open(rPF.c_str(),std::ios::in);

  if(!file.is_open()) {
    G4cout << "\nCould not open recoil excitation probability file " << rPF
	   << "! No probabilities will be built!" << G4endl;
    return;
  }

  G4cout << "\nBuilding recoil excitation probabilities from " << rPF << G4endl;

  std::vector<G4double> thetas;
  std::vector< std::vector<G4double> > probs;

  std::string line, word;

  int line_num = 0;
  while(std::getline(file,line)) {
    
    std::stringstream linestream1(line);
    
    int word_num = 0;
    while(linestream1 >> word) {

      G4double temp;
      std::stringstream ss(word);
      ss >> temp;

      if(!line_num && word_num) {
	std::vector<G4double> tmpvec;
	probs.push_back(tmpvec);
      }

      if(!word_num) {
	thetas.push_back(temp);
      }
      else {
	probs.at(word_num-1).push_back(temp);
      }
      
      word_num++;
    }
    line_num++;
  }

  if(rConsidered) {

    G4cout << " Renormalizing for recoil considered state " << rConsidered << "..." << G4endl; 
    
    std::vector<G4double> tmp1 = probs.at(rConsidered);
    probs.clear();

    G4double max = 0;
    for(unsigned int i=0;i<tmp1.size();i++) {
      if(tmp1.at(i) > max) {
	max = tmp1.at(i);
      }
    }

    for(unsigned int i=0;i<tmp1.size();i++) {
      tmp1.at(i) /= max;
    }

    std::vector<G4double> tmp0;
    for(unsigned int i=0;i<tmp1.size();i++) {
      tmp0.push_back(1.0-tmp1.at(i));
    }

    probs.push_back(tmp0);
    probs.push_back(tmp1);
    
  }

  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<probs.size();i++) {

    G4int state = i;
    if(rConsidered && i == 1) {
      state = rConsidered;
    }

    if(probs.at(i).size() != size) {
      G4cout << "Recoil state " << state << " has " << probs.at(i).size()
	     << " probability splines points, while there are " << size
	     << " thetaCM spline points! This is likely a mistake and could break things." << G4endl; 
    }
    
    rSplines.push_back(new G4DataInterpolation(&thetas[0],&(probs.at(i))[0],size,0,0));
    
  }

  if(!rConsidered) {
    if(rSplines.size() == rLevels.size()) {
      G4cout << "All " << rSplines.size() << " recoil probability splines successfully built!" << G4endl;
    }
    else {
      G4cout << "Created " << rSplines.size() << " recoil probability splines for " << rLevels.size()
	     << " levels! The simulation might not behave properly!" << G4endl;
    }
  }
  else {
    if(rSplines.size() == 2) {
      G4cout << "Both recoil probability splines successfully built!" << G4endl;
    }
    else {
      G4cout << "Created " << rSplines.size() << " recoil probability splines instead of 2!"; 
      G4cout << "The simulation might not behave properly!" << G4endl;
    }
  }

  return;
  
}

G4int Excitation::ChooseProjectileState(const G4double th) {

  if(pSimple) {
    return 1;
  }

  if(pSelected > -1) {
    return pSelected;
  }

  if(!pSplines.size()) {
    return 0;
  }
  
  G4double sum = 0.0;
  std::vector<double> probs;

  for(unsigned int i=0;i<pSplines.size();i++) {
    
    G4double val = pSplines.at(i)->CubicSplineInterpolation(th);
    
    sum += val;
    probs.push_back(val);
    
  }
  
  for(unsigned int i=0;i<probs.size();i++) {
    probs.at(i) /= sum;
  }

  G4double sumBR = 0.0;
  G4double num = G4UniformRand();
  G4int index = 0;

  for(unsigned int i=0;i<probs.size();i++) {
    sumBR += probs.at(i);
    if(num < sumBR) {
      index = i;
      break;
    }
  }

  if(pConsidered && index) {
    return pConsidered;
  }
  
  return index;
  
}

G4int Excitation::ChooseRecoilState(const G4double th) {

  if(rSimple) {
    return 1;
  }

  if(rSelected > -1) {
    return rSelected;
  }
  
  if(!rSplines.size()) {
    return 0;
  }
  
  G4double sum = 0.0;
  std::vector<double> probs;

  for(unsigned int i=0;i<rSplines.size();i++) {
    
    G4double val = rSplines.at(i)->CubicSplineInterpolation(th);
    
    sum += val;
    probs.push_back(val);
    
  }
  
  for(unsigned int i=0;i<probs.size();i++) {
    probs.at(i) /= sum;
  }

  G4double sumBR = 0.0;
  G4double num = G4UniformRand();
  G4int index = 0;

  for(unsigned int i=0;i<probs.size();i++) {
    sumBR += probs.at(i);
    if(num < sumBR) {
      index = i;
      break;
    }
  }

  if(rConsidered && index) {
    return rConsidered;
  }
  
  return index;

}

G4double Excitation::GetProjectileExcitation(G4int index) {

  if(pSimple) {
    return pSimpleEn;
  }
  
  if(!index) {
    return 0.0*MeV;
  }

  return pLevels.at(index)->GetDefinition()->GetPDGMass() - pLevels.at(0)->GetDefinition()->GetPDGMass();
  
}

G4double Excitation::GetRecoilExcitation(G4int index) {

  if(rSimple) {
    return rSimpleEn;
  }

  if(!index) {
    return 0.0*MeV;
  }

  return rLevels.at(index)->GetDefinition()->GetPDGMass() - rLevels.at(0)->GetDefinition()->GetPDGMass();
  
}

void Excitation::Polarize(G4int pI, G4int rI, G4double th, G4double ph) {

  if(pI) {
    pLevels.at(pI)->SetPolarization(polar->GetProjectilePolarization(pI,th,ph));
  }

  if(rI) {
    rLevels.at(rI)->SetPolarization(polar->GetRecoilPolarization(rI,th,ph));
  }
  
  return;
}

void Excitation::Unpolarize() {

  for(auto p : pLevels) {
    p->Unpolarize();
  }

  for(auto r : rLevels) {
    r->Unpolarize();
  }
  
  return;
}
