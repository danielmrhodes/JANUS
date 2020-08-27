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
  
  pFN = "";
  pPF = "";

  pSimple = false;
  pSimpleEn = 0.0*MeV;;
  pSimpleLt = 0.0*ps;
  pSelected = -1;
    
  rFN = "";
  rPF = "";

  rSimple = false;
  rSimpleEn = 0.0*MeV;;
  rSimpleLt = 0.0*ps;
  rSelected = -1;

  //popGS = false;
  
}

Excitation::~Excitation() {

  delete messenger;
}

void Excitation::BuildLevelSchemes(int pZ, int pA, int rZ, int rA) {

  BuildProjectileLS(pZ,pA);
  BuildRecoilLS(rZ,rA);
  
  return;
}

void Excitation::BuildProbabilities() {

  BuildProjectileSplines();
  BuildRecoilSplines();
  
  return;
}

void Excitation::BuildProjectileLS(int Z, int A) {

  G4int J = 2;
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition* projGS = table->GetIon(Z,A,0.0*MeV);
  projGS->SetPDGStable(true);

  Polarized_Particle* polGS = new Polarized_Particle(projGS,Z,A,0,0.0*MeV);
  pLevels.push_back(polGS);

  if(pSimple) {

    G4cout << "\nBuilding simple projectile level scheme!\n";
    
    G4ParticleDefinition* part = table->GetIon(Z,A,pSimpleEn);
    part->SetPDGStable(false);
    part->SetPDGLifeTime(pSimpleLt);

    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,J,pSimpleEn);
    part->GetDecayTable()->Insert(new Gamma_Decay(ppart,pLevels.at(0),1));

    pLevels.push_back(ppart);

    G4cout << " 1 " << pSimpleEn/keV << " " << pSimpleLt/ps << " 1\n  0 1\nSuccess!" << G4endl;  

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

  unsigned int state_index = 0;
  G4double energy = 1.*MeV;
  G4double lifetime = 1.*ps;
  
  std::string line, word;
  while(std::getline(file,line)) {
    
    G4ParticleDefinition* part;
    int nbr = 0;
    
    std::stringstream linestream1(line);
    G4int word_num = 0;
    while(linestream1 >> word) {

      G4double temp;
      std::stringstream ss(word);
      ss >> temp;

      switch (word_num) {

        case 0: { //Index
	  state_index = (int)temp;
	  break;
        }

        case 1: { //State energy
	  energy = temp*keV;
	  part = table->GetIon(Z,A,energy);
	  break;
        }

        case 2: { //State lifetime
	  lifetime = temp*ps;
	  part->SetPDGStable(false);
	  part->SetPDGLifeTime(lifetime);
	  break;
        }

        case 3: { //Number of branches
	  nbr = (int)temp;
	  if(nbr == 0) {
	    G4cout << "Probem reading projectile level scheme file " << pFN
		   << "! No decay braches declared for state " << state_index << G4endl;
	  }
	  break;
        }

        default: {
	  G4cout << "Probem reading projectile level scheme file " << pFN
		 << "! Too many entries for state " <<
		 state_index << G4endl;
	  break;
        }
	
      }

      word_num++;
    }

    G4cout << " " << state_index << " " << energy/keV << " " << lifetime/ps << " " << nbr << G4endl;

    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,J,energy);
    for(int i=0;i<nbr;i++) {

      std::getline(file,line);
      std::stringstream linestream2(line);
      
      linestream2 >> word;

      int index;
      std::stringstream ss1(word);
      ss1 >> index;

      linestream2 >> word;
      
      double BR;
      std::stringstream ss2(word);
      ss2 >> BR;

      G4cout << "  " << index << " " << BR << G4endl;

      part->GetDecayTable()->Insert(new Gamma_Decay(ppart,pLevels.at(index),BR));
	
    }

    if(state_index != pLevels.size()) {
      G4cout << "States are out of order in projectile level scheme file " << pFN << "!" << G4endl;
    }

    pLevels.push_back(ppart);
    
  }

  G4cout << pLevels.size()-1  << " excited states built for the projectile!" << G4endl;

  return;
  
}

void Excitation::BuildRecoilLS(int Z, int A) {

  G4int J = 2;
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition* recGS = table->GetIon(Z,A,0.0*MeV);
  recGS->SetPDGStable(true);

  Polarized_Particle* polGS = new Polarized_Particle(recGS,Z,A,0,0.0*MeV);
  rLevels.push_back(polGS);

  if(rSimple) {

    G4cout << "\nBuilding simple recoil level scheme!\n";
    
    G4ParticleDefinition* part = table->GetIon(Z,A,rSimpleEn);
    part->SetPDGStable(false);
    part->SetPDGLifeTime(rSimpleLt);

    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,J,rSimpleEn);
    part->GetDecayTable()->Insert(new Gamma_Decay(ppart,rLevels.at(0),1));
    
    rLevels.push_back(ppart);

    G4cout << " 1 " << rSimpleEn/keV << " " << rSimpleLt/ps << " 1\n  0 1\nSuccess!" << G4endl;  

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

  unsigned int state_index = 0;
  G4double energy = 1.*MeV;
  G4double lifetime = 1.*ps;
  
  std::string line, word;
  while(std::getline(file,line)) {

    G4ParticleDefinition* part;
    int nbr = 0;
    
    std::stringstream linestream1(line);
    G4int word_num = 0;
    while(linestream1 >> word) {

      G4double temp;
      std::stringstream ss(word);
      ss >> temp;

      switch (word_num) {

        case 0: { //Index
	  state_index = (int)temp;
	  break;
        }

        case 1: { //State energy
	  energy = temp*keV;
	  part = table->GetIon(Z,A,energy);
	  break;
        }

        case 2: { //State lifetime
	  lifetime = temp*ps;
	  part->SetPDGStable(false);
	  part->SetPDGLifeTime(lifetime);
	  break;
        }

        case 3: { //Number of branches
	  nbr = (int)temp;
	  if(nbr == 0) {
	    G4cout << "Probem reading recoil level scheme file " << rFN
		   << "! No decay braches declared for state " << state_index << G4endl;
	  }
	  break;
        }

        default: {
	  G4cout << "Probem reading recoil level scheme file " << rFN
		 << "! Too many entries for state " <<
		 state_index << G4endl;
	  break;
        }
	
      }

      word_num++;
    }

    G4cout << " " << state_index << " " << energy/keV << " " << lifetime/ps << " " << nbr << G4endl;

    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,J,energy);
    for(int i=0;i<nbr;i++) {

      std::getline(file,line);
      std::stringstream linestream2(line);
      
      linestream2 >> word;

      int index;
      std::stringstream ss1(word);
      ss1 >> index;

      linestream2 >> word;
      
      double BR;
      std::stringstream ss2(word);
      ss2 >> BR;

      G4cout << "  " << index << " " << BR << G4endl;

      part->GetDecayTable()->Insert(new Gamma_Decay(ppart,rLevels.at(index),BR));
	
    }

    if(state_index != rLevels.size()) {
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
  std::vector<std::vector<G4double>> probs;

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

  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<probs.size();i++) {

    if(probs.at(i).size() != size) {
      G4cout << "Projectile state " << i << " has " << probs.at(i).size()
	     << " probability splines points, while there are " << size
	     << " thetaCM spline points! This is likely a mistake and could break things." << G4endl; 
    }
    
    pSplines.push_back(new G4DataInterpolation(&thetas[0],&(probs.at(i))[0],size,0,0));
    
  }

  if(pSplines.size() == pLevels.size()) {
    G4cout << "All " << pSplines.size() << " projectile probability splines successfully built!" << G4endl;
  }
  else {
    G4cout << "Created " << pSplines.size() << " projectile probability splines for " << pLevels.size()
	   << " levels! The simulation might not behave properly!" << G4endl;
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
  std::vector<std::vector<G4double>> probs;

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

  const unsigned int size = thetas.size();
  for(unsigned int i=0;i<probs.size();i++) {

    if(probs.at(i).size() != size) {
      G4cout << "Recoil state " << i << " has " << probs.at(i).size()
	     << " probability splines points, while there are " << size
	     << " thetaCM spline points! This is likely a mistake and could break things." << G4endl; 
    }
    
    rSplines.push_back(new G4DataInterpolation(&thetas[0],&(probs.at(i))[0],size,0,0));
    
  }

  if(rSplines.size() == rLevels.size()) {
    G4cout << "All " << rSplines.size() << " recoil probability splines successfully built!" << G4endl;
  }
  else {
    G4cout << "Created " << rSplines.size() << " recoil probability splines for " << rLevels.size()
	   << " levels! The simulation might not behave properly!" << G4endl;
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

  for(unsigned int i=0;i<probs.size();i++) {
    sumBR += probs.at(i);
    if(num < sumBR) {
      return i;
    }
  }

  G4cout << "No decay channel selected for the projectile! Defaulting to the ground state." << G4endl;
  
  return 0;
  
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

  for(unsigned int i=0;i<probs.size();i++) {
    sumBR += probs.at(i);
    if(num < sumBR) {
      return i;
    }
  }

  G4cout << "No decay channel selected for the recoil! Defaulting to the ground state." << G4endl;
  
  return 0;
}

G4double Excitation::GetProjectileExcitation(int index) {

  if(pSimple) {
    return pSimpleEn;
  }
  
  if(!index) {
    return 0.0*MeV;
  }

  return pLevels.at(index)->GetDefinition()->GetPDGMass() - pLevels.at(0)->GetDefinition()->GetPDGMass();
  
}

G4double Excitation::GetRecoilExcitation(int index) {

  if(rSimple) {
    return rSimpleEn;
  }

  if(!index) {
    return 0.0*MeV;
  }

  return rLevels.at(index)->GetDefinition()->GetPDGMass() - rLevels.at(0)->GetDefinition()->GetPDGMass();
  
}
