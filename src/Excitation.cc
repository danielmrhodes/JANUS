#include "Excitation.hh"
#include "Gamma_Decay.hh"
#include "Primary_Generator.hh"

#include "G4RunManager.hh"
#include "G4ProcessManager.hh"
#include "G4IonTable.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

Excitation::Excitation(G4bool prj) : proj(prj) {

  messenger = new Excitation_Messenger(this, proj);
  polar = new Polarization(proj);
  
  lfn = "";
  pfn = "";
  
  selected = -1;
  considered = 0;
  gss = 0.0;
  
  xacc = gsl_interp_accel_alloc();
  yacc = gsl_interp_accel_alloc();
  
}

Excitation::~Excitation() {

  delete messenger;
  delete polar;

  for(unsigned int i=0;i<interps.size();i++)
    gsl_spline2d_free(interps.at(i));
  
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  
}

void Excitation::BuildStatisticalTensors() {

  if(lfn == "")
    return;
  
  polar->BuildStatisticalTensors();

  return;
}

void Excitation::BuildLevelScheme() {

  G4String nuc;
  if(proj)
    nuc = "projectile";
  else
    nuc = "recoil";

  G4RunManager* Rman = G4RunManager::GetRunManager();
  Primary_Generator* gen = (Primary_Generator*)Rman->GetUserPrimaryGeneratorAction();
  
  G4int Z = gen->GetZ(proj);
  G4int A = gen->GetA(proj);
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4ParticleDefinition* GS = table->GetIon(Z,A,0.0*MeV);
  GS->SetPDGLifeTime(-1.0);
  
  Polarized_Particle* polGS = new Polarized_Particle(GS,Z,A,gss,0.0*MeV);
  levels.push_back(polGS);

  if(lfn == "") {
    G4cout << "\nNo " << nuc << " excitations." << G4endl;
    return;
  }
  G4cout << "\nBuilding " << nuc << " level scheme from " << lfn << G4endl;
  
  ReadLevelSchemeFile(Z,A);
  G4cout << levels.size()-1  << " excited states built for the " << nuc << "!" << G4endl;

  return;
  
}

void Excitation::ReadLevelSchemeFile(G4int Z, G4int A) {

  G4String nuc;
  if(proj)
    nuc = "projectile";
  else
    nuc = "recoil";
  
  std::ifstream file;
  file.open(lfn.c_str(),std::ios::in);
  
  if(!file.is_open()) {
    G4cout << " \033[1;31mCould not open " << nuc << " level scheme file " << lfn
	   << "! No levels will be built! \033[m" << G4endl;
    return;
  }
  
  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  
  std::string line, word;
  while(std::getline(file,line)) {

    G4int state_index, nbr;
    G4double energy, spin, lifetime; 
    
    std::stringstream ss1(line);
    ss1 >> state_index >> energy >> spin >> lifetime >> nbr;

    energy *= keV;
    lifetime *= ps;

    G4cout << " " << state_index << " " << energy/keV << " " << spin << " " << lifetime/ps
	   << " " << nbr;
    
    G4ParticleDefinition* part = table->GetIon(Z,A,energy);
    if((!considered || state_index == considered) && nbr) {
      part->SetPDGLifeTime(lifetime);
      part->SetDecayTable(new G4DecayTable());
      part->GetProcessManager()->SetParticleType(part);
      part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);
    }
    if(!nbr) {
      G4cout << " \033[1;36m Warning: " << nuc << " state " << state_index
	     << " has no decay branches.\033[m";
      part->SetPDGLifeTime(-1.0);
    }
    G4cout << G4endl;
    
    Polarized_Particle* ppart = new Polarized_Particle(part,Z,A,spin,energy);
    for(int i=0;i<nbr;i++) {

      std::getline(file,line);
      std::stringstream ss2(line);

      G4int index, L0, Lp;
      G4double BR, del, cc;
      ss2 >> index >> BR >> L0 >> Lp >> del >> cc;

      G4cout << "  " << index << " " << BR << " " << L0 << " " << Lp << " " << del << " " << cc
	     << G4endl;
      
      if(!considered || state_index == considered)
	part->GetDecayTable()->Insert(new Gamma_Decay(ppart,levels.at(index),BR,L0,Lp,del,cc));
      
    }

    if((unsigned int)state_index != levels.size())
      G4cout << "\033[1;31m" << nuc << "states are out of order in level scheme file " << lfn
	     << "\033[m" << G4endl;

    levels.push_back(ppart);
    
  }
  
  return;
}

void Excitation::BuildProbabilities() {
  
  if(lfn == "")
    return;

  G4String nuc;
  if(proj)
    nuc = "projectile";
  else
    nuc = "recoil";
  
  if(selected > -1) {
    G4cout << "\n" << nuc << " state " << selected << " will always be populated" << G4endl;
    return;
  }
  
  G4cout << "\nBuilding " << nuc << " excitation probabilities from " << pfn << G4endl;
  ReadProbFile();
  
  if(considered) {
    G4cout << " Renormalizing for " << nuc << " considered state " << considered << "..." << G4endl;
    Renormalize();    
  }
    
  G4int numE = energies.size();
  G4int numT = thetas.size();
  G4int numS = levels.size();
  G4int numP = probs.size();

  for(G4int i=0;i<numS;i++) {
    interps.push_back(gsl_spline2d_alloc(gsl_interp2d_bicubic,numE,numT));  
    gsl_spline2d_init(interps.back(),&energies[0],&thetas[0],&probs[i*numE*numT],numE,numT);
  }
  
  if(!considered)
    if(numP == numS*numE*numT)
      G4cout << "All " << nuc << " excitation probability grids successfully built!" << G4endl;
    else
      G4cout << "\033[1;31mThere was a problem building the excitation grids for the " << nuc
	     << "!\033[m" << G4endl;
  else
    if(numP == 2*numE*numT)
      G4cout << "Both " << nuc << " excitation probability grids successfully built!" << G4endl;
    else
      G4cout << "\033[1;31mThere was a problem building the excitation grids for the " << nuc
	     << "!\033[m" << G4endl;
  
  return;
  
}

void Excitation::ReadProbFile() {

  energies.clear();
  thetas.clear();
  probs.clear();
  
  std::ifstream file;
  file.open(pfn.c_str(),std::ios::in);

  G4String nuc;
  if(proj)
    nuc = "projectile";
  else
    nuc = "recoil";
  
  if(!file.is_open()) {
    G4cout << "\n\033[1;31mCould not open " << nuc << " excitation probability file " << pfn
	   << "! No probabilities will be built!\033[m" << G4endl;
    return;
  }

  std::string line;
  std::getline(file,line);

  G4double en;
  std::stringstream ss1(line);

  while(ss1 >> en)
    energies.push_back(en);
  G4int numE = energies.size();

  std::getline(file,line);

  G4double th;
  std::stringstream ss2(line);

  while(ss2 >> th)
    thetas.push_back(th);
  G4int numT = thetas.size();

  std::getline(file,line);
  
  G4int indexE = 0;
  G4int indexT = 0;
  G4int indexS = 0;
  
  while(std::getline(file,line)) {

    if(line.empty()) {
      indexS = 0;
      indexE++;
      continue;
    }

    G4int numP = probs.size();
    if(numP == indexS*numE*numT)
      probs.resize(numP + numE*numT);
    
    G4double prb;
    std::stringstream ss3(line);
    while(ss3 >> prb) {
      probs.at(indexS*numE*numT + indexT*numE + indexE) = prb;
      indexT++;
    }

    indexT = 0;
    indexS++;
  }
  
  return;
}

void Excitation::Renormalize() {

  G4String nuc;
  if(proj)
    nuc = "projectile";
  else
    nuc = "recoil";
  
  G4int numE = energies.size();
  G4int numT = thetas.size();
  
  std::vector<G4double> tmp1;
  tmp1.resize(numE*numT);
  for(int i=0;i<numE;i++) {
    for(int j=0;j<numT;j++) {
      tmp1.at(i*numT + j) = probs.at(considered*numE*numT + i*numT + j);
    }
  }
  probs.clear();
    
  G4double max = 0.0;
  for(auto p : tmp1)
    if(p > max)
      max = p;

  for(unsigned int i=0;i<tmp1.size();i++)
    tmp1.at(i) /= max;
  
  G4cout << nuc << " state " << considered << " excitation probabilites rescaled by " << 1.0/max
	 << G4endl;
  
  std::vector<G4double> tmp0;
  for(auto p : tmp1)
    tmp0.push_back(1.0-p);

  for(auto p : tmp0)
    probs.push_back(p);

  for(auto p : tmp1)
    probs.push_back(p);
  
  return;
}

G4int Excitation::ChooseState(G4double en, G4double th) {
  
  if(selected > -1)
    return selected;

  if(!probs.size())
    return 0;
  
  G4int numE = energies.size();
  G4int numT = thetas.size();
  G4int numS = levels.size();

  if(en < energies[0]) {

    static G4int counter = 0;
    if(counter < 5)
      G4cout << "\033[1;36mYou need to go to a lower energy when making the "
	     << "probability input file (E = " << en << ", min = " << energies[0] << ")\033[m"
	     << G4endl;
    counter++;
    
    en = energies[0];
  }
  else if(en > energies[numE-1]) {

    static G4int counter = 0;
    if(counter < 5)
      G4cout << "\033[1;36mYou need to go to a higher energy when making the "
	     << "probability input file (E = " << en << ", max = " << energies[numE-1] << ")\033[m"
	     << G4endl;
    counter++;
    
    en = energies[numE-1];
  }

  if(th < thetas[0]) {

    static G4int counter = 0;
    if(counter < 5)
      G4cout << "\033[1;36mYou need to go to a lower theta when making the "
	     << "probability input file (th = " << th << ", min = " << thetas[0] << ")\033[m"
	     << G4endl;
    counter++;
    
    th = thetas[0]; 
  }
  else if(th > thetas[numT-1]) {

    static G4int counter = 0;
    if(counter < 5)
      G4cout << "\033[1;36mYou need to go to a higher theta when making the "
	     << "probability input file (th = " << th << ", max = " << thetas[numT-1] << ")\033[m"
	     << G4endl;
    counter++;
    
    th = thetas[numT-1];
  }

  G4double sum = 0.0;
  std::vector<G4double> prbs;
  prbs.resize(numS);
  
  for(G4int i=0;i<numS;i++) {
    
    G4double val = gsl_spline2d_eval(interps[i],en,th,xacc,yacc);
    sum += val;
    prbs[i] = val;
    
  }

  for(G4int i=0;i<numS;i++)
    prbs[i] /= sum;

  G4double sumBR = 0.0;
  G4double num = G4UniformRand();
  G4int index = 0;

  for(G4int i=0;i<numS;i++) {
    sumBR += prbs[i];
    if(num < sumBR) {
      index = i;
      break;
    }
  }

  if(considered && index)
    return considered;
  
  return index;
  
}


G4double Excitation::GetExcitation(G4int index) {
  
  if(!index)
    return 0.0*MeV;

  return levels[index]->GetDefinition()->GetPDGMass() - levels[0]->GetDefinition()->GetPDGMass();
  
}

void Excitation::Polarize(G4int index, G4double en, G4double th, G4double ph) {

  if(index)
    levels[index]->SetPolarization(polar->GetPolarization(index,en,th,ph));
  
  return;
}

void Excitation::Unpolarize() {

  for(auto& lvl : levels)
    lvl->Unpolarize();
  
  return;
}
