#include "Gamma_Source.hh"
#include "Gamma_Decay.hh"

#include "G4ProcessManager.hh"
#include "G4IonTable.hh"
#include "G4Decay.hh"
#include "G4DecayTable.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

Gamma_Source::Gamma_Source() {

  messenger = new Gamma_Source_Messenger(this);

  source_energy = -1.0*MeV;
  GSS = 0.0;
  file_name = "";
  
}

Gamma_Source::~Gamma_Source() {

  delete messenger;
  
}

void Gamma_Source::BuildLevelScheme() {

  if(source_energy > 0.0*MeV) {
    G4cout << "\nSimple isotropic gamma-ray of " << source_energy/keV << " keV will be emitted each event"
	   << G4endl;

    return;
  }

  if(file_name == "") {
    G4cout << "Neither the source level scheme file nor the source energy is set! The simulation will break..." << G4endl;

    return;
  }

  G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());

  G4ParticleDefinition* ground_state = table->GetIon(82,208,0.0*MeV);
  ground_state->SetPDGLifeTime(-1.0);
  
  Polarized_Particle* polGS = new Polarized_Particle(ground_state,82,208,GSS,0.0*MeV);
  levels.push_back(polGS);
 
  std::ifstream file;
  file.open(file_name.c_str(),std::ios::in);
  
  if(!file.is_open()) {
    G4cout << "\nCould not open source level scheme file " << file_name << "! No levels will be built!"
	   << G4endl;
    return;
  }

  G4cout << "\nBuilding source level scheme from " << file_name << G4endl;
  
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

    G4double prob;
    std::stringstream ss4(word);
    ss4 >> prob;
    probs.push_back(prob);

    linestream1 >> word;

    G4int nbr;
    std::stringstream ss11(word);
    ss11 >> nbr;
    
    if(nbr == 0) {
      G4cout << "Problem reading source level scheme file " << file_name
	     << "! No decay branches declared for state " << state_index << G4endl;
    }
	
    G4ParticleDefinition* part = table->GetIon(82,208,energy);
    part->SetPDGLifeTime(lifetime);
    part->SetDecayTable(new G4DecayTable());
    part->GetProcessManager()->SetParticleType(part);
    part->GetProcessManager()->AddProcess(new G4Decay(),0,-1,0);

    G4cout << " " << state_index << " " << energy/keV << " " << spin << " " << lifetime/ps << " "
	   << prob << " " << nbr << G4endl;
    
    Polarized_Particle* ppart = new Polarized_Particle(part,82,208,spin,energy);
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

      linestream2 >> word;
      
      G4double cc;
      std::stringstream ss10(word);
      ss10 >> cc;

      G4cout << "  " << index << " " << BR << " " << L0 << " " << Lp << " " << del << " " << cc << G4endl;
      
      part->GetDecayTable()->Insert(new Gamma_Decay(ppart,levels.at(index),BR,L0,Lp,del,cc));
      
    }

    if((unsigned int)state_index != levels.size()) {
      G4cout << "States are out of order in source level scheme file " << file_name << "!" << G4endl;
    }

    levels.push_back(ppart);
    
  }

  if(levels.size()-1 == probs.size()) {
    G4cout << levels.size()-1  << " excited states built for the source!" << G4endl;
  }
  else {
    G4cout << "There are " << levels.size()-1 << " excited states but " << probs.size()
	   << " population probabilities!" << G4endl;
  }

  G4double sum = 0.0;
  for(unsigned int i=0;i<probs.size();i++) {
    sum += probs.at(i);
  }

  for(unsigned int i=0;i<probs.size();i++) {
    probs.at(i) /= sum;
  }
  
  return;
  
}

G4int Gamma_Source::ChooseState() {

  G4double sumBR = 0.0;
  G4double num = G4UniformRand();

  for(unsigned int i=0;i<probs.size();i++) {
    sumBR += probs.at(i);
    if(num < sumBR) {
      return i+1;
    }
  }

  return 0;
}

void Gamma_Source::Unpolarize() {

  for(auto lvl : levels) {
    lvl->Unpolarize();
  }
  
  return;
}
