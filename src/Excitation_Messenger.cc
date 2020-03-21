#include "Excitation_Messenger.hh"

Excitation_Messenger::Excitation_Messenger(Excitation* exc) : excitation(exc) {

  //All information about level schemes and excitation probabilities
  excitation_dir = new G4UIdirectory("/Excitation/");

  //Projectile directory
  proj_dir = new G4UIdirectory("/Excitation/Projectile/");

  //Level scheme file
  pLS_cmd = new G4UIcmdWithAString("/Excitation/Projectile/LevelScheme",this);
  pLS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pLS_cmd->SetGuidance("Set name of projectile level scheme file");

  //Probabilities file
  pPF_cmd = new G4UIcmdWithAString("/Excitation/Projectile/Probabilities",this);
  pPF_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pPF_cmd->SetGuidance("Set name of projectile excitation probabilities file");

  //Simple excitation
  pSim_cmd = new G4UIcmdWithoutParameter("/Excitation/Projectile/Simple",this);
  pSim_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pSim_cmd->SetGuidance("Declare a simple excitation scheme for the projectile");

  //Simple state energy
  pSEn_cmd = new G4UIcmdWithADoubleAndUnit("/Excitation/Projectile/SimpleEnergy",this);
  pSEn_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pSEn_cmd->SetGuidance("Set energy of simple state in the projectile");

  //Simple state lifetime
  pSLt_cmd = new G4UIcmdWithADoubleAndUnit("/Excitation/Projectile/SimpleLifetime",this);
  pSLt_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pSLt_cmd->SetGuidance("Set lifetime of simple state in the projectile");

  pSel_cmd = new G4UIcmdWithAnInteger("/Excitation/Projectile/PopulateState",this);
  pSel_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pSel_cmd->SetGuidance("Choose one state to always populate in the projectile");
  
  //Recoil directory
  rec_dir = new G4UIdirectory("/Excitation/Recoil/");

  //Level scheme file
  rLS_cmd = new G4UIcmdWithAString("/Excitation/Recoil/LevelScheme",this);
  rLS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rLS_cmd->SetGuidance("Set name of recoil level scheme file");

  //Probabilities file
  rPF_cmd = new G4UIcmdWithAString("/Excitation/Recoil/Probabilites",this);
  rPF_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rPF_cmd->SetGuidance("Set name of recoil excitation probabilities file");

  //Simple excitation
  rSim_cmd = new G4UIcmdWithoutParameter("/Excitation/Recoil/Simple",this);
  rSim_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rSim_cmd->SetGuidance("Declare a simple excitation scheme for the recoil");

  //Simple state energy
  rSEn_cmd = new G4UIcmdWithADoubleAndUnit("/Excitation/Recoil/SimpleEnergy",this);
  rSEn_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rSEn_cmd->SetGuidance("Set energy of simple state in the recoil");

  //Simple state lifetime
  rSLt_cmd = new G4UIcmdWithADoubleAndUnit("/Excitation/Recoil/SimpleLifetime",this);
  rSLt_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rSLt_cmd->SetGuidance("Set lifetime of simple state in the recoil");

  rSel_cmd = new G4UIcmdWithAnInteger("/Excitation/Recoil/PopulateState",this);
  rSel_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rSel_cmd->SetGuidance("Choose one state to always populate in the recoil");
  
}

Excitation_Messenger::~Excitation_Messenger() {

  delete excitation_dir;
  delete proj_dir;
  delete rec_dir;
  
  delete pLS_cmd;
  delete pPF_cmd;
  
  delete pSim_cmd;
  delete pSEn_cmd;
  delete pSLt_cmd;
  delete pSel_cmd;
  
  delete rLS_cmd;
  delete rPF_cmd;

  delete rSim_cmd;
  delete rSEn_cmd;
  delete rSLt_cmd;
  delete rSel_cmd;
  
}

void Excitation_Messenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  //////////Projectile commands//////////
  if(command == pLS_cmd) {
    excitation->SetProjLSFile(newValue);
    G4cout << "Setting projectile level scheme file to " << newValue << G4endl;
  }

  else if(command == pPF_cmd) {
    excitation->SetProjPrbFile(newValue);
    G4cout << "Setting projectile excitation probabilites file to " << newValue << G4endl;
  }

  else if(command == pSim_cmd) {
    excitation->SetSimpleProj();
    G4cout << "Setting a simple excitation scheme for the projectile!" << G4endl;
  }

  else if(command == pSEn_cmd) {
    excitation->SetSimpleProjEn(pSEn_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting energy of simple projectile state to " << newValue << G4endl;
  }

  else if(command == pSLt_cmd) {
    excitation->SetSimpleProjLt(pSLt_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting lifetime of simple projectile state to " << newValue << G4endl;
  }

  else if(command == pSel_cmd) {
    excitation->FixProjState(pSel_cmd->GetNewIntValue(newValue));
    G4cout << "Selecting state " << newValue << " to populate in the projectile" << G4endl;
  }
  ///////////////////////////////////////

  ////////////Recoil commands////////////
  else if(command == rLS_cmd) {
    excitation->SetRecLSFile(newValue);
    G4cout << "Setting recoil level scheme file to " << newValue << G4endl;
  }

  else if(command == rPF_cmd) {
    excitation->SetRecPrbFile(newValue);
    G4cout << "Setting recoil excitation probabilites file to " << newValue << G4endl;
  }

  else if(command == rSim_cmd) {
    excitation->SetSimpleRec();
    G4cout << "Setting a simple excitation scheme for the recoil!" << G4endl;
  }

  else if(command == rSEn_cmd) {
    excitation->SetSimpleRecEn(rSEn_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting energy of simple recoil state to " << newValue << G4endl;
  }

  else if(command == rSLt_cmd) {
    excitation->SetSimpleRecLt(rSLt_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting lifetime of simple recoil state to " << newValue << G4endl;
  }

  else if(command == rSel_cmd) {
    excitation->FixRecState(rSel_cmd->GetNewIntValue(newValue));
    G4cout << "Selecting state " << newValue << " to populate in the recoil" << G4endl;
  }
  ///////////////////////////////////////

  return;
  
}
