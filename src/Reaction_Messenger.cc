#include "Reaction_Messenger.hh"

Reaction_Messenger::Reaction_Messenger(Reaction* reac) : reaction(reac) {

  //All info about the scattering process accessed through this directory
  reaction_dir = new G4UIdirectory("/Reaction/");

  //Projectile
  beamZ_cmd = new G4UIcmdWithAnInteger("/Reaction/ProjectileZ",this);
  beamZ_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  beamZ_cmd->SetGuidance("Set Z of projectile nucleus");

  beamA_cmd = new G4UIcmdWithAnInteger("/Reaction/ProjectileA",this);
  beamA_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  beamA_cmd->SetGuidance("Set A of projectile nucleus");

  onlyP_cmd = new G4UIcmdWithoutParameter("/Reaction/OnlyProjectiles",this);
  onlyP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  onlyP_cmd->SetGuidance("Only consider the projectile when defining desired LAB scattering angle ranges");

  //Recoil
  recoilZ_cmd = new G4UIcmdWithAnInteger("/Reaction/RecoilZ",this);
  recoilZ_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  recoilZ_cmd->SetGuidance("Set Z of recoil nucleus");

  recoilA_cmd = new G4UIcmdWithAnInteger("/Reaction/RecoilA",this);
  recoilA_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  recoilA_cmd->SetGuidance("Set A of recoil nucleus");

  onlyR_cmd = new G4UIcmdWithoutParameter("/Reaction/OnlyRecoils",this);
  onlyR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  onlyR_cmd->SetGuidance("Only consider the recoil when defining desired LAB scattering angle ranges");

  //Scattering angle commands
  toJanus_cmd = new G4UIcmdWithoutParameter("/Reaction/SendToJanus",this);
  toJanus_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  toJanus_cmd->SetGuidance("Only sample parts of the Rutherford distribution which result in a particle entering Bambino2");

  toUS_cmd = new G4UIcmdWithoutParameter("/Reaction/SendToUpstreamJanus",this);
  toUS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  toUS_cmd->SetGuidance("Only sample parts of the Rutherford distribution which result in the projectile entering the upstream Bambino2 detector");

  toDS_cmd = new G4UIcmdWithoutParameter("/Reaction/SendToDownstreamJanus",this);
  toDS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  toDS_cmd->SetGuidance("Only sample parts of the Rutherford distribution which result in a particle entering the downstream Bambino2 detector");

  addTheta_cmd = new G4UIcmdWithADoubleAndUnit("/Reaction/AddThetaLAB",this);
  addTheta_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  addTheta_cmd->SetGuidance("Add an angle to define desired LAB scattering angle range. This command must always be used two at a time, with the smaller angle coming first. Otherwise it doesn't work.");
  
}
Reaction_Messenger::~Reaction_Messenger() {

  delete reaction_dir;
  
  delete beamZ_cmd;
  delete beamA_cmd;
  delete onlyP_cmd;

  delete recoilZ_cmd;
  delete recoilA_cmd;
  delete onlyR_cmd;

  delete toJanus_cmd;
  delete toUS_cmd;
  delete toDS_cmd;
  delete addTheta_cmd;
  
}

void Reaction_Messenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if(command == beamZ_cmd) {
    reaction->SetBeamZ(beamZ_cmd->GetNewIntValue(newValue));
    G4cout << "Setting projectile nucleus Z to " << newValue << G4endl;
  }

  else if(command == beamA_cmd) {
    reaction->SetBeamA(beamA_cmd->GetNewIntValue(newValue));
    G4cout << "Setting projectile nucleus A to " << newValue << G4endl;
  }

  else if(command == onlyP_cmd) {
    reaction->SetOnlyP();
    G4cout << "Only considering the projectile nucleus when determining the desired scattering angle ranges!" << G4endl;
  }
  
  else if(command == recoilZ_cmd) {
    reaction->SetRecoilZ(recoilZ_cmd->GetNewIntValue(newValue));
    G4cout << "Setting recoil nucleus Z to " << newValue << G4endl;
  }
  
  else if(command == recoilA_cmd) {
    reaction->SetRecoilA(recoilA_cmd->GetNewIntValue(newValue));
    G4cout << "Setting recoil nucleus A to " << newValue << G4endl;
  }

  else if(command == onlyR_cmd) {
    reaction->SetOnlyR();
    G4cout << "Only considering the recoil nucleus when determining the desired scattering angle ranges!" << G4endl;
  }

  else if(command == toJanus_cmd) {
    reaction->Bambino2Thetas();
    G4cout << "Ensuring a particle will always enter Bambino2!" << G4endl;
  }

  else if(command == toUS_cmd) {
    reaction->UpstreamThetas();
    G4cout << "Ensuring the projectile will always enter the upstream Bambino2 detector!" << G4endl;
  }

  else if(command == toDS_cmd) {
    reaction->DownstreamThetas();
    G4cout << "Ensuring a particle will always enter the downstream Bambino2 detector!" << G4endl;
  }

  else if(command == addTheta_cmd) {
    reaction->AddThetaLAB(addTheta_cmd->GetNewDoubleValue(newValue));
    G4cout << "Adding theta = " << newValue << " to list of desired thetas!" << G4endl;
  }
  
  return;
}

