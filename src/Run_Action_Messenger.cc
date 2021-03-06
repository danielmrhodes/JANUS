#include "Run_Action_Messenger.hh"

Run_Action_Messenger::Run_Action_Messenger(Run_Action* rAct) : runAction(rAct) {

  //All info about the output accessed through this directory
  output_dir = new G4UIdirectory("/Output/");

  name_cmd = new G4UIcmdWithAString("/Output/FileName",this);
  name_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  name_cmd->SetGuidance("Set output file name.");
  
}

Run_Action_Messenger::~Run_Action_Messenger() {

  delete output_dir;
  delete name_cmd;
  
}

void Run_Action_Messenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if(command == name_cmd) {
    runAction->SetFileName(newValue);
    G4cout << "Setting output file name to " << newValue << G4endl;
  }
  
  return;
}

