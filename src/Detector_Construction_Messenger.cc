#include "Detector_Construction_Messenger.hh"
#include "Detector_Construction.hh"

Detector_Construction_Messenger::Detector_Construction_Messenger(Detector_Construction* con) : construction(con) {

  //All geometries accessed through this directory
  geometry_dir = new G4UIdirectory("/Geometry/");

  update_cmd = new G4UIcmdWithoutParameter("/Geometry/Update",this);
  update_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  update_cmd->SetGuidance("Mandatory command for updating all geometry");

  //Bambino2 directory
  bambino2_dir = new G4UIdirectory("/Geometry/Bambino2/");

  placeSi_cmd = new G4UIcmdWithoutParameter("/Geometry/Bambino2/Construct",this);
  placeSi_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  placeSi_cmd->SetGuidance("Place the silicon detectors");
  
  offsetUS_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Bambino2/UpstreamOffset",this);
  offsetUS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  offsetUS_cmd->SetGuidance("Set (positive) z-offset of upstream detector (Default: 3 cm)");

  offsetDS_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Bambino2/DownstreamOffset",this);
  offsetDS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  offsetDS_cmd->SetGuidance("Set (positive) z-offset of downstream detector (Default: 3 cm)");

  //SeGA directory
  sega_dir = new G4UIdirectory("/Geometry/SeGA/");

  offsetS_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/SeGA/Offset",this);
  offsetS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  offsetS_cmd->SetGuidance("Set z-offset of SeGA (Default: 0 cm)");

  //Target directory
  target_dir = new G4UIdirectory("/Geometry/Target/");

  placeTarg_cmd = new G4UIcmdWithoutParameter("/Geometry/Target/Construct",this);
  placeTarg_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  placeTarg_cmd->SetGuidance("Place the target");
  
  Z_cmd = new G4UIcmdWithAnInteger("/Geometry/Target/Z",this);
  Z_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  Z_cmd->SetGuidance("Set target Z");

  N_cmd = new G4UIcmdWithAnInteger("/Geometry/Target/N",this);
  N_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  N_cmd->SetGuidance("Set target N");
  
  density_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Target/Density",this);
  density_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  density_cmd->SetGuidance("Set target density");
  
  mass_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Target/Mass",this);
  mass_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  mass_cmd->SetGuidance("Set target mass");
  
  thickness_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Target/Thickness",this);
  thickness_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  thickness_cmd->SetGuidance("Set target linear thickness (length)");
  
  radius_cmd = new G4UIcmdWithADoubleAndUnit("/Geometry/Target/Radius",this);
  radius_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  radius_cmd->SetGuidance("Set target radius");
  
  target_cmd = new G4UIcmdWithAString("/Geometry/Target/StandardTarget",this);
  target_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  target_cmd->SetCandidates("48Ti Ti48 208Pb Pb208 196Pt Pt196");
  target_cmd->SetGuidance("Construct a standard target: 208Pb, 48Ti, or 196Pt");

  print_targ_cmd = new G4UIcmdWithoutParameter("/Geometry/Target/Print",this);
  print_targ_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  print_targ_cmd->SetGuidance("Print target parameters");
  
}

Detector_Construction_Messenger::~Detector_Construction_Messenger() {

  delete geometry_dir;
  delete update_cmd;
  
  delete bambino2_dir;
  delete placeSi_cmd;
  delete offsetUS_cmd;
  delete offsetDS_cmd;

  delete sega_dir;
  delete offsetS_cmd;

  delete target_dir;
  delete placeTarg_cmd;
  delete Z_cmd;
  delete N_cmd;
  delete density_cmd;
  delete mass_cmd;
  delete thickness_cmd;
  delete radius_cmd;
  delete target_cmd;
  
}

void Detector_Construction_Messenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  /////Bambino2 commands/////
  if(command == placeSi_cmd) {
    construction->SetPlaceSilicon();
    G4cout << "Simulation will include the silicon detectors" << G4endl;
  }
  
  else if(command == offsetUS_cmd) {
    construction->SetUS_Offset(offsetUS_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting upstream Bambino2 offset to " << newValue << G4endl;
  }
  
  else if(command == offsetDS_cmd) {
    construction->SetDS_Offset(offsetDS_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting downstream Bambino2 offset to " << newValue << G4endl;
  }
  ///////////////////////////

  /////SeGA commands/////
  else if(command == offsetS_cmd) {
    construction->SetSeGA_Offset(offsetS_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting SeGA offset to " << newValue << G4endl;
  }
  ///////////////////////

  /////Target commands/////
  else if(command == placeTarg_cmd) {
    construction->SetPlaceTarget();
    G4cout << "Simulation will include the target" << G4endl;
  }
  
  else if(command == Z_cmd) {
    construction->SetTargetZ(Z_cmd->GetNewIntValue(newValue));
    G4cout << "Setting target Z to " << newValue << G4endl;
  }

  else if(command == N_cmd) {
    construction->SetTargetN(N_cmd->GetNewIntValue(newValue));
    G4cout << "Setting target N to " << newValue << G4endl;
  }

  else if(command == density_cmd) {
    construction->SetTargetDensity(density_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting target density to " << newValue << G4endl;
  }

  else if(command == mass_cmd) {
    construction->SetTargetMass(mass_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting target mass to " << newValue << G4endl;
  }

  else if(command == thickness_cmd) {
    construction->SetTargetThickness(thickness_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting target thickness to " << newValue << G4endl;
  }

  else if(command == radius_cmd) {
    construction->SetTargetRadius(radius_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting target radius to " << newValue << G4endl;
  }

  else if(command == target_cmd) {
    G4cout << "Setting parameters for a standard " << newValue << " target" << G4endl;
    construction->SetTarget(newValue);
  }

  else if(command == print_targ_cmd) {
    construction->PrintTarget();
  }
  /////////////////////////


  /////Update command/////
  else if(command == update_cmd) {
    construction->Update();
  }
  ////////////////////////

  return;
}
