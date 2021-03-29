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

  //Statistical tensor file
  pTF_cmd = new G4UIcmdWithAString("/Excitation/Projectile/StatisticalTensors",this);
  pTF_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pTF_cmd->SetGuidance("Set name of projectile statistical tensor file");
  
  //Selected state
  pSel_cmd = new G4UIcmdWithAnInteger("/Excitation/Projectile/PopulateState",this);
  pSel_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pSel_cmd->SetGuidance("Choose one state to always populate in the projectile");

  //Considered state
  pCon_cmd = new G4UIcmdWithAnInteger("/Excitation/Projectile/OnlyConsiderState",this);
  pCon_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pCon_cmd->SetGuidance("Only allow one excited state to be populated in the projectile, but use real probabilities");

  //Ground state spin
  pGSS_cmd = new G4UIcmdWithADouble("/Excitation/Projectile/GroundStateSpin",this);
  pGSS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pGSS_cmd->SetGuidance("Set spin of projectile ground state");
  
  //Recoil directory
  rec_dir = new G4UIdirectory("/Excitation/Recoil/");

  //Level scheme file
  rLS_cmd = new G4UIcmdWithAString("/Excitation/Recoil/LevelScheme",this);
  rLS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rLS_cmd->SetGuidance("Set name of recoil level scheme file");

  //Probabilities file
  rPF_cmd = new G4UIcmdWithAString("/Excitation/Recoil/Probabilities",this);
  rPF_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rPF_cmd->SetGuidance("Set name of recoil excitation probabilities file");

  //Statistical tensor file
  rTF_cmd = new G4UIcmdWithAString("/Excitation/Recoil/StatisticalTensors",this);
  rTF_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rTF_cmd->SetGuidance("Set name of recoil statistical tensor file");

  //Selected state
  rSel_cmd = new G4UIcmdWithAnInteger("/Excitation/Recoil/PopulateState",this);
  rSel_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rSel_cmd->SetGuidance("Choose one state to always populate in the recoil");

  //Considered state
  rCon_cmd = new G4UIcmdWithAnInteger("/Excitation/Recoil/OnlyConsiderState",this);
  rCon_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rCon_cmd->SetGuidance("Only allow one excited state to be populated in the recoil, but use real probabilities");

  //Ground state spin
  rGSS_cmd = new G4UIcmdWithADouble("/Excitation/Recoil/GroundStateSpin",this);
  rGSS_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rGSS_cmd->SetGuidance("Set spin of recoil ground state");

  //Deorentation effect paramter directory for projectile
  deoP_dir = new G4UIdirectory("/DeorientationEffect/Projectile/");

  //Gk coefficients
  pCGk_cmd = new G4UIcmdWithABool("/DeorientationEffect/Projectile/CalculateGk",this);
  pCGk_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pCGk_cmd->SetGuidance("Calculate Gk coefficiencts and apply them to the projectile statistical tensors");

  //Average J
  avjP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/AverageJ",this);
  avjP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  avjP_cmd->SetGuidance("Set the average atomic spin of the projectile");

  //Gamma
  gamP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/Gamma",this);
  gamP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gamP_cmd->SetGuidance("Set the FWHM of the frequency distribution (ps^-1 ) in the projectile");

  //Lambda
  lamP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/Lambda",this);
  lamP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  lamP_cmd->SetGuidance("Set the transition rate (ps^-1 ) between static and fluctuating states for the projectile");

  //TauC
  tauP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/TauC",this);
  tauP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  tauP_cmd->SetGuidance("Set the correlation time (ps) of the projectile");

  //g-factor
  gfcP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/GFactor",this);
  gfcP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gfcP_cmd->SetGuidance("Set the nuclear gyromagnetic ratio (g-factor) for the projectile");

  //Field coefficient
  fldP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/FieldCoefficient",this);
  fldP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  fldP_cmd->SetGuidance("Set the hyperfine field coefficient (10^8 T) for the projectile");

  //Field exponent
  expP_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Projectile/FieldExponent",this);
  expP_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  expP_cmd->SetGuidance("Set the hyperfine field exponent for the projectile");

  //Deorentation effect paramter directory for recoil
  deoR_dir = new G4UIdirectory("/DeorientationEffect/Recoil/");

  //Gk coefficients
  rCGk_cmd = new G4UIcmdWithABool("/DeorientationEffect/Recoil/CalculateGk",this);
  rCGk_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rCGk_cmd->SetGuidance("Calculate Gk coefficiencts and apply them to the recoil statistical tensors");

  //Average J
  avjR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/AverageJ",this);
  avjR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  avjR_cmd->SetGuidance("Set the average atomic spin of the recoil");

  //Gamma
  gamR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/Gamma",this);
  gamR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gamR_cmd->SetGuidance("Set the FWHM of the frequency distribution (ps^-1 ) in the recoil");

  //Lambda
  lamR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/Lambda",this);
  lamR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  lamR_cmd->SetGuidance("Set the transition rate (ps^-1 ) between static and fluctuating states for the recoil");

  //TauC
  tauR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/TauC",this);
  tauR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  tauR_cmd->SetGuidance("Set the correlation time (ps) of the recoil");

  //g-factor scaling
  gfcR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/GFactor",this);
  gfcR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gfcR_cmd->SetGuidance("Set the nuclear gyromagnetic ratio (g-factor) for the recoil");

  //Field coefficient
  fldR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/FieldCoefficient",this);
  fldR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  fldR_cmd->SetGuidance("Set the hyperfine field coefficient (10^8 T) for the recoil");

  //Field exponent
  expR_cmd = new G4UIcmdWithADouble("/DeorientationEffect/Recoil/FieldExponent",this);
  expR_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  expR_cmd->SetGuidance("Set the hyperfine field exponent for the recoil");
  
}

Excitation_Messenger::~Excitation_Messenger() {

  delete excitation_dir;
  delete proj_dir;
  delete rec_dir;
  
  delete pLS_cmd;
  delete pPF_cmd;
  delete pTF_cmd;
  
  delete pSel_cmd;
  delete pCon_cmd;
  delete pGSS_cmd;
  delete pCGk_cmd;
  
  delete rLS_cmd;
  delete rPF_cmd;
  delete rTF_cmd;
  
  delete rSel_cmd;
  delete rCon_cmd;
  delete rGSS_cmd;
  delete rCGk_cmd;

  delete deoP_dir;
  delete avjP_cmd;
  delete gamP_cmd;
  delete lamP_cmd;
  delete tauP_cmd;
  delete gfcP_cmd;
  delete fldP_cmd;
  delete expP_cmd;

  delete deoR_dir;
  delete avjR_cmd;
  delete gamR_cmd;
  delete lamR_cmd;
  delete tauR_cmd;
  delete gfcR_cmd;
  delete fldR_cmd;
  delete expR_cmd;
  
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

  else if(command == pTF_cmd) {
    excitation->SetProjTensorFile(newValue);
    G4cout << "Setting projectile statistical tensor file to " << newValue << G4endl;
  }

  else if(command == pSel_cmd) {
    excitation->FixProjState(pSel_cmd->GetNewIntValue(newValue));
    G4cout << "Selecting state " << newValue << " to populate in the projectile" << G4endl;
  }

  else if(command == pCon_cmd) {
    excitation->OnlyConsiderProjState(pCon_cmd->GetNewIntValue(newValue));
    G4cout << "Will only consider state " << newValue << " when populating excited states in the projectile"
	   << G4endl;
  }

  else if(command == pGSS_cmd) {
    excitation->SetProjGSS(pGSS_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting ground state spin of projectile to " << newValue << G4endl;
  }

  else if(command == pCGk_cmd) {
    excitation->SetProjCalcGk(pCGk_cmd->GetNewBoolValue(newValue));
    G4cout << "Setting flag for projectile Gk calculation to " << newValue << G4endl;
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

  else if(command == rTF_cmd) {
    excitation->SetRecTensorFile(newValue);
    G4cout << "Setting recoil statistical tensor file to " << newValue << G4endl;
  }

  else if(command == rSel_cmd) {
    excitation->FixRecState(rSel_cmd->GetNewIntValue(newValue));
    G4cout << "Selecting state " << newValue << " to populate in the recoil" << G4endl;
  }

  else if(command == rCon_cmd) {
    excitation->OnlyConsiderRecState(rCon_cmd->GetNewIntValue(newValue));
    G4cout << "Will only consider state " << newValue << " when populating excited states in the recoil"
	   << G4endl;
  }

  else if(command == rGSS_cmd) {
    excitation->SetRecGSS(rGSS_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting ground state spin of recoil to " << newValue << G4endl;
  }

  else if(command == rCGk_cmd) {
    excitation->SetRecCalcGk(rCGk_cmd->GetNewBoolValue(newValue));
    G4cout << "Setting flag for recoil Gk calculation to " << newValue << G4endl;
  }
  ///////////////////////////////////////

  ////////////Projectile deorientation effect commands////////////
  else if(command == avjP_cmd) {
    excitation->SetProjAverageJ(avjP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting average atomic spin of the projectle to " << newValue << G4endl;
  }

  else if(command == gamP_cmd) {
    excitation->SetProjGamma(gamP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting the FWHM of the frequency distribution in the projectle to " << newValue << " ps^-1"
	   << G4endl;
  }

  else if(command == lamP_cmd) {
    excitation->SetProjLambdaStar(lamP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting state fluctuation rate in the projectle to " << newValue << " ps^-1" << G4endl;
  }

  else if(command == tauP_cmd) {
    excitation->SetProjTauC(tauP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting the correlation time in the projectle to " << newValue << " ps" << G4endl;
  }

  else if(command == gfcP_cmd) {
    excitation->SetProjGFac(gfcP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting g-factor for the projectle to " << newValue << G4endl;
  }

  else if(command == fldP_cmd) {
    excitation->SetProjFieldCoef(fldP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field coefficient in the projectle to " << newValue << "*10^8 T" << G4endl;
  }

  else if(command == expP_cmd) {
    excitation->SetProjFieldExp(expP_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field exponent in the projectle to " << newValue << G4endl;
  }
  ////////////////////////////////////////////////////////////////

  ////////////Recoil deorientation effect commands////////////
  else if(command == avjR_cmd) {
    excitation->SetProjAverageJ(avjR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting average atomic spin of the recoil to " << newValue << G4endl;
  }

  else if(command == gamR_cmd) {
    excitation->SetProjGamma(gamR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting FWHM of frequency distribution in the recoil to " << newValue << " ps^-1" << G4endl;
  }

  else if(command == lamR_cmd) {
    excitation->SetProjLambdaStar(lamR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting state fluctuation rate in the recoil to " << newValue << " ps^-1" << G4endl;
  }

  else if(command == tauR_cmd) {
    excitation->SetProjTauC(tauR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting correlation time in the recoil to " << newValue << " ps" << G4endl;
  }

  else if(command == gfcR_cmd) {
    excitation->SetProjGFac(gfcR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting g-factor for the recoil to " << newValue << G4endl;
  }

  else if(command == fldR_cmd) {
    excitation->SetProjFieldCoef(fldR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field coefficient in the recoil to " << newValue << "*10^8 T" << G4endl;
  }

  else if(command == expR_cmd) {
    excitation->SetProjFieldExp(expR_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field exponent in the recoil to " << newValue << G4endl;
  }
  ////////////////////////////////////////////////////////////
  
  return;
  
}
