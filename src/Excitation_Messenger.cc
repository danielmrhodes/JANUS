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

  //Gk coefficients
  pCGk_cmd = new G4UIcmdWithABool("/Excitation/Projectile/CalculateGk",this);
  pCGk_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  pCGk_cmd->SetGuidance("Calculate Gk coefficiencts and apply them to the projectile statistical tensor");
  
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

  //Gk coefficients
  rCGk_cmd = new G4UIcmdWithABool("/Excitation/Recoil/CalculateGk",this);
  rCGk_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  rCGk_cmd->SetGuidance("Calculate Gk coefficiencts and apply them to the recoil statistical tensor");

  //Deorentation effect paramter directory
  deo_dir = new G4UIdirectory("/Excitation/DeorientationEffect/");

  //Average J
  avj_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/AverageJ",this);
  avj_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  avj_cmd->SetGuidance("Set the average atomic spin");

  //Gamma
  gam_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/Gamma",this);
  gam_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gam_cmd->SetGuidance("Set the FWHM of the frequency distribution (ps^-1 )");

  //Lambda
  lam_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/Lambda",this);
  lam_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  lam_cmd->SetGuidance("Set the transition rate (ps^-1 ) between static and fluctuating states");

  //TauC
  tau_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/TauC",this);
  tau_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  tau_cmd->SetGuidance("Set the correlation time (ps)");

  //g-factor scaling
  gfm_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/GfactorScaling",this);
  gfm_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  gfm_cmd->SetGuidance("Set the scaling factor applied to the nuclear gyromagnetic ratio g = Z/A");

  //Field coefficient
  fld_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/FieldCoefficient",this);
  fld_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  fld_cmd->SetGuidance("Set the hyperfine field coefficient (10^8 T)");

  //Field exponent
  exp_cmd = new G4UIcmdWithADouble("/Excitation/DeorientationEffect/FieldExponent",this);
  exp_cmd->AvailableForStates(G4ApplicationState::G4State_Idle);
  exp_cmd->SetGuidance("Set the hyperfine field exponent");
  
}

Excitation_Messenger::~Excitation_Messenger() {

  delete excitation_dir;
  delete proj_dir;
  delete rec_dir;
  
  delete pLS_cmd;
  delete pPF_cmd;
  delete pTF_cmd;
  
  delete pSim_cmd;
  delete pSEn_cmd;
  delete pSLt_cmd;
  delete pSel_cmd;
  delete pCon_cmd;
  delete pGSS_cmd;
  delete pCGk_cmd;
  
  delete rLS_cmd;
  delete rPF_cmd;
  delete rTF_cmd;

  delete rSim_cmd;
  delete rSEn_cmd;
  delete rSLt_cmd;
  delete rSel_cmd;
  delete rCon_cmd;
  delete rGSS_cmd;
  delete rCGk_cmd;

  delete deo_dir;
  delete avj_cmd;
  delete gam_cmd;
  delete lam_cmd;
  delete tau_cmd;
  delete gfm_cmd;
  delete fld_cmd;
  delete exp_cmd;
  
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

  ////////////Deorientation effect commands////////////
  else if(command == avj_cmd) {
    excitation->SetAverageJ(avj_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting average atomic spin to " << newValue << G4endl;
  }

  else if(command == gam_cmd) {
    excitation->SetGamma(gam_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting FWHM of frequency distribution to " << newValue << " ps^-1" << G4endl;
  }

  else if(command == lam_cmd) {
    excitation->SetLambdaStar(lam_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting state fluctuation rate to " << newValue << " ps^-1" << G4endl;
  }

  else if(command == tau_cmd) {
    excitation->SetTauC(tau_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting correlation time to " << newValue << " ps" << G4endl;
  }

  else if(command == gfm_cmd) {
    excitation->SetGFacMult(gfm_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting g-factor scaling to " << newValue << G4endl;
  }

  else if(command == fld_cmd) {
    excitation->SetFieldCoef(fld_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field coefficient to " << newValue << "*10^8 T" << G4endl;
  }

  else if(command == exp_cmd) {
    excitation->SetFieldExp(exp_cmd->GetNewDoubleValue(newValue));
    G4cout << "Setting hyperfine field exponent to " << newValue << G4endl;
  }
  /////////////////////////////////////////////////////
  
  return;
  
}
