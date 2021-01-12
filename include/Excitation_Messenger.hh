#ifndef Excitation_Messenger_h
#define Excitation_Messenger_h 1

#include "Excitation.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

class Excitation;
class Excitation_Messenger : public G4UImessenger {

public:

  Excitation_Messenger(Excitation* exc);
  ~Excitation_Messenger();
  
private:

  Excitation* excitation;

  void SetNewValue(G4UIcommand* command, G4String newValue);

  G4UIdirectory* excitation_dir; //Excitation directory

  G4UIdirectory* proj_dir; //Projectile directory
  
  G4UIcmdWithAString* pLS_cmd; //Level scheme
  G4UIcmdWithAString* pPF_cmd; //Probabilities
  G4UIcmdWithAString* pTF_cmd; //Statistical tensor

  G4UIcmdWithAnInteger* pSel_cmd; //Selected state to populate every event
  G4UIcmdWithAnInteger* pCon_cmd; //Only consider this excited state
  G4UIcmdWithADouble* pGSS_cmd; //Ground state spin
  G4UIcmdWithABool* pCGk_cmd; //Calculate Gk coefficients
  
  G4UIdirectory* rec_dir; //Recoil directory
  
  G4UIcmdWithAString* rLS_cmd; //Level scheme
  G4UIcmdWithAString* rPF_cmd; //Probabilities
  G4UIcmdWithAString* rTF_cmd; //Statistical tensor
  
  G4UIcmdWithAnInteger* rSel_cmd; //Selected state to populate every event
  G4UIcmdWithAnInteger* rCon_cmd; //Only consider this excited state
  G4UIcmdWithADouble* rGSS_cmd; //Ground state spin
  G4UIcmdWithABool* rCGk_cmd; //Calculate Gk coefficients

  G4UIdirectory* deo_dir; //Deorientation effect directory

  //Control deorientation effect paramerts
  G4UIcmdWithADouble* avj_cmd; //Average J
  G4UIcmdWithADouble* gam_cmd; //Gamma
  G4UIcmdWithADouble* lam_cmd; //Lambda star
  G4UIcmdWithADouble* tau_cmd; //Tau_C
  G4UIcmdWithADouble* gfm_cmd; //g-factor multiplier
  G4UIcmdWithADouble* fld_cmd; //Field coefficient
  G4UIcmdWithADouble* exp_cmd; //Field exponent
  
};

#endif
