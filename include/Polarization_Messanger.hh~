#ifndef Excitation_Messenger_h
#define Excitation_Messenger_h 1

#include "Excitation.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
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
  
  G4UIdirectory* rec_dir; //Recoil directory
  
  G4UIcmdWithAString* rLS_cmd; //Level scheme
  G4UIcmdWithAString* rPF_cmd; //Probabilities
  G4UIcmdWithAString* rTF_cmd; //Statistical tensor
  
  G4UIcmdWithAnInteger* rSel_cmd; //Selected state to populate every event
  G4UIcmdWithAnInteger* rCon_cmd; //Only consider this excited state
  G4UIcmdWithADouble* rGSS_cmd; //Ground state spin

  //Control deorientation effect parameters for projectile
  G4UIdirectory* deoP_dir; //Directory
  G4UIcmdWithABool* pCGk_cmd; //Calculate Gk coefficients
  G4UIcmdWithADouble* avjP_cmd; //Average J
  G4UIcmdWithADouble* gamP_cmd; //Gamma
  G4UIcmdWithADouble* lamP_cmd; //Lambda star
  G4UIcmdWithADouble* tauP_cmd; //Tau_C
  G4UIcmdWithADouble* gfcP_cmd; //g-factor
  G4UIcmdWithADouble* fldP_cmd; //Field coefficient
  G4UIcmdWithADouble* expP_cmd; //Field exponent

  //Control deorientation effect parameters for recoil
  G4UIdirectory* deoR_dir; //Directory
  G4UIcmdWithABool* rCGk_cmd; //Calculate Gk coefficients
  G4UIcmdWithADouble* avjR_cmd; //Average J
  G4UIcmdWithADouble* gamR_cmd; //Gamma
  G4UIcmdWithADouble* lamR_cmd; //Lambda star
  G4UIcmdWithADouble* tauR_cmd; //Tau_C
  G4UIcmdWithADouble* gfcR_cmd; //g-factor
  G4UIcmdWithADouble* fldR_cmd; //Field coefficient
  G4UIcmdWithADouble* expR_cmd; //Field exponent
  
};

#endif
