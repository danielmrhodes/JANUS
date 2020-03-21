#ifndef Excitation_Messenger_h
#define Excitation_Messenger_h 1

#include "Excitation.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

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

  G4UIcmdWithoutParameter* pSim_cmd; //Simple level
  G4UIcmdWithADoubleAndUnit* pSEn_cmd; //Simple level energy
  G4UIcmdWithADoubleAndUnit* pSLt_cmd; //Simple level lifetime
  G4UIcmdWithAnInteger* pSel_cmd; //Selected state to populate
  
  G4UIdirectory* rec_dir; //Recoil directory
  
  G4UIcmdWithAString* rLS_cmd; //Level scheme
  G4UIcmdWithAString* rPF_cmd; //Probabilities

  G4UIcmdWithoutParameter* rSim_cmd; //Simple level
  G4UIcmdWithADoubleAndUnit* rSEn_cmd; //Simple level energy
  G4UIcmdWithADoubleAndUnit* rSLt_cmd; //Simple level lifetime
  G4UIcmdWithAnInteger* rSel_cmd; //Selected state to populate
  
};

#endif
