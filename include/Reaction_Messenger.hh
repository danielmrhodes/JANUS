#ifndef Reaction_Messenger_h
#define Reaction_Messenger_h 1

#include "Reaction.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

class Reaction;
class Reaction_Messenger : public G4UImessenger {

public:

  Reaction_Messenger(Reaction* reac);
  ~Reaction_Messenger();
  
private:

  Reaction* reaction;

  void SetNewValue(G4UIcommand* command, G4String newValue);

  //Reaction directory
  G4UIdirectory* reaction_dir;

  //Projectile commands
  G4UIcmdWithAnInteger* beamZ_cmd;
  G4UIcmdWithAnInteger* beamA_cmd;
  G4UIcmdWithoutParameter* onlyP_cmd;

  //Recoil commands
  G4UIcmdWithAnInteger* recoilZ_cmd;
  G4UIcmdWithAnInteger* recoilA_cmd;
  G4UIcmdWithoutParameter* onlyR_cmd;
  G4UIcmdWithADoubleAndUnit* recoilThresh_cmd;

  //Scattering angle commands
  G4UIcmdWithoutParameter* toJanus_cmd;
  G4UIcmdWithoutParameter* toUS_cmd;
  G4UIcmdWithoutParameter* toDS_cmd;
  G4UIcmdWithoutParameter* rDSpUS_cmd;
  
  G4UIcmdWithADoubleAndUnit* addTheta_cmd;
  
};

#endif
