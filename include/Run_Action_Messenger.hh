#ifndef Run_Action_Messenger_h
#define Run_Action_Messenger_h 1

#include "Run_Action.hh"
#include "G4UImessenger.hh"
//#include "G4UIcmdWithoutParameter.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"
//#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

class Run_Action;

class Run_Action_Messenger : public G4UImessenger {

public:
  
  Run_Action_Messenger(Run_Action* rAct);
  ~Run_Action_Messenger();
  
private:

  Run_Action* runAction;

  void SetNewValue(G4UIcommand* command, G4String newValue);

  //Output directory
  G4UIdirectory* output_dir;
  
  G4UIcmdWithAString* name_cmd;
  
};

#endif
