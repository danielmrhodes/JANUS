#ifndef Run_Action_h
#define Run_Action_h 1

#include "G4UserRunAction.hh"
#include "Run_Action_Messenger.hh"
#include "G4String.hh"

class Run_Action_Messenger;

class Run_Action : public G4UserRunAction {
  
public:

  Run_Action();
  ~Run_Action();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void SetFileName(G4String name) {fname = name;}

private:

  Run_Action_Messenger* messenger;
  
  G4String fname;
  FILE* output;

};
  
#endif
