#ifndef Run_Action_h
#define Run_Action_h 1

#include "G4UserRunAction.hh"
#include "G4String.hh"

class Run_Action : public G4UserRunAction {
  
public:

  Run_Action();
  ~Run_Action();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

};
  
#endif
