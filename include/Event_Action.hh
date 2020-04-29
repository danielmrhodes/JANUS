#ifndef Event_Action_h
#define Event_Action_h 1

#include "G4UserEventAction.hh"
#include "stdio.h"

class Event_Action : public G4UserEventAction {
  
public:
   
  Event_Action();
  ~Event_Action();

  void BeginOfEventAction(const G4Event* evt);
  void EndOfEventAction(const G4Event* evt);

  void SetPerEvent(const int nEvents);
  void SetOutputFile(FILE* f) {output = f;}
  
private:
  
  int perEvt;
  FILE* output;
  
};

#endif
