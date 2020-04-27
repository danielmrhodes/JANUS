#ifndef Event_Action_h
#define Event_Action_h 1

#include "G4UserEventAction.hh"
#include "stdio.h"
#include "Tracking_Action.hh"

class Event_Action : public G4UserEventAction {
  
public:
   
  Event_Action();
  ~Event_Action();

  void BeginOfEventAction(const G4Event* evt);
  void EndOfEventAction(const G4Event* evt);

  void SetPerEvent(const int nEvents);
  void SetOutputFile(FILE* f) {output = f;}

  void SetTrackingAction(Tracking_Action* tr) {trkAct = tr;}
  
private:

  Tracking_Action* trkAct;
  
  int perEvt;
  FILE* output;
  
};

#endif
