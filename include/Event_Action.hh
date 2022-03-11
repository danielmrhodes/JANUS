#ifndef Event_Action_h
#define Event_Action_h 1

#include "Event_Action_Messenger.hh"
#include "G4UserEventAction.hh"
#include "G4String.hh"
#include "stdio.h"

class Event_Action_Messenger;
class Event_Action : public G4UserEventAction {
  
public:
   
  Event_Action();
  ~Event_Action();

  void BeginOfEventAction(const G4Event* evt);
  void EndOfEventAction(const G4Event* evt);

  G4String GetOutputFileName() {return fname;}
  FILE* GetOutputFile() {return output;}
  
  void SetNEvents(G4int n) {nEvents = n;} 
  void SetOutputFile(FILE* f) {output = f;}
  void SetOutputFileName(G4String n) {fname = n;}
  void OWC() {owc = true;}
  
private:

  Event_Action_Messenger* messenger;
  
  G4int nEvents;
  G4bool owc;

  G4String fname;
  FILE* output;
  
};

#endif
