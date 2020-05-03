#include "Primary_Generator.hh"
#include "Run_Action.hh"
#include "Event_Action.hh"
#include "Tracking_Action.hh"
#include "IonSD.hh"
#include "GammaSD.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"

Run_Action::Run_Action() : fname("output.dat"), output(NULL) {

  messenger = new Run_Action_Messenger(this);
    
}

Run_Action::~Run_Action() {

  delete messenger;
  
}

void Run_Action::BeginOfRunAction(const G4Run* run) {

  G4cout << "\nStarting run!" << G4endl; 
  
  output = fopen(fname.c_str(),"wb");
  int num = run->GetNumberOfEventToBeProcessed();

  G4RunManager* Rman = G4RunManager::GetRunManager();

  Primary_Generator* gen = (Primary_Generator*)Rman->GetUserPrimaryGeneratorAction();
  gen->Update();

  Tracking_Action* trkAct = (Tracking_Action*)Rman->GetUserTrackingAction();
  trkAct->SetMode(gen->GetMode());

  Event_Action* evtAct = (Event_Action*)Rman->GetUserEventAction();
  evtAct->SetOutputFile(output);
  evtAct->SetPerEvent(num);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  switch(gen->GetMode()) {
    case Primary_Generator::MODE::Scattering: {

      IonSD* iSD = (IonSD*)SDman->FindSensitiveDetector("IonTracker");
      iSD->SetProjectileName(gen->GetProjectileName());
      iSD->SetRecoilName(gen->GetRecoilName());

      G4cout << "Simulating " << num << " two-body scattering events" << G4endl;
      break;
  
    }
    case Primary_Generator::MODE::Source: {

      GammaSD* gSD = (GammaSD*)SDman->FindSensitiveDetector("GammaTracker");
      gSD->SetTrackingAction(trkAct);
      
      G4cout << "Simulating " << num << " source gamma-ray events" << G4endl;
      break;

    }
    case Primary_Generator::MODE::Full: {

      IonSD* iSD = (IonSD*)SDman->FindSensitiveDetector("IonTracker");
      iSD->SetProjectileName(gen->GetProjectileName());
      iSD->SetRecoilName(gen->GetRecoilName());
      
      GammaSD* gSD = (GammaSD*)SDman->FindSensitiveDetector("GammaTracker");
      gSD->SetTrackingAction(trkAct);

      trkAct->SetProjectileName(gen->GetProjectileName());
      
      G4cout << "Simulating " << num << " full excitation events" << G4endl;
      break;

    }
  }

  return;
}

void Run_Action::EndOfRunAction(const G4Run*) {

  fclose(output);

  G4cout << "Run Complete!" << G4endl;

  return;
}
