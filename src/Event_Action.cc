#include "Event_Action.hh"
#include "G4Event.hh"
#include "Ion_Hit.hh"
#include "Gamma_Hit.hh"
#include "G4SystemOfUnits.hh"
#include "Data_Format.hh"

Event_Action::Event_Action() : nEvents(0), owc(false), fname("output.dat"), output(NULL) {

  messenger = new Event_Action_Messenger(this);
  
}
Event_Action::~Event_Action() {

  delete messenger;
  
}

void Event_Action::BeginOfEventAction(const G4Event* evt) {

  G4int id = evt->GetEventID();
  G4cout << "Event " << id+1 << " (" << 100*(id+1)/nEvents << "%)\r" << std::flush;

  return;
}

void Event_Action::EndOfEventAction(const G4Event* evt) {
  
  JANUSData data;
  G4int nB = 0;
  G4int nS = 0;
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  for(G4int i=0;i<HCE->GetNumberOfCollections();i++) {
    if(HCE->GetHC(i)->GetName() == "ionCollection") {

      Ion_Hit_Collection* iHC = (Ion_Hit_Collection*)HCE->GetHC(i);  
      for(G4int j=0;j<iHC->entries();j++) {

	if(nB > 9) {
	  G4cout << "Too many ion hits!" << G4endl;
	  break;
	}
    
        Ion_Hit* hit = (Ion_Hit*)iHC->GetHit(j);
        G4ThreeVector pos = hit->GetPos();

	data.bData[nB]= {hit->GetDetector(),hit->GetRing(),hit->GetSector(),
		         hit->GetEdep()/MeV,pos.x()/cm,pos.y()/cm,pos.z()/cm,
	                 hit->IsProjectile(),hit->IsRecoil()};
	 
	nB++;

      }
    }
    else if(HCE->GetHC(i)->GetName() == "gammaCollection") {

      Gamma_Hit_Collection* gHC = (Gamma_Hit_Collection*)HCE->GetHC(i);
      for(G4int j=0;j<gHC->entries();j++) {

	if(nS > 49) {
	  G4cout << "Too many gamma hits!" << G4endl;
	  break;
	}
    
        Gamma_Hit* hit = (Gamma_Hit*)gHC->GetHit(j);
        G4ThreeVector pos = hit->GetPos();

	data.sData[nS] = {hit->GetDetector(),hit->GetSegment(),hit->GetEdep()/keV,
			  pos.x()/cm,pos.y()/cm,pos.z()/cm,hit->IsFEP(),hit->IsProjFEP()};

	nS++;

      }
    }
  }
  
  if(nB == 0 && nS == 0)
    return;

  if(owc && (nB == 0 || nS == 0))
    return;
    
  Header header;
  header.evtNum = evt->GetEventID();
  header.nBdata = nB;
  header.nSdata = nS;
  
  fwrite(&header,header.bytes(),1,output);
  fwrite(&data.bData,sizeof(Bambino2Data),nB,output);
  fwrite(&data.sData,sizeof(SegaData),nS,output);
  
  return;
}
