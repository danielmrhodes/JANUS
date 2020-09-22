#include "Event_Action.hh"
#include "G4Event.hh"
#include "Ion_Hit.hh"
#include "Gamma_Hit.hh"
#include "G4SystemOfUnits.hh"
#include "Data_Format.hh"

Event_Action::Event_Action() {}
Event_Action::~Event_Action() {}

void Event_Action::BeginOfEventAction(const G4Event* evt) {

  int id = evt->GetEventID();
  if(!(id%perEvt)) {
    G4cout << "Event " << id << G4endl;
  }

  return;
}

void Event_Action::EndOfEventAction(const G4Event* evt) {
  
  JANUSData data;
  int nB = 0;
  int nS = 0;
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  for(int i=0;i<HCE->GetNumberOfCollections();i++) {
    if(HCE->GetHC(i)->GetName() == "ionCollection") {

      Ion_Hit_Collection* iHC = (Ion_Hit_Collection*)HCE->GetHC(i);  
      for(int j=0;j<iHC->entries();j++) {

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
      for(int j=0;j<gHC->entries();j++) {

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
  
  if(nB > 0 || nS > 0) {

    Header header;
    header.evtNum = evt->GetEventID();
    header.nBdata = nB;
    header.nSdata = nS;
  
    fwrite(&header,header.bytes(),1,output);
    fwrite(&data.bData,nB*sizeof(Bambino2Data),1,output);
    fwrite(&data.sData,nS*sizeof(SegaData),1,output);

  }
  
  return;
}

void Event_Action::SetPerEvent(const int nEvents) {

  if(nEvents > 2000000) {
    perEvt = 500000;
  }
  else if(nEvents > 1000000) {
    perEvt = 200000;
  }
  else if(nEvents > 500000) {
    perEvt = 100000;
  }
  else if(nEvents > 200000) {
    perEvt = 50000;
  }
  else if(nEvents > 100000) {
    perEvt = 20000;
  }
  else if(nEvents > 50000) {
    perEvt = 10000;
  }
  else if(nEvents > 20000) {
    perEvt = 5000;
  }
  else if(nEvents > 10000) {
    perEvt = 2000;
  }
  else if(nEvents > 5000) {
    perEvt = 1000;
  }
  else if(nEvents > 2000) {
    perEvt = 500;
  }
  else if(nEvents > 1000) {
    perEvt = 200;
  }
  else if(nEvents > 500) {
    perEvt = 100;
  }
  else if(nEvents > 200) {
    perEvt = 50;
  }
  else if(nEvents > 100) {
    perEvt = 20;
  }
  else if(nEvents > 50) {
    perEvt = 10;
  }
  else if(nEvents > 20) {
    perEvt = 5;
  }
  else if(nEvents > 10) {
    perEvt = 2;
  }
  else {
    perEvt = 1;
  }
  
  return;
}
