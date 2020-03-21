#include "GammaSD.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4VProcess.hh"

GammaSD::GammaSD(G4String name) : G4VSensitiveDetector(name) {
  collectionName.insert("gammaCollection");
}

GammaSD::~GammaSD() {}

void GammaSD::Initialize(G4HCofThisEvent*) {

  HC = new Gamma_Hit_Collection(SensitiveDetectorName,collectionName[0]);
  
  return;
}

G4bool GammaSD::ProcessHits(G4Step* step, G4TouchableHistory*) {

  if(step->GetTotalEnergyDeposit() > 0.0) {

    Gamma_Hit* hit = new Gamma_Hit();
    hit->SetPos(step->GetPreStepPoint()->GetPosition());
    hit->SetDetSeg(step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
    hit->SetEdep(step->GetTotalEnergyDeposit());
    
    HC->insert(hit);
    
  }
 
  return true;
}

void GammaSD::EndOfEvent(G4HCofThisEvent* HCE) {

 label:

  for(int i=0;i<HC->entries();i++) {
    for(int j=i+1;j<HC->entries();j++) {

      Gamma_Hit* hit1 = (Gamma_Hit*)HC->GetHit(i);
      Gamma_Hit* hit2 = (Gamma_Hit*)HC->GetHit(j);

      if(hit1->GetSegment() == hit2->GetSegment() && hit1->GetDetector() == hit2->GetDetector()) {
	
	hit1->SetEdep(hit1->GetEdep()+hit2->GetEdep());
	
	std::vector<Gamma_Hit*>* vec = HC->GetVector();
	vec->erase(vec->begin()+j);

	goto label;
	
      }
    }
  }

  G4double cores[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(int i=0;i<HC->entries();i++) {
    
    Gamma_Hit* hit = (Gamma_Hit*)HC->GetHit(i);	
    cores[hit->GetDetector()-1] += hit->GetEdep();
    
  }

  for(int i=0;i<16;i++) {
    if(cores[i] > 0) {

      Gamma_Hit* hit = new Gamma_Hit();
      hit->SetPos(G4ThreeVector());
      hit->SetDetSeg(100*(i+1));
      hit->SetEdep(cores[i]);

      HC->insert(hit);
      
    }
  }
  
  HCE->AddHitsCollection(HCE->GetNumberOfCollections(),HC);

  return;
  
}
