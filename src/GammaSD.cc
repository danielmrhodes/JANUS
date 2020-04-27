#include "GammaSD.hh"
#include "G4SystemOfUnits.hh"

GammaSD::GammaSD(G4String name) : G4VSensitiveDetector(name) {
  collectionName.insert("gammaCollection");
}

GammaSD::~GammaSD() {}

void GammaSD::Initialize(G4HCofThisEvent*) {

  HC = new Gamma_Hit_Collection(SensitiveDetectorName,collectionName[0]);
  detMap.clear();
  
  return;
}

G4bool GammaSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
  
  if(step->GetTotalEnergyDeposit() > 0.0) {

    Gamma_Hit* hit = new Gamma_Hit();
    hit->SetPos(step->GetPreStepPoint()->GetPosition());
    hit->SetDetSeg(step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo());
    hit->SetEdep(step->GetTotalEnergyDeposit());
    HC->insert(hit);
      
    int det = hit->GetDetector();
    int id = step->GetTrack()->GetTrackID();
    if(std::find(detMap[det].begin(),detMap[det].end(),id) == detMap[det].end()) {
      detMap[det].push_back(id);
    }
      
    
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

  std::map<G4int,std::vector<G4int>> idMap = trkAct->GetIDMap();
  std::map<G4int,G4double> enMap = trkAct->GetEnergyMap();
  for(int i=0;i<16;i++) {
    if(cores[i] > 0) {

      Gamma_Hit* hit = new Gamma_Hit();
      hit->SetPos(G4ThreeVector());
      hit->SetDetSeg(100*(i+1));
      hit->SetEdep(cores[i]);

      G4bool broken = false;
      for(auto it = idMap.begin();it != idMap.end();++it) {
	
	if(std::find(detMap[i+1].begin(),detMap[i+1].end(),it->first) == detMap[i+1].end()) {
	  continue;
	}

	broken = false;
	for(unsigned int id : detMap[i+1]) {
	  if(std::find(it->second.begin(),it->second.end(),id) == it->second.end()) {
	    broken = true;
	    break;
	  }
	}

	if(broken) {
	  break;
	}

	G4double diff = enMap[it->first] - cores[i];
	if(diff*diff < (0.01*keV)*(0.01*keV)) {
	  hit->SetFEP();
	}
	  
      }
      if(broken) {
	HC->insert(hit);
	continue;
      }
      
      HC->insert(hit);
      
    }
  }
  
  HCE->AddHitsCollection(HCE->GetNumberOfCollections(),HC);

  return;
  
}
