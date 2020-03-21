#include "Bambino2.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "IonSD.hh"
#include "G4SDManager.hh"

Bambino2::Bambino2() {

  innerRadius=1.1*cm;
  outerRadius=3.5*cm;
  thickness=300.0*um;

  nRings=24;
  nSectors=32;
  
}

void Bambino2::Placement(G4LogicalVolume* world, G4double USoff, G4double DSoff) {

  //Visualization
  G4VisAttributes* vis1 = new G4VisAttributes(G4Colour::Yellow());
  G4VisAttributes* vis2 = new G4VisAttributes(G4Colour::Red());
  vis1->SetVisibility(true);
  vis1->SetForceSolid(true);
  vis2->SetVisibility(true);
  vis2->SetForceSolid(true);

  //Sensitive Detector
  IonSD* TrackerIon = new IonSD("IonTracker");
  G4SDManager::GetSDMpointer()->AddNewDetector(TrackerIon);

  //Bambino2 material
  G4Material* mat = new G4Material("Si",14,28.0855*g/mole,2.329*g/cm3);

  //Radial and angular spacing
  G4double dr = (outerRadius - innerRadius)/(double)nRings;
  G4double dphi = 2*pi/(double)nSectors;
  
  for(int det=0;det<2;det++) {
    for(int ring=0;ring<nRings;ring++) {
      for(int sec=0;sec<nSectors;sec++) {
      
        G4double in = innerRadius + ring*dr;

        G4Tubs* segs = new G4Tubs("SegS",in,in+dr,thickness/2.0,(sec-0.5)*dphi,dphi);
        G4LogicalVolume* segl = new G4LogicalVolume(segs,mat,"SegL",0,TrackerIon);
	
	if(ring%2) {
	  if(sec%2) {
	    segl->SetVisAttributes(vis1);
	  }
	  else {
	    segl->SetVisAttributes(vis2);
	  }
	}
	else {
	  if(sec%2) {
	    segl->SetVisAttributes(vis2);
	  }
	  else {
	    segl->SetVisAttributes(vis1);
	  }
	}

	G4int copy = det*10000 + (ring+1)*100;
	if(sec<9) {
	    copy += 9-sec; 
	}
	else {
	  copy += 32-(sec-9);
	}
	
	if(det==0) { //Upstream 
	  new G4PVPlacement(0,G4ThreeVector(0,0,-USoff),segl,"Bambino2",world,false,copy,false);
	}
	else { //Downstream
	  new G4PVPlacement(0,G4ThreeVector(0,0,DSoff),segl,"Bambino2",world,false,copy,false);
	}
      
      }
    }
  }
  
  return;
}
