#ifndef Bambino2_h
#define Bambino2_h 1

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "IonSD.hh"

class Bambino2 {

public:

  Bambino2(G4bool make_sensitive);
  ~Bambino2();

  void Placement(G4LogicalVolume* world, G4double USoff, G4double DSoff);
  

private:

  int nRings;
  int nSectors;
  
  G4double innerRadius;
  G4double outerRadius;
  G4double thickness;
  //G4double DL[4];

  IonSD* TrackerIon;
  
};

#endif
