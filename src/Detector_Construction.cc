#include "Detector_Construction.hh"
#include "Detector_Construction_Messenger.hh"
#include "Bambino2.hh"
#include "SeGA.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

Detector_Construction::Detector_Construction() {

  messenger = new Detector_Construction_Messenger(this);

  US_Offset = 3.0*cm;
  DS_Offset = 3.0*cm;
  SeGA_Offset = 0.0*cm;

  target_Z = 82;
  target_N = 126;
  target_density = 11.382*g/cm3;
  target_mass = 207.97665*g/mole;
  target_thickness = 882*nm;
  target_radius = 0.5*cm;
  
}

Detector_Construction::~Detector_Construction() {

  delete messenger;
  
}

G4VPhysicalVolume* Detector_Construction::Construct() {

  //NIST database for materials
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic"); //World material
  
  //Make the world
  G4Box* solid_world = new G4Box("World_Solid",1.0*m,1.0*m,1.0*m);
  logic_world = new G4LogicalVolume(solid_world,world_mat,"World_Logical");
  world = new G4PVPlacement(0,G4ThreeVector(),logic_world,"World",0,false,0,false);

  return world;
}

G4VPhysicalVolume* Detector_Construction::PlaceVolumes() {
  
  //Target material (isotopically pure)
  target_mat = new G4Material("target_mat",target_density,1); //Bulk material (1 component)
  G4Element* target_ele = new G4Element("target_ele","target_symbol",1); //Element (1 isoptope)
  G4Isotope* target_iso = new G4Isotope("target_iso",target_Z,target_N,target_mass); //The isotope
  target_ele->AddIsotope(target_iso,1.0);
  target_mat->AddElement(target_ele,1.0);

  //Make the target
  G4Tubs* solid_target = new G4Tubs("Target_Sol",0*cm,target_radius,target_thickness/2.0,0.0*deg,360.0*deg);
  G4LogicalVolume* logic_target = new G4LogicalVolume(solid_target,target_mat,"Target_Logical");
  new G4PVPlacement(0,G4ThreeVector(),logic_target,"Target",logic_world,false,0,false);
  
  //Make Bambino2
  Bambino2* bam = new Bambino2();
  bam->Placement(logic_world,US_Offset,DS_Offset);

  //Make SeGA
  SeGA* seg = new SeGA();
  seg->Placement(logic_world,SeGA_Offset);

  //Beam tube and fram material (aluminium)
  G4Material* Al = new G4Material("Aluminum",13,26.98*g/mole,2.7*g/cm3);

  //Make the beam tube
  G4Tubs* solid_BT = new G4Tubs("BT_Sol",7.366*cm,7.62*cm,44*cm,0*deg,360*deg);
  G4LogicalVolume* logic_BT = new G4LogicalVolume(solid_BT,Al,"BT_Logical");
  new G4PVPlacement(0,G4ThreeVector(),logic_BT,"Beam_Tube",logic_world,false,0,false);

  //Make frame
  /*
  G4Tubs* solid_frame = new G4Tubs("frame_Sol",47.6*cm,65.7*cm,0.635*cm,0*deg,360*deg);
  G4LogicalVolume* logic_frame = new G4LogicalVolume(solid_frame,Al,"frame_Logical");
  new G4PVPlacement(0,G4ThreeVector(0,0,16*cm),logic_frame,"frame1",logic_world,false,0,false);
  new G4PVPlacement(0,G4ThreeVector(0,0,-16*cm),logic_frame,"frame2",logic_world,false,0,false);
  */
  
  G4Box* solid_plate = new G4Box("plate_Sol",0.75*m,0.5*m,0.635*cm);
  G4Tubs* tub = new G4Tubs("tubs",0*cm,33*cm,1.0*cm,0*deg,360*deg);
  G4SubtractionSolid* solid_face = new G4SubtractionSolid("face_Sol",solid_plate,tub);
  G4LogicalVolume* logic_face = new G4LogicalVolume(solid_face,Al,"face_Logical");
  new G4PVPlacement(0,G4ThreeVector(0,0,8*cm+SeGA_Offset),logic_face,"face1",logic_world,false,0,false);
  new G4PVPlacement(0,G4ThreeVector(0,0,-8*cm+SeGA_Offset),logic_face,"face2",logic_world,false,0,false);

  G4RotationMatrix* Rot = new G4RotationMatrix();
  Rot->rotateX(90*deg);
  
  G4LogicalVolume* logic_plate = new G4LogicalVolume(solid_plate,Al,"plate_Logical");
  new G4PVPlacement(Rot,G4ThreeVector(0,-0.65*m,SeGA_Offset),logic_plate,"plate1",logic_world,false,0,false);
  new G4PVPlacement(Rot,G4ThreeVector(0,-0.75*m,SeGA_Offset),logic_plate,"plate2",logic_world,false,0,false);
  new G4PVPlacement(Rot,G4ThreeVector(0,-0.85*m,SeGA_Offset),logic_plate,"plate3",logic_world,false,0,false);
  new G4PVPlacement(Rot,G4ThreeVector(0,-0.95*m,SeGA_Offset),logic_plate,"plate4",logic_world,false,0,false);

  G4Box* solid_arm = new G4Box("arm_Sol",0.635*cm,0.65*m,4*cm);
  G4LogicalVolume* logic_arm = new G4LogicalVolume(solid_arm,Al,"arm_Logical");
  new G4PVPlacement(0,G4ThreeVector(0.75*m,0,SeGA_Offset),logic_arm,"arm1",logic_world,false,0,false);
  new G4PVPlacement(0,G4ThreeVector(-0.75*m,0,SeGA_Offset),logic_arm,"arm2",logic_world,false,0,false);

  G4Box* solid_top = new G4Box("top_Sol",0.75*m,0.635*cm,4*cm);
  G4LogicalVolume* logic_top = new G4LogicalVolume(solid_top,Al,"top_Logical");
  new G4PVPlacement(0,G4ThreeVector(0,0.65*m,SeGA_Offset),logic_top,"top1",logic_world,false,0,false);
  
  //Gate valve material
  G4Element* C  = new G4Element("Carbon","C",6,12.011*g/mole);
  G4Element* Co = new G4Element("Cobalt","Co",27,58.9332*g/mole);
  G4Element* Fe = new G4Element("Iron","Fe",26,55.85*g/mole);
  G4Material* GV_mat = new G4Material("ssteel",7.7*g/cm3,3);
  GV_mat->AddElement(C,0.04);
  GV_mat->AddElement(Fe,0.88);
  GV_mat->AddElement(Co,0.08);

  //Make gate valve
  G4Tubs* solid_GV = new G4Tubs("GV_Sol",9.3*cm,17*cm,5.9*cm,0*deg,360*deg);
  G4LogicalVolume* logic_GV = new G4LogicalVolume(solid_GV,GV_mat,"GV_Logical");
  new G4PVPlacement(0,G4ThreeVector(0,0,50*cm),logic_GV,"GV",logic_world,false,0,false);

  //Visualization
  G4VisAttributes* vis = new G4VisAttributes(G4Colour::Grey());
  vis->SetVisibility(true);
  vis->SetForceSolid(false);
  //vis->SetLineStyle(G4VisAttributes::LineStyle::unbroken);
  
  logic_target->SetVisAttributes(vis);
  logic_BT->SetVisAttributes(vis);
  logic_face->SetVisAttributes(vis);
  logic_GV->SetVisAttributes(vis);
  logic_plate->SetVisAttributes(vis);
  logic_arm->SetVisAttributes(vis);
  
  return world;
}

void Detector_Construction::Update() {

  static bool volumes_placed = false;
  if(!volumes_placed) {
    
    G4cout << "Updating geometry with input parameters." << G4endl;
    G4RunManager::GetRunManager()->DefineWorldVolume(PlaceVolumes());
    volumes_placed = true;
    
  }
  else {
    G4cout << "Geometry has already been updated. Cannot make additional changes." << G4endl;
  }

  return;

}

void Detector_Construction::SetTarget(G4String target) {

  //Need to get exact parameters for every target...
  if(target == "48Ti" || target == "Ti48") {
    target_Z = 22;
    target_N = 26;
    target_density = 4.515*g/cm3;
    target_mass = 47.9475*g/mole;
    target_thickness = 2.20*um;
    target_radius = 0.5*cm;
  }
  else if(target == "208Pb" || target == "Pb208") {
    target_Z = 82;
    target_N = 126;
    target_density = 11.382*g/cm3;
    target_mass = 207.97665*g/mole;
    target_thickness = 882*nm;
    target_radius = 0.5*cm;
  }
  else if(target == "196Pt" || target == "Pt196") {
    target_Z = 78;
    target_N = 118;
    target_density = 21.547*g/cm3;
    target_mass = 195.9650*g/mole;
    target_thickness = 738*nm;
    target_radius = 0.5*cm;
  }

  PrintTarget();
  
  return;
}

void Detector_Construction::PrintTarget() {

  G4cout << "\t Z: " << target_Z << "\n\t N: " << target_N << "\n\t Mass: " << G4BestUnit(target_mass,"Mass")
	 <<  "\n\t Density: " << G4BestUnit(target_density,"Volumic Mass") << "\n\t Thickness: "
	 << G4BestUnit(target_thickness,"Length") << "\n\t Radius: " << G4BestUnit(target_radius,"Length")
	 << G4endl;

  return;
}
