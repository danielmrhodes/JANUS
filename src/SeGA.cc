#include "SeGA.hh"
#include "G4Transform3D.hh"
#include "G4AssemblyVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"

SeGA::SeGA(G4bool make_sensitive) {
  
  HpGe = new G4Material("HpGe",32.0,72.61*g/mole,5.323*g/cm3);
  Al =   new G4Material("Al",13.0,26.982*g/mole,2.70*g/cm3); //LR
  //preampMat = new G4Material("preampMat", 13.0,26.982*g/mole,1.35*g/cm3); //LR (Air,copper,and aluminum?)
  //vacuum = new G4Material("vacuum",1,1*g/mole,1e-5*g/cm3);

  // crystal dimensions
  Length = 4.025*cm; //LR;
  outerRadius = 3.165*cm; //LR (Lew)
  fingerRadius = 0.5*cm;

  // dead layer dimensions
  DLinnerRadius = fingerRadius; //LR
  DLouterRadius = DLinnerRadius + 0.03*cm; //LR (Dirk) //0.03*cm original

  // can dimensions
  iCanOuterRadius = 3.73*cm; //LR
  iCanInnerRadius = iCanOuterRadius - 0.05*cm; //LR
  iCanLength = Length; //LR

  oCanThickness = 0.05*cm;
  oCanOuterRadius = 4.325*cm; //LR (Dirk)
  oCanInnerRadius = oCanOuterRadius - oCanThickness; //LR 
  oCanLength = 10.0*cm;
  oCanOffset.setZ(oCanLength - Length - oCanThickness); //LR (0.7*cm gives overlaps in the 37 degree ring)

  preampRadius = oCanInnerRadius - 50*um; //LR
  preampLength = oCanLength - Length - (0.6/2.0)*cm; //LR
  preampOffset.setZ(oCanOffset.z() + Length);

  neckRadius = 3.33/2.0*cm; //LR (approx. from drawing)
  neckLength = 14.47/2.0*cm; //LR (approx. from drawing)
  neckOffset.setX((oCanOuterRadius + neckLength*sin(45.0*deg) + neckRadius*cos(45.0*deg))); //LR 
  neckOffset.setZ(oCanLength + neckLength*cos(45.0*deg)); //LR (approx. from drawing)

  cryoThickness = 0.4*cm; //LR (guess)
  //cryoOuterRadius = (23.06/2.0)*cm; //LR (approx. from drawing)
  cryoOuterRadius = (23.06/2.0)*cm - 1.1*mm; //DMR try to avoid overlaps
  cryoInnerRadius = cryoOuterRadius - cryoThickness; //LR (approx. from drawing)
  cryoBaseThickness = (1.5/2.0)*cm; //LR (guess)
  cryoBaseOffset.setX(neckOffset.x() + neckLength*sin(45.0*deg) + cryoBaseThickness*sin(45.0*deg)); //LR
  cryoBaseOffset.setZ(neckOffset.z() + neckLength*cos(45.0*deg) + cryoBaseThickness*cos(45.0*deg)); //LR

  cryoLength = 34.60/2.0*cm; //LR (approx. from drawing)
  cryoOffset.setX(cryoBaseOffset.x() + cryoBaseThickness*sin(45.0*deg) + cryoLength*sin(45.0*deg)); //LR
  cryoOffset.setZ(cryoBaseOffset.z() + cryoBaseThickness*cos(45.0*deg) + cryoLength*cos(45.0*deg)); //LR

  phisegs=4;
  zsegs=8;

  //Visualization attributes
  Vis1 = new G4VisAttributes(G4Colour::Green());
  Vis1->SetVisibility(true);
  Vis1->SetForceSolid(true);

  Vis2 = new G4VisAttributes(G4Colour::Blue());
  Vis2->SetVisibility(true);
  Vis2->SetForceSolid(true);

  Vis3 = new G4VisAttributes(G4Colour::Magenta());
  Vis3->SetVisibility(true);
  Vis3->SetForceSolid(true);

  G4Colour transGrey (0.8,0.8,0.8,0.3);
  Vis4 = new G4VisAttributes(transGrey);
  Vis4->SetVisibility(true);
  Vis4->SetForceSolid(false);

  G4Colour greyGreen (0.5,0.6,0.65,0.1);
  Vis5 = new G4VisAttributes(greyGreen);
  Vis5->SetVisibility(true);
  Vis5->SetForceSolid(false);

  //Sensitive Detector
  TrackerGamma = NULL;
  if(make_sensitive) {
    TrackerGamma = new GammaSD("GammaTracker");
    G4SDManager::GetSDMpointer()->AddNewDetector(TrackerGamma);
  }

}

SeGA::~SeGA() {}

void SeGA::PlaceDetector(G4LogicalVolume* expHall_log, G4int detNum, G4double Zoffset) { 

  // Material surrounding the crystal
  G4Tubs* iCan = new G4Tubs("iCan",iCanInnerRadius,iCanOuterRadius,iCanLength,0*deg,360*deg); //LR
  G4LogicalVolume* iCan_log = new G4LogicalVolume(iCan,Al,"iCan_log",0,0,0); //LR

  G4Tubs* oCanFull = new G4Tubs("oCanFull",0.0,oCanOuterRadius,oCanLength,0*deg,360*deg); //LR
  G4Tubs* oCanHollow = new G4Tubs("oCanHollow",0.0,oCanInnerRadius,oCanLength,0*deg,360*deg); //LR
  
  G4SubtractionSolid* oCan = new G4SubtractionSolid("oCan",oCanFull,oCanHollow,0,G4ThreeVector(0.*cm,0.*cm,oCanThickness));
  G4LogicalVolume* oCan_log = new G4LogicalVolume(oCan,Al,"oCan_log",0,0,0); //LR

  G4Tubs* preamp = new G4Tubs("preamp",0.,preampRadius,preampLength,0*deg,360*deg); //LR
  G4LogicalVolume* preamp_log = new G4LogicalVolume(preamp,Al,"preamp_log",0,0,0); //LR

  G4Tubs* neck = new G4Tubs("neck",0.,neckRadius,neckLength,0*deg,360*deg); //LR
  G4LogicalVolume* neck_log = new G4LogicalVolume(neck,Al,"neck_log",0,0,0); //LR

  G4Tubs* cryoBase = new G4Tubs("cryoBase",0.,cryoOuterRadius,cryoBaseThickness,0*deg,360*deg); //LR
  G4LogicalVolume* cryoBase_log = new G4LogicalVolume(cryoBase,Al,"cryoBase_log",0,0,0); //LR

  G4Tubs* cryoFull = new G4Tubs("cryoFull",0.0*cm,cryoOuterRadius,cryoLength,0*deg,360*deg); //LR
  G4Tubs* cryoHollow = new G4Tubs("cryoHollow",0.0*cm,cryoInnerRadius,cryoLength-cryoThickness,0*deg,360*deg); //LR
  
  G4SubtractionSolid* cryo = new G4SubtractionSolid("cryo",cryoFull,cryoHollow);
  G4LogicalVolume* cryo_log = new G4LogicalVolume(cryo,Al,"cryo_log",0,0,0); //LR

  G4double rd = 12.975*cm;
  G4double phid = detNum*(360/8)*deg;
  G4double zd = Length + 2*oCanThickness + 0.6*cm;

  if(detNum > 7) {
    zd*=-1;
  }

  G4bool check = false;
  
  G4int seg = 1;
  for(G4int i=0;i<phisegs;i++) {
    for(G4int j=0;j<zsegs;j++) {

      G4double Inner;
      if(j==7) {
        Inner=0.0*cm;
      }
      else {
	Inner=DLouterRadius;
      }
      
      G4Tubs* seg_Tub = new G4Tubs("segTub",Inner,outerRadius,Length/zsegs,360/phisegs*i*deg,360/phisegs*deg);
      G4LogicalVolume* seg_log = new G4LogicalVolume(seg_Tub,HpGe,"segLog",0,TrackerGamma);

      if(j%2) {
	if(i%2) {
          seg_log->SetVisAttributes(Vis1);
	}
	else {
	  seg_log->SetVisAttributes(Vis2);
	}
      }
      else {
	if(i%2) {
          seg_log->SetVisAttributes(Vis2);
	}
	else {
	  seg_log->SetVisAttributes(Vis1);
	}
      }
     
      G4int Copy = (detNum+1)*100 + seg;	

      G4double zshift = j*2*Length/zsegs - double(zsegs-1)/double(zsegs)*Length;
      G4ThreeVector vec(rd*cos(phid),rd*sin(phid),zd+zshift+Zoffset);

      new G4PVPlacement(0,vec,seg_log,"SeGA",expHall_log,false,Copy,check);

      seg++;
	
    }
  }

  if(DLouterRadius > fingerRadius) {
    G4Tubs* DL = new G4Tubs("DL",DLinnerRadius,DLouterRadius,Length*double(zsegs-1)/double(zsegs),0*deg,360*deg);
    G4LogicalVolume* DL_log = new G4LogicalVolume(DL,HpGe,"DL_log",0,0,0);
    DL_log->SetVisAttributes(Vis3);
  
    G4ThreeVector vec1(rd*cos(phid),rd*sin(phid),zd - Length/double(zsegs) + Zoffset);
    new G4PVPlacement(0,vec1,DL_log,"SeGA_DL",expHall_log,false,detNum,check);
  }
  
  oCan_log->SetVisAttributes(Vis4);
  iCan_log->SetVisAttributes(Vis4);
  preamp_log->SetVisAttributes(Vis4);
  neck_log->SetVisAttributes(Vis4);
  cryoBase_log->SetVisAttributes(Vis5);
  cryo_log->SetVisAttributes(Vis5);

  // Place Detector
  G4RotationMatrix assmRot = G4RotationMatrix::IDENTITY;
  G4RotationMatrix cryoRot = G4RotationMatrix::IDENTITY; //LR
  cryoRot.rotateY(45.*deg); //LR

  G4ThreeVector assmPos;
  G4Transform3D assmTrans = G4Transform3D(assmRot,assmPos);
  G4Transform3D oCanTrans = G4Transform3D(assmRot,oCanOffset);
  G4Transform3D preampTrans = G4Transform3D(assmRot,preampOffset);
  G4Transform3D neckTrans = G4Transform3D(cryoRot,neckOffset);
  G4Transform3D cryoBaseTrans = G4Transform3D(cryoRot,cryoBaseOffset);
  G4Transform3D cryoTrans = G4Transform3D(cryoRot,cryoOffset);

  G4AssemblyVolume* SeGA_assembly = new G4AssemblyVolume();
  SeGA_assembly->AddPlacedVolume(iCan_log,assmTrans);
  SeGA_assembly->AddPlacedVolume(oCan_log,oCanTrans);
  SeGA_assembly->AddPlacedVolume(preamp_log,preampTrans);
  SeGA_assembly->AddPlacedVolume(neck_log,neckTrans);
  SeGA_assembly->AddPlacedVolume(cryoBase_log,cryoBaseTrans);
  SeGA_assembly->AddPlacedVolume(cryo_log,cryoTrans);
  
  if(detNum>7) {
    assmRot.rotateX(180*deg);
  }
  assmRot.rotateZ(phid);
  
  assmPos.setX(rd*cos(phid));
  assmPos.setY(rd*sin(phid));
  assmPos.setZ(zd+Zoffset);
  
  G4Transform3D SeGATrans = G4Transform3D(assmRot,assmPos);
  SeGA_assembly->MakeImprint(expHall_log,SeGATrans,detNum,check);

}

void SeGA::Placement(G4LogicalVolume* expHall_log, G4double ZOff) {

  for(G4int i=0;i<16;i++) {
    PlaceDetector(expHall_log,i,ZOff); 
  }
  
  return;

}
