#include <iostream>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTree.h"

#include "Data_Format.hh"

inline void print_data(const SegaData& data) {

  std::cout << "\n  Detector " << data.det << " Segment " << data.seg << "\n  " << data.en
	    << " keV deposited at (" << data.x << "," << data.y << "," << data.z
	    << ")\n";
  
}

inline void print_data(const Bambino2Data& data) {
  
  std::cout << "\n  Detector " << data.det << " Ring " << data.ring
	    << " Sector " << data.sector << "\n  " << data.en
	    << " MeV deposited at (" << data.x << "," << data.y << "," << data.z
	    << ")\n  Proj(" << data.proj << ") Rec("
	    << data.rec << " )\n";
  
}

void print_data(const Header& head, const JANUSData& data) {

  std::cout << "\n---------------------------------\n";
  std::cout << "JANUS Event " << head.evtNum << "\n";
  std::cout << head.nBdata << " Bambino2 Events\n";
  std::cout << head.nSdata << " SeGA Events\n\n";

  if(head.nBdata) {
    std::cout << "Bambino2 Info:\n";
    for(int i=0;i<head.nBdata;i++) {

      std::cout << " Entry " << i;
      print_data(data.bData[i]);
    
    }
  }

  if(head.nSdata) {
    std::cout << "SeGA Info:\n";
    for(int i=0;i<head.nSdata;i++) {

      std::cout << " Entry " << i;
      print_data(data.sData[i]);
    
    }
  }
  std::cout << "---------------------------------\n" << std::flush;
  
}

int main(int argc, char** argv) {
  
  if(argc < 3) {
    std::cerr << "Usage: unpacker INPUT_FILE OUTPUT_FILE" << std::endl;
    return 1;
  }

  const char* input_filename = argv[1];
  const char* output_filename = argv[2];

  FILE* input_file = fopen(input_filename,"rb");

  TH1* num = new TH1D("EvtNum","Event Number",5000,0,1000000);

  //Bambino2 Singles
  TH1* rEn = new TH1D("Ring_Energy","Janus Ring Energy",10000,0,1000);
  TH1* sEn = new TH1D("Sector_Energy","Janus Sector Energy",10000,0,1000);

  TH2* bSum = new TH2D("Summary","Janus Summary",122,0,122,1000,0,1000);
  
  TH2* kin = new TH2D("Kinematics","Kinematic Curve",180,0,180,1000,0,1000);
  TH2* pKin = new TH2D("Projectile_Kinematics","Projectile Kinematic Curve",180,0,180,1000,0,1000);
  TH2* rKin = new TH2D("Recoil_Kinematics","Recoil Kinematic Curve",180,0,180,1000,0,1000);
  TH2* bKin = new TH2D("Both_Kinematics","Both Kinematic Curve",180,0,180,1000,0,1000);
  
  TH2* pid1 = new TH2D("Pid_Det1","DS PID",26,0,26,1000,0,1000);
  TH2* sec1 = new TH2D("Sectors_Det1","DS Sectors",34,0,34,1000,0,1000);
  
  TH2* pid0 = new TH2D("Pid_Det0","US PID",26,0,26,1000,0,1000);
  TH2* sec0 = new TH2D("Sectors_Det0","US Sectors",34,0,34,1000,0,1000);

  //SeGA Singles
  TH1* sEnergy = new TH1D("Segment_Energy","SeGA Segment Energy",3000,0,3000);
  TH1* cEnergy = new TH1D("Core_Energy","SeGA Core Energy",3000,0,3000);
  TH1* cEnergy_FEP = new TH1D("Core_Energy_FEP","SeGA Core Energy FEP",3000,0,3000);
  TH1* cEnergy_nFEP = new TH1D("Core_Energy_nFEP","SeGA Core Energy Not FEP",3000,0,3000);

  TH2* segSum = new TH2D("Segment_Summary","SeGA Segments",512,1,513,2000,0,4000);
  TH2* coreSum = new TH2D("Core_Summary","SeGA Cores",16,1,17,2000,0,4000);

  Header header;
  JANUSData data;
  while(fread(&header,header.bytes(),1,input_file)) {
  
    const int nB = header.nBdata;
    const int nS = header.nSdata;
    const int nE = header.evtNum;
    num->Fill(nE);
    
    fread(&data.bData,nB*sizeof(Bambino2Data),1,input_file);
    fread(&data.sData,nS*sizeof(SegaData),1,input_file);
    
    //Bambino2 Singles
    for(int i=0;i<nB;i++) {

      double en = data.bData[i].en;
	  
      TVector3 pos(data.bData[i].x,data.bData[i].y,data.bData[i].z);
      double th = pos.Theta()*TMath::RadToDeg();
	
      int det = data.bData[i].det;
      int ring = data.bData[i].ring;
      int sec = data.bData[i].sector;

      bool pj = data.bData[i].proj;
      bool rc = data.bData[i].rec;
	
      if(ring) {

	rEn->Fill(en);
	kin->Fill(th,en);

	if(pj) {
	  pKin->Fill(th,en);
	}
	if(rc) {
	  rKin->Fill(th,en);
	}
	if(pj && rc) {
	  bKin->Fill(th,en);
	}

	if(det==0) {
	  bSum->Fill(ring+32,en);
	  pid0->Fill(ring,en);
	}
	else {
	  bSum->Fill(ring+96,en);
	  pid1->Fill(ring,en);
	}
 
      }
      else {

	sEn->Fill(en);

	if(det==0) {
	  bSum->Fill(sec,en);
	  sec0->Fill(sec,en);
	}
	else {
	  bSum->Fill(sec+64,en);
	  sec1->Fill(sec,en);
	}

      }
	
    } //End Bambino2 Singles

    //SeGA Singles
    for(int i=0;i<nS;i++) {

      int det = data.sData[i].det;
      int seg = data.sData[i].seg;
      double energy = data.sData[i].en;
      bool fep = data.sData[i].fep;

      if(bool(seg)) {
	  
	sEnergy->Fill(energy);
	  
	int num = 32*(det-1) + seg;
	segSum->Fill(num,energy);
	  
      }
      else {

	coreSum->Fill(det,energy);
	cEnergy->Fill(energy);
      }

      if(fep) {
	cEnergy_FEP->Fill(energy);
      }
      else {
	cEnergy_nFEP->Fill(energy);
      }
	
    } //End SeGA Singles
    
  } //End while(true)
   
  fclose(input_file);

  TFile* output_file = new TFile(output_filename,"RECREATE");
  output_file->mkdir("SeGA");
  output_file->mkdir("Bambino2");

  num->Write();

  output_file->cd("SeGA");
  
  sEnergy->Write();
  cEnergy->Write();
  cEnergy_FEP->Write();
  cEnergy_nFEP->Write();
  
  segSum->Write();
  coreSum->Write();

  output_file->cd("Bambino2");
  
  rEn->Write();
  sEn->Write();
  
  kin->Write();
  pKin->Write();
  rKin->Write();
  bKin->Write();
  
  pid0->Write();
  pid1->Write();

  sec0->Write();
  sec1->Write();

  bSum->Write();

  output_file->Close();
  
}

