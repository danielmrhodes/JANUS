#include <iostream>
#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMath.h"

#include "Data_Format.hh"

//These should match the parameters defined in the simulation input
int beamZ = 48;
int beamA = 106;
double beam_mass = 98626.9; // MeV/c^2
  
int targZ = 82;
int targA = 208;
double targ_mass = 193688; // MeV/c^2
//double targ_mass = 44652.;

double beam_en = 450.0; // MeV
//double beam_en = 291.5;

////Kinematics////
double Theta_CM_FP(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += TMath::Pi();
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }
  else {
    return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + TMath::Pi();
  }

}

double Theta_CM_FR(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += TMath::Pi();
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return TMath::Pi() - (std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB);
  }
  else {
    return -std::asin(tau*std::sin(-ThetaLAB)) - ThetaLAB;
  }

}

double Theta_LAB_Max(double Ep, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  if(tau < 1.0) {
    return TMath::Pi();
  }
  else {
    return std::asin(1.0/tau);
  }
  
}

double Theta_LAB(double thetaCM, double Ep, double Ex=0.) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(thetaCM)/(std::cos(thetaCM) + tau);

  if(tanTheta > 0) {
    return std::atan(tanTheta);
  }
  else {
    return std::atan(tanTheta) + TMath::Pi();
  }
  
}

double Recoil_Theta_LAB(double thetaCM, double Ep, double Ex=0.) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(TMath::Pi() - thetaCM)/(std::cos(TMath::Pi() - thetaCM) + tau);
  
  return std::atan(tanTheta);
  
}

double KE_LAB(double thetaCM, double Ep, double Ex=0.) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = std::pow(targ_mass/(beam_mass + targ_mass),2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

double Recoil_KE_LAB(double thetaCM, double Ep, double Ex=0.) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = beam_mass*targ_mass/std::pow(beam_mass + targ_mass,2);
  double term2 = 1 + tau*tau + 2*tau*std::cos(TMath::Pi() - thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}
//////////////////

////Positions////
//Bambino2 segment position
TVector3 GetPos(const int det, const int ring, const int sec) {

   if(det < 0 || det > 1 || ring < 1 || ring > 24 || sec < 1 || sec > 32) {
     std::cout << "Bad det, ring, sec (" << det << "," << ring << "," << sec << ")"<< std::endl;
     return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  const double PI = TMath::Pi();
  
  double phi_offset = 0.5*PI; // Phi of sector 1 of downstream detector 
  bool clockwise; // Winding direction of sectors.
  if(det==0) {
    clockwise = false;
  }
  else {
    clockwise = true;
  }

  double janus_outer_radius = 3.5;
  double janus_inner_radius = 1.1;

  TVector3 pos(1.,0,0);
  double rad_slope = (janus_outer_radius - janus_inner_radius)/24.;
  double rad_offset = janus_inner_radius;
   
  pos.SetPerp((ring - 0.5)*rad_slope + rad_offset);
  
  double phi = phi_offset + (clockwise ? -1 : 1) * 2.*PI/32. * (sec - 1);
  pos.SetPhi(phi);

  double zoff;
  if(det == 1) {
    zoff = 2.6;
  }
  else {
    zoff = 3.4;
  }
  pos.SetZ(zoff);

  if(det == 0) {
    pos.RotateY(PI);
  }

  return pos;
  
}

//SeGA segment position
TVector3 GetPos(const int det, const int seg) {

  if(det < 0 || det > 16 || seg < 1 || seg > 32) {
    std::cout << "Bad det, seg (" << det << "," << seg << ")" << std::endl;
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  const double PI = TMath::Pi();

  int quad = int(seg-1)/(int)8;

  int slice;
  if(seg < 9) {
    slice = seg-1;
  }
  else if(seg < 17) {
    slice = seg-9;
  }
  else if(seg < 25) {
    slice = seg-17;
  }
  else {
    slice = seg-25;
  }

  double length = 4.025;
  double outerRadius = 3.165;
  double innerRadius;
  
  if(slice == 7) {
    innerRadius = 0.0;
  }
  else {
    innerRadius = 0.5 + 0.03; //fingerRadius + DL thickness
  }

  TVector3 origin(1,0,0);
  origin.SetPerp(outerRadius - innerRadius);
  origin.SetPhi((quad+0.5)*2*PI/4.);
  origin.SetZ((length/8.0)*(2.0*slice - 7.0));

  double rd = 12.975;
  double phid = (det-1)*(2.*PI/8.);
  double zd = length + 2*0.05;
  if(det > 8) {
    zd*=-1;
  }
  
  double Zoffset = 3.1;
  TVector3 pos(rd*TMath::Cos(phid),rd*TMath::Sin(phid),zd+Zoffset);

  return origin+pos;

}
/////////////////

struct BAM2 {

  //realistic info
  int det, ring, sector;
  double rEn, sEn;

  //perfect info
  TVector3 rPos, sPos;
  bool rP, rR, sP, sR;
  
};

struct SEGA {

  int det, nsegs, segs[32];
  double cEn, sEn[32];

  int MainSeg() {
    
    int index=0;
    for(int i=0;i<nsegs;i++) {
      if(sEn[i] > sEn[index]) {
	index=i;
      }
    }
    
    return segs[index];
  }
  
};

struct BuiltData {

  int evt;
  int nBa;
  int nSe;
  
  BAM2 bam2[50];
  SEGA sega[16];

  void MakeHit(const Bambino2Data ringHit, const Bambino2Data secHit) {

    bam2[nBa].det = ringHit.det; //either one works
    bam2[nBa].ring = ringHit.ring;
    bam2[nBa].sector = secHit.sector;
    
    bam2[nBa].rEn = ringHit.en;
    bam2[nBa].sEn = secHit.en;

    bam2[nBa].rPos = TVector3(ringHit.x,ringHit.y,ringHit.z);
    bam2[nBa].sPos = TVector3(secHit.x,secHit.y,secHit.z);
    
    bam2[nBa].rP = ringHit.proj;
    bam2[nBa].rR = ringHit.rec;

    bam2[nBa].sP = secHit.proj;
    bam2[nBa].sR = secHit.rec;

    nBa++;
  }
  
};

BuiltData BuildData(const Header& head, const JANUSData& dat) {

  const int nBch = head.nBdata;
  const int nSch = head.nSdata;
  const int nE = head.evtNum;
  
  BuiltData data;
  data.evt = nE;
  data.nBa = 0;
  data.nSe = 0;

  //Correlate Bambino2 Data
  std::vector<Bambino2Data> rings;
  std::vector<Bambino2Data> sectors;
  for(int i=0;i<nBch;i++) {

    Bambino2Data bChan(dat.bData[i]);
    
    if(bChan.IsRing()) {
      rings.push_back(bChan);
    }
    else {
      sectors.push_back(bChan);
    }
  }

  std::sort(rings.begin(),rings.end());
  std::sort(sectors.begin(),sectors.end());

  std::vector<bool> used_rings;
  std::vector<bool> used_sectors;
  
  used_rings.resize(rings.size());
  std::fill(used_rings.begin(),used_rings.end(),false);

  used_sectors.resize(sectors.size());
  std::fill(used_sectors.begin(),used_sectors.end(),false);

  for(unsigned int i=0;i<sectors.size();i++) {
    if(used_sectors.at(i)) {
      continue;
    }
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j)) {
        continue;
      }

         //Same detector
      if((sectors.at(i).det == rings.at(j).det) &&

	 //Same energy
	 (TMath::Abs(sectors.at(i).en - rings.at(j).en) < 1.)) {

	data.MakeHit(rings.at(j),sectors.at(i));

	used_sectors.at(i) = true;
	used_rings.at(j) = true;
	break;

      }
    }
  }

  bool broken = false;
  for(unsigned int i=0;i<sectors.size();i++) {
    broken=false;
    if(used_sectors.at(i)) {
      continue;
    }
    
    for(unsigned int j=0;j<rings.size();j++) {
      if(used_rings.at(j)) {
        continue;
      }

      for(unsigned int k=0;k<sectors.size();k++) {
	if(used_sectors.at(k)) {
	  continue;
	}

	  //Same detector
	if((sectors.at(i).det == rings.at(j).det) && (sectors.at(k).det == rings.at(j).det) &&

	   //Sector energies add to ring energy
	   (TMath::Abs(sectors.at(i).en + sectors.at(k).en - rings.at(j).en) < 1.)) {

	  data.MakeHit(rings.at(j),sectors.at(i));
	  data.MakeHit(rings.at(j),sectors.at(k));
	  
	  used_sectors.at(i) = true;
	  used_rings.at(j) = true;
	  used_sectors.at(k) = true;
	  broken = true;
	  break;

	}
      }
      if(broken) {
	break;
      }
    }
  }

  //Organize SeGA data
  std::vector<bool> exists;
  exists.resize(16);
  std::fill(exists.begin(),exists.end(),false);

  int nS = 0;
  for(int i=0;i<nSch;i++) {

    int detect = dat.sData[i].det;
    int segment = dat.sData[i].seg;
    double energy = dat.sData[i].en;
    
    if(!exists.at(detect-1)) {
      data.sega[nS].det = detect;

      if((bool)segment) {
        data.sega[nS].nsegs = 1;
	data.sega[nS].segs[0] = segment;
	data.sega[nS].sEn[0] = energy;
      }
      else {
	data.sega[nS].nsegs = 0;
	data.sega[nS].cEn = energy;
      }

      nS++;
      exists.at(detect-1) = true;
    }
    else {

      int index = 0;
      for(int j=0;j<nS;j++) {
	if(data.sega[j].det == detect) {
	  index = j;
	  break;
	}
      }
      
      int Nsegs = data.sega[index].nsegs;
      
      if((bool)segment) {
	data.sega[index].segs[Nsegs] = segment;
	data.sega[index].sEn[Nsegs] = energy;
	data.sega[index].nsegs++;
      }
      else {
	data.sega[index].cEn = energy;
      }
      
    }
    
  }
  data.nSe = nS;

  return data;
}

int main(int argc, char** argv) {
  
  if(argc < 3) {
    std::cerr << "Usage: correlator INPUT_FILE OUTPUT_FILE" << std::endl;
    return 1;
  }

  const char* input_filename = argv[1];
  const char* output_filename = argv[2];

  FILE* input_file = fopen(input_filename,"rb");
  
  //Bambino2 singles
  TH2* bSum = new TH2D("Summary","Janus Summary",120,1,121,500,0,500);
  
  TH1* secD0 = new TH1D("Sectors_Det0","US Sectors",32,1,33);
  TH1* secD1 = new TH1D("Sectors_Det1","DS Sectors",32,1,33);
  
  TH2* rPid0 = new TH2D("RingPID_Det0","US RingEn PID",24,1,25,500,0,500);
  TH2* rPid1 = new TH2D("RingPID_Det1","DS RingEn PID",24,1,25,500,0,500);
  TH2* rPid1m1 = new TH2D("RingPID_Det1m1","DS RingEn PID Mult1",24,1,25,500,0,500);
  TH2* rPid1m2 = new TH2D("RingPID_Det1m2","DS RingEn PID Mult2",24,1,25,500,0,500);

  TH2* sPid0 = new TH2D("SecPID_Det0","US SectorEn PID",24,1,25,500,0,500);
  TH2* sPid1 = new TH2D("SecPID_Det1","DS SectorEn PID",24,1,25,500,0,500);
  TH2* sPid1m1 = new TH2D("SecPID_Det1m1","DS SectorEn PID Mult1",24,1,25,500,0,500);
  TH2* sPid1m2 = new TH2D("SecPID_Det1m2","DS SectorEn PID Mult2",24,1,25,500,0,500);

  TH2* sPid1_p = new TH2D("SecPID_Det1_proj","DS SectorEn Projectile PID",24,1,25,500,0,500);
  TH2* sPid1_r = new TH2D("SecPID_Det1_rec","DS SectorEn Recoil PID",24,1,25,500,0,500);
  
  //SeGA singles
  TH1* coreEnergy = new TH1D("Core_Energy","SeGA Core Energy",3000,0,3000);
  TH1* segEnergy = new TH1D("Seg_Energy","SeGA Segment Energy",3000,0,3000);

  TH2* coreSum = new TH2D("Core_Summary","Core Energy Summary",16,1,17,3000,0,3000);
  TH2* segSum = new TH2D("Seg_Summary","Segment Energy Summary",512,1,513,3000,0,3000);

  //Coincidences
  //Projectile DS
  TH2* sPidDS = new TH2D("SecPID_DS","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidDS = new TH2D("RingPID_DS","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* pCoreEnergyDS = new TH1D("Core_EnergyDS","SeGA Core Energy",3000,0,3000);
  TH2* pCoreSumDS = new TH2D("Core_SummaryDS","Core Energy Summary",16,1,17,3000,0,3000);
  
  TH1* pDopEnergyDS = new TH1D("Dop_EnergyDS","Doppler Energy",6000,0,3000);
  TH2* pDopSumDS = new TH2D("Dop_SummaryDS","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH2* pDopvPartDS = new TH2D("DopEn_v_PartEn_DS","Doppler Energy vs Particle Energy",3000,0,3000,500,0,500);
  TH2* pDopvPartNS2DS = new TH2D("DopEn_v_PartEn_NoS2_DS","Doppler Energy vs Particle Energy (No Sol2)",3000,0,3000,500,0,500);

  TH2* pThCorDS = new TH2D("Theta_CorrDS","Theta Correlation",3000,0,3000,90,0,180);
  TH2* pThCrtDS = new TH2D("Theta_CrctDS","Theta Correction",6000,0,3000,90,0,180);

  double thing1 = 65.*180./32.0;
  TH2* pPhCorDS = new TH2D("Phi_CorrDS","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pPhCrtDS = new TH2D("Phi_CrctDS","Phi Correction",6000,0,3000,32,0,thing1);

  TH1* pReconEnergyDS = new TH1D("Recon_EnergyDS","Recon Energy",6000,0,3000);
  TH2* pReconSumDS = new TH2D("Recon_SummaryDS","Recon Energy Summary",16,1,17,6000,0,3000);

  TH2* pReconvPartDS = new TH2D("ReconEn_v_partEn_DS","Recon Energy vs Particle Energy",3000,0,3000,500,0,500);
  TH2* pReconvPartNS2DS = new TH2D("ReconEn_v_partEn_NoS2_DS","Recon Energy vs Particle Energy (No Sol2)",3000,0,3000,500,0,500);

  TH2* pReconThCorDS = new TH2D("ReconTheta_CorrDS","Recon Theta Correlation",3000,0,3000,90,0,180);
  TH2* pReconThCrtDS = new TH2D("ReconTheta_CrctDS","Recon Theta Correction",6000,0,3000,90,0,180);

  TH2* pReconPhCorDS = new TH2D("ReconPhi_CorrDS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pReconPhCrtDS = new TH2D("ReconPhi_CrctDS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Projectile US
  TH2* sPidUS = new TH2D("SecPID_US","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidUS = new TH2D("RingPID_US","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* pCoreEnergyUS = new TH1D("Core_EnergyUS","SeGA Core Energy",3000,0,3000);
  TH2* pCoreSumUS = new TH2D("Core_SummaryUS","Core Energy Summary",16,1,17,3000,0,3000);
  
  TH1* pDopEnergyUS = new TH1D("Dop_EnergyUS","Doppler Energy",6000,0,3000);
  TH2* pDopSumUS = new TH2D("Dop_SummaryUS","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH2* pDopvPartUS = new TH2D("Dop_PartEn_US","Doppler Energy vs Particle Energy",3000,0,3000,500,0,500);

  TH2* pThCorUS = new TH2D("Theta_CorrUS","Theta Correlation",3000,0,3000,90,0,180);
  TH2* pThCrtUS = new TH2D("Theta_CrctUS","Theta Correction",6000,0,3000,90,0,180);
  
  TH2* pPhCorUS = new TH2D("Phi_CorrUS","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pPhCrtUS = new TH2D("Phi_CrctUS","Phi Correction",6000,0,3000,32,0,thing1);

  //Recoil
  TH2* sPidRec = new TH2D("SecPID_Rec","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidRec = new TH2D("RingPID_Rec","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* rCoreEnergy = new TH1D("Core_EnergyRec","SeGA Core Energy",3000,0,3000);
  TH2* rCoreSum = new TH2D("Core_SummaryRec","Core Energy Summary",16,1,17,3000,0,3000);
  
  TH1* rDopEnergy = new TH1D("Dop_EnergyRec","Doppler Energy",6000,0,3000);
  TH2* rDopSum = new TH2D("Dop_SummaryRec","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH2* rDopvPart = new TH2D("DopEv_v_PartEn_Rec","Doppler Energy vs Particle Energy",3000,0,3000,500,0,500);
  //TH2* rDopvPartNS2 = new TH2D("DopEv_v_PartEn_NoS2_Rec","Doppler Energy vs Particle Energy (No Sol2)",3000,0,3000,500,0,500);

  TH2* rThCor = new TH2D("Theta_CorrRec","Theta Correlation",3000,0,3000,90,0,180);
  TH2* rThCrt = new TH2D("Theta_CrctRec","Theta Correction",6000,0,3000,90,0,180);;
  TH2* rPhCor = new TH2D("Phi_CorrRec","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* rPhCrt = new TH2D("Phi_CrctRec","Phi Correction",6000,0,3000,32,0,thing1);

  TH1* rReconEnergy = new TH1D("Recon_EnergyRec","Recon Energy",6000,0,3000);
  TH2* rReconSum = new TH2D("Recon_SummaryRec","Recon Energy Summary",16,1,17,6000,0,3000);

  TH2* rReconvPart = new TH2D("ReconEn_v_PartEn_Rec","Recon Energy vs Particle Energy",3000,0,3000,500,0,500);
  //TH2* rReconvPartNS2 = new TH2D("ReconEv_v_PartEn_NoS2_Rec","Recon Energy vs Particle Energy (No Sol2)",3000,0,3000,500,0,500);

  TH2* rReconThCor = new TH2D("ReconTheta_CorrRec","Recon Theta Correlation",3000,0,3000,90,0,180);
  TH2* rReconThCrt = new TH2D("ReconTheta_CrctRec","Recon Theta Correction",6000,0,3000,90,0,180);

  TH2* rReconPhCor = new TH2D("ReconPhi_CorrRec","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* rReconPhCrt = new TH2D("ReconPhi_CrctRec","Recon Phi Correction",6000,0,3000,32,0,thing1);

  /*
  //SeGA dets
  std::vector<TH1*> pDetDopEnDS;
  std::vector<TH2*> pDetDopSumDS;
  std::vector<TH2*> pDetPhCorDS;
  std::vector<TH2*> pDetPhCrtDS;
  std::vector<TH2*> pDetThCorDS;
  std::vector<TH2*> pDetThCrtDS;

  std::vector<TH1*> pDetDopEnUS;
  std::vector<TH2*> pDetDopSumUS;
  std::vector<TH2*> pDetPhCorUS;
  std::vector<TH2*> pDetPhCrtUS;
  std::vector<TH2*> pDetThCorUS;
  std::vector<TH2*> pDetThCrtUS;

  std::vector<TH1*> rDetDopEn;
  std::vector<TH2*> rDetDopSum;
  std::vector<TH2*> rDetPhCor;
  std::vector<TH2*> rDetPhCrt;
  std::vector<TH2*> rDetThCor;
  std::vector<TH2*> rDetThCrt;
  
  for(int i=0;i<16;i++) {

    //Downstream
    pDetDopEnDS.push_back(new TH1D(Form("Dop_EnergyDS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				6000,0,3000));

    pDetDopSumDS.push_back(new TH2D(Form("Dop_SegSumDS_Det%02i",i+1),Form("Det%02i Doppler Seg Summary",i+1),
				 32,1,33,6000,0,3000));

    pDetPhCorDS.push_back(new TH2D(Form("PhiCorDS_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				3000,0,3000,32,0,thing1));

    pDetPhCrtDS.push_back(new TH2D(Form("PhiCrtDS_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				6000,0,3000,32,0,thing1));

    pDetThCorDS.push_back(new TH2D(Form("ThetaCorDS_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				3000,0,3000,90,0,180));

    pDetThCrtDS.push_back(new TH2D(Form("ThetaCrtDS_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				6000,0,3000,90,0,180));

    //Upstream
    pDetDopEnUS.push_back(new TH1D(Form("Dop_EnergyUS_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				   6000,0,3000));

    pDetDopSumUS.push_back(new TH2D(Form("Dop_SegSumUS_Det%02i",i+1),Form("Det%02i Doppler Seg Summary",i+1),
				    32,1,33,6000,0,3000));

    pDetPhCorUS.push_back(new TH2D(Form("PhiCorUS_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				   3000,0,3000,32,0,thing1));

    pDetPhCrtUS.push_back(new TH2D(Form("PhiCrtUS_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				   6000,0,3000,32,0,thing1));

    pDetThCorUS.push_back(new TH2D(Form("ThetaCorUS_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				   3000,0,3000,90,0,180));

    pDetThCrtUS.push_back(new TH2D(Form("ThetaCrtUS_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				   6000,0,3000,90,0,180));

    //Recoil
    rDetDopEn.push_back(new TH1D(Form("Dop_EnergyRec_Det%02i",i+1),Form("Det%02i Gamma Energy",i+1),
				   6000,0,3000));

    rDetDopSum.push_back(new TH2D(Form("Dop_SegSumRec_Det%02i",i+1),Form("Det%02i Doppler Seg Summary",i+1),
				    32,1,33,6000,0,3000));

    rDetPhCor.push_back(new TH2D(Form("PhiCorRec_Det%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				   3000,0,3000,32,0,thing1));

    rDetPhCrt.push_back(new TH2D(Form("PhiCrtRec_Det%02i",i+1),Form("Det%02i Phi Correction",i+1),
				   6000,0,3000,32,0,thing1));

    rDetThCor.push_back(new TH2D(Form("ThetaCorRec_Det%02i",i+1),Form("Det%02i Theta Correlation",i+1),
				   3000,0,3000,90,0,180));

    rDetThCrt.push_back(new TH2D(Form("ThetaCrtRec_Det%02i",i+1),Form("Det%02i Theta Correction",i+1),
				   6000,0,3000,90,0,180));
    
  }
  */

  //Bambino2 rings
  std::vector<TH1*> pRingCoreEnDS;
  std::vector<TH1*> pRingDopEnDS;
  std::vector<TH1*> pRingRecEnDS;
  //std::vector<TH2*> pRingPhCorDS;
  //std::vector<TH2*> pRingPhCrtDS;
  //std::vector<TH2*> pRingThCorDS;
  //std::vector<TH2*> pRingThCrtDS;

  std::vector<TH1*> pRingCoreEnUS;
  std::vector<TH1*> pRingDopEnUS;
  //std::vector<TH2*> pRingPhCorUS;
  //std::vector<TH2*> pRingPhCrtUS;
  //std::vector<TH2*> pRingThCorUS;
  //std::vector<TH2*> pRingThCrtUS;

  std::vector<TH1*> rRingCoreEn;
  std::vector<TH1*> rRingDopEn;
  std::vector<TH1*> rRingRecEn;
  //std::vector<TH2*> rRingPhCor;
  //std::vector<TH2*> rRingPhCrt;
  //std::vector<TH2*> rRingThCor;
  //std::vector<TH2*> rRingThCrt;
  std::vector<TH2*> rRingRecThCor;
  std::vector<TH2*> rRingRecThCrt;
  
  for(int i=0;i<24;i++) {

    //Downstream
    pRingCoreEnDS.push_back(new TH1D(Form("Core_EnergyDS_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				     3000,0,3000));

    pRingDopEnDS.push_back(new TH1D(Form("Dop_EnergyDS_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				    6000,0,3000));

    pRingRecEnDS.push_back(new TH1D(Form("Rec_EnergyDS_R%02i",i+1),Form("Ring%02i Recon Energy",i+1),
				    6000,0,3000));

    /*
    pRingPhCorDS.push_back(new TH2D(Form("PhiCorDS_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				3000,0,3000,32,0,thing1));

    pRingPhCrtDS.push_back(new TH2D(Form("PhiCrtDS_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				6000,0,3000,32,0,thing1));

    pRingThCorDS.push_back(new TH2D(Form("ThetaCorDS_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				3000,0,3000,90,0,180));

    pRingThCrtDS.push_back(new TH2D(Form("ThetaCrtDS_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				6000,0,3000,90,0,180));
    */

    //Upstream
    pRingCoreEnUS.push_back(new TH1D(Form("Core_EnergyUS_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				     3000,0,3000));

    pRingDopEnUS.push_back(new TH1D(Form("Dop_EnergyUS_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				    6000,0,3000));

    /*
    pRingPhCorUS.push_back(new TH2D(Form("PhiCorUS_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				    3000,0,3000,32,0,thing1));

    pRingPhCrtUS.push_back(new TH2D(Form("PhiCrtUS_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				    6000,0,3000,32,0,thing1));

    pRingThCorUS.push_back(new TH2D(Form("ThetaCorUS_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				    3000,0,3000,90,0,180));

    pRingThCrtUS.push_back(new TH2D(Form("ThetaCrtUS_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				    6000,0,3000,90,0,180));
    */

    //Recoil
    rRingCoreEn.push_back(new TH1D(Form("Core_EnergyRec_R%02i",i+1),Form("Ring%02i Core Energy",i+1),
				   3000,0,3000));

    rRingDopEn.push_back(new TH1D(Form("Dop_EnergyRec_R%02i",i+1),Form("Ring%02i Doppler Energy",i+1),
				  6000,0,3000));

    rRingRecEn.push_back(new TH1D(Form("Rec_EnergyRec_R%02i",i+1),Form("Ring%02i Recon Energy",i+1),
				  6000,0,3000));

    /*
    rRingPhCor.push_back(new TH2D(Form("PhiCorRec_R%02i",i+1),Form("Det%02i Phi Correlation",i+1),
				    3000,0,3000,32,0,thing1));

    rRingPhCrt.push_back(new TH2D(Form("PhiCrtRec_R%02i",i+1),Form("Ring%02i Phi Correction",i+1),
				    6000,0,3000,32,0,thing1));

    rRingThCor.push_back(new TH2D(Form("ThetaCorRec_R%02i",i+1),Form("Ring%02i Theta Correlation",i+1),
				    3000,0,3000,90,0,180));

    rRingThCrt.push_back(new TH2D(Form("ThetaCrtRec_R%02i",i+1),Form("Ring%02i Theta Correction",i+1),
				    6000,0,3000,90,0,180));
    */

    rRingRecThCor.push_back(new TH2D(Form("rRecThetaCor_R%02i",i+1),Form("Ring%02i Recon Theta Correlation",i+1),
				     3000,0,3000,90,0,180));

    rRingRecThCrt.push_back(new TH2D(Form("rRecThetaCrt_R%02i",i+1),Form("Ring%02i Recon Theta Correction",i+1),
				     6000,0,3000,90,0,180));
  }

  std::cout << "Correlating and histograming data..." << std::endl;
  
  TVector3 incBeam = TVector3(0.0,0.0,1.0);
  double Sol2_En = KE_LAB(Theta_CM_FP(Theta_LAB_Max(beam_en),beam_en),beam_en);
  
  Header header;
  JANUSData jData;
  while(fread(&header,header.bytes(),1,input_file)) {

    const int nB = header.nBdata;
    const int nS = header.nSdata;
    //const int nE = header.evtNum;
  
    fread(&jData.bData,nB*sizeof(Bambino2Data),1,input_file);
    fread(&jData.sData,nS*sizeof(SegaData),1,input_file);

    BuiltData data = BuildData(header,jData);

    //Bambino2 singles
    for(int i=0;i<data.nBa;i++) {

      int det = data.bam2[i].det;
      int ring = data.bam2[i].ring;
      int sec = data.bam2[i].sector;

      double ring_en = data.bam2[i].rEn;
      double sec_en = data.bam2[i].sEn;

      TVector3 segPos = GetPos(det,ring,sec);
      
      if(!det) {

	bSum->Fill(sec,sec_en);
	bSum->Fill(ring+32,ring_en);
	
	rPid0->Fill(ring,ring_en);
	sPid0->Fill(ring,sec_en);
	secD0->Fill(sec);

      }
      else {
	
	rPid1->Fill(ring,ring_en);
	sPid1->Fill(ring,sec_en);
	secD1->Fill(sec);

	bSum->Fill(sec+64,sec_en);
	bSum->Fill(ring+96,ring_en);

	if(data.bam2[i].rP && data.bam2[i].sP) {
	  sPid1_p->Fill(ring,sec_en);
	}

	if(data.bam2[i].rR && data.bam2[i].sR) {
	  sPid1_r->Fill(ring,sec_en);
	}

	if(data.nBa == 1) {
	  rPid1m1->Fill(ring,ring_en);
	  sPid1m1->Fill(ring,sec_en);
	}
	else if(data.nBa == 2) {
	  rPid1m2->Fill(ring,ring_en);
	  sPid1m2->Fill(ring,sec_en);
	}

      }
      
    } //End Bambino2 singles

    //SeGA singles
    for(int i=0;i<data.nSe;i++) {

      int det = data.sega[i].det;
      double core_en = data.sega[i].cEn;
      
      coreEnergy->Fill(core_en);
      coreSum->Fill(det,core_en);

      for(int j=0;j<data.sega[i].nsegs;j++) {
	
	double seg_en = data.sega[i].sEn[j];
	int num = data.sega[i].segs[j] + 32*(det-1);
	
	segEnergy->Fill(seg_en);
	segSum->Fill(num,seg_en);
	
      }
      
    } //End SeGA singles

    //Coincidences
    if(data.nBa > 0 && data.nSe > 0) {
      
      for(int i=0;i<data.nBa;i++) {

	int bDet = data.bam2[i].det;
	int ring = data.bam2[i].ring;

	double ring_en = data.bam2[i].rEn;
        double sec_en = data.bam2[i].sEn;

	TVector3 bPos = GetPos(data.bam2[i].det,data.bam2[i].ring,data.bam2[i].sector);
	
	if(bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile DS gate

	  sPidDS->Fill(ring,sec_en);
	  rPidDS->Fill(ring,ring_en);

	  bool sol2 = false;
	  if(sec_en < Sol2_En) {
	    sol2 = true;
	  }
	  
	  double thetaCM = Theta_CM_FP(bPos.Theta(),beam_en,sol2);
	  double thetaCM_ns2 = Theta_CM_FP(bPos.Theta(),beam_en,false);
	  
	  double energy = KE_LAB(thetaCM,beam_en);
	  double gam = (energy - 12.0)/beam_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  double energy_ns2 = KE_LAB(thetaCM_ns2,beam_en);
	  double gam_ns2 = (energy_ns2 - 12.0)/beam_mass + 1.0;
	  double beta_ns2 = TMath::Sqrt(1.0 - 1.0/(gam_ns2*gam_ns2));

	  double recon_energy = Recoil_KE_LAB(thetaCM,beam_en);
	  double recon_gam = (recon_energy - 12.0)/targ_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));

	  double recon_energy_ns2 = Recoil_KE_LAB(thetaCM_ns2,beam_en);
	  double recon_gam_ns2 = (recon_energy_ns2 - 12.0)/targ_mass + 1.0;
	  double recon_beta_ns2 = TMath::Sqrt(1.0 - 1.0/(recon_gam_ns2*recon_gam_ns2));
	  
	  TVector3 rPos(0,0,1);
	  rPos.SetTheta(Recoil_Theta_LAB(thetaCM,beam_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  TVector3 rPos_ns2(0,0,1);
	  rPos_ns2.SetTheta(Recoil_Theta_LAB(thetaCM_ns2,beam_en));
	  rPos_ns2.SetPhi(bPos.Phi() - TMath::Pi());
	  
	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double coreEn = data.sega[i].cEn;
	  
	    TVector3 sPos = GetPos(det,seg);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	    double dopEn_ns2 = gam_ns2*(1 - beta_ns2*TMath::Cos(theta))*coreEn;
	    
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    double recon_theta_ns2 = rPos_ns2.Angle(sPos);
	    double recon_en_ns2 = recon_gam_ns2*(1 - recon_beta_ns2*TMath::Cos(recon_theta_ns2))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0) {
	      recon_phi += TMath::TwoPi();
	    }

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0) {
	      reconAng += TMath::TwoPi();
	    }
	    
	    pCoreEnergyDS->Fill(coreEn);
	    pCoreSumDS->Fill(det,coreEn);
	    
	    pDopEnergyDS->Fill(dopEn);
	    pDopSumDS->Fill(det,dopEn);

	    pDopvPartDS->Fill(dopEn,sec_en);
	    pDopvPartNS2DS->Fill(dopEn_ns2,sec_en);
	    
	    pThCorDS->Fill(coreEn,theta*TMath::RadToDeg());
	    pThCrtDS->Fill(dopEn,theta*TMath::RadToDeg());

	    pPhCorDS->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pPhCrtDS->Fill(dopEn,planeAng*TMath::RadToDeg());  

	    /*
	    pDetDopEnDS.at(det-1)->Fill(dopEn);
	    pDetDopSumDS.at(det-1)->Fill(seg,dopEn);

	    pDetThCorDS.at(det-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    pDetThCrtDS.at(det-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    pDetPhCorDS.at(det-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pDetPhCrtDS.at(det-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */

	    pRingCoreEnDS.at(ring-1)->Fill(coreEn);
	    pRingDopEnDS.at(ring-1)->Fill(dopEn);

	    /*
	    pRingThCorDS.at(ring-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    pRingThCrtDS.at(ring-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    pRingPhCorDS.at(ring-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pRingPhCrtDS.at(ring-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */
	    

	    pReconEnergyDS->Fill(recon_en);
	    pReconSumDS->Fill(det,recon_en);

	    pReconvPartDS->Fill(recon_en,sec_en);
	    pReconvPartNS2DS->Fill(recon_en_ns2,sec_en);

	    pReconThCorDS->Fill(coreEn,recon_theta*TMath::RadToDeg());
	    pReconThCrtDS->Fill(recon_en,recon_theta*TMath::RadToDeg());

	    pReconPhCorDS->Fill(coreEn,reconAng*TMath::RadToDeg());
	    pReconPhCrtDS->Fill(recon_en,reconAng*TMath::RadToDeg());

	    pRingRecEnDS.at(ring-1)->Fill(recon_en);
	    
	  } //End SeGA loop
	} //End DS projectile gate

	else if(!bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile US gate

	  double energy = KE_LAB(Theta_CM_FP(bPos.Theta(),beam_en),beam_en);
	  double gam = (energy - 12.0)/beam_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  sPidUS->Fill(ring,sec_en);
	  rPidUS->Fill(ring,ring_en);
	  
	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double coreEn = data.sega[i].cEn;
	  
	    TVector3 sPos = GetPos(det,seg);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    pCoreEnergyUS->Fill(coreEn);
	    pCoreSumUS->Fill(det,coreEn);
	    
	    pDopEnergyUS->Fill(dopEn);
	    pDopSumUS->Fill(det,dopEn);

	    pDopvPartUS->Fill(dopEn,sec_en);

	    pThCorUS->Fill(coreEn,theta*TMath::RadToDeg());
	    pThCrtUS->Fill(dopEn,theta*TMath::RadToDeg());

	    pPhCorUS->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pPhCrtUS->Fill(dopEn,planeAng*TMath::RadToDeg());  

	    /*
	    pDetDopEnUS.at(det-1)->Fill(dopEn);
	    pDetDopSumUS.at(det-1)->Fill(seg,dopEn);

	    pDetThCorUS.at(det-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    pDetThCrtUS.at(det-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    pDetPhCorUS.at(det-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pDetPhCrtUS.at(det-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */

	    pRingCoreEnUS.at(ring-1)->Fill(coreEn);
	    pRingDopEnUS.at(ring-1)->Fill(dopEn);

	    /*
	    pRingThCorUS.at(ring-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    pRingThCrtUS.at(ring-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    pRingPhCorUS.at(ring-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    pRingPhCrtUS.at(ring-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */
	    
	  
	  } //End SeGA loop
	} ////End US projectile gate

	else if(data.bam2[i].rR && data.bam2[i].sR) { //Recoil gate

	  sPidRec->Fill(ring,sec_en);
	  rPidRec->Fill(ring,ring_en);

	  double thetaCM = Theta_CM_FR(bPos.Theta(),beam_en);
	  
	  double energy = Recoil_KE_LAB(thetaCM,beam_en);
	  double gam = (energy - 12.0)/targ_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  //double recon_theta = Theta_LAB(thetaCM,beam_en);
	  double recon_energy = KE_LAB(thetaCM,beam_en);
	  double recon_gam = (recon_energy - 12.0)/beam_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));

	  TVector3 rPos(0,0,1);
	  rPos.SetTheta(Theta_LAB(thetaCM,beam_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double coreEn = data.sega[i].cEn;
	  
	    TVector3 sPos = GetPos(det,seg);

	    double theta = bPos.Angle(sPos);
	    double dopEn = gam*(1 - beta*TMath::Cos(theta))*coreEn;
	  
	    TVector3 reacPlane = bPos.Cross(incBeam);
	    TVector3 detPlane = sPos.Cross(incBeam);

	    double reac_phi = reacPlane.Phi();
	    if(reac_phi < 0) {
	      reac_phi += TMath::TwoPi();
	    }

	    double det_phi = detPlane.Phi();
	    if(det_phi < 0) {
	      det_phi += TMath::TwoPi();
	    }

	    double planeAng = reac_phi - det_phi;
	    if(planeAng < 0) {
	      planeAng += TMath::TwoPi();
	    }

	    double recon_theta = rPos.Angle(sPos);
	    double recon_en = recon_gam*(1 - recon_beta*TMath::Cos(recon_theta))*coreEn;

	    TVector3 reconPlane = rPos.Cross(incBeam);

	    double recon_phi = reconPlane.Phi();
	    if(recon_phi < 0) {
	      recon_phi += TMath::TwoPi();
	    }

	    double reconAng = recon_phi - det_phi;
	    if(reconAng < 0) {
	      reconAng += TMath::TwoPi();
	    }
	    
	    rCoreEnergy->Fill(coreEn);
	    rCoreSum->Fill(det,coreEn);
	    
	    rDopEnergy->Fill(dopEn);
	    rDopSum->Fill(det,dopEn);

	    rDopvPart->Fill(dopEn,sec_en);

	    rThCor->Fill(coreEn,theta*TMath::RadToDeg());
	    rThCrt->Fill(dopEn,theta*TMath::RadToDeg());

	    rPhCor->Fill(coreEn,planeAng*TMath::RadToDeg());
	    rPhCrt->Fill(dopEn,planeAng*TMath::RadToDeg());

	    /*
	    rDetDopEn.at(det-1)->Fill(dopEn);
	    rDetDopSum.at(det-1)->Fill(seg,dopEn);

	    rDetThCor.at(det-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    rDetThCrt.at(det-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    rDetPhCor.at(det-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    rDetPhCrt.at(det-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */

	    rRingCoreEn.at(ring-1)->Fill(coreEn);
	    rRingDopEn.at(ring-1)->Fill(dopEn);

	    /*
	    rRingThCor.at(ring-1)->Fill(coreEn,theta*TMath::RadToDeg());
	    rRingThCrt.at(ring-1)->Fill(dopEn,theta*TMath::RadToDeg());

	    rRingPhCor.at(ring-1)->Fill(coreEn,planeAng*TMath::RadToDeg());
	    rRingPhCrt.at(ring-1)->Fill(dopEn,planeAng*TMath::RadToDeg());
	    */

	    rReconEnergy->Fill(recon_en);
	    rReconSum->Fill(det,recon_en);

	    rReconvPart->Fill(recon_en,sec_en);

	    rReconThCor->Fill(coreEn,recon_theta*TMath::RadToDeg());
	    rReconThCrt->Fill(recon_en,recon_theta*TMath::RadToDeg());

	    rReconPhCor->Fill(coreEn,reconAng*TMath::RadToDeg());
	    rReconPhCrt->Fill(recon_en,reconAng*TMath::RadToDeg());

	    rRingRecEn.at(ring-1)->Fill(recon_en);

	    rRingRecThCor.at(ring-1)->Fill(coreEn,recon_theta*TMath::RadToDeg());
	    rRingRecThCrt.at(ring-1)->Fill(recon_en,recon_theta*TMath::RadToDeg());
	    
	  } //End SeGA loop  
	} //End recoil gate
	
      } //End Bambino2 loop
    } //End coincidences

    
  } //End while loop
  fclose(input_file);

  std::cout << "Writing histograms to file..." << std::endl;

  TFile* outFile = new TFile(output_filename,"RECREATE");
  outFile->mkdir("SeGA");
  outFile->mkdir("Bambino2");
  
  outFile->mkdir("Coincidence/ProjectileDS");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler");
  outFile->mkdir("Coincidence/ProjectileDS/Doppler/Rings");
  outFile->mkdir("Coincidence/ProjectileDS/Recon");
  outFile->mkdir("Coincidence/ProjectileDS/Recon/Rings");
  /*
  outFile->mkdir("Coincidence/ProjectileDS/Bambino2Rings");
  outFile->mkdir("Coincidence/ProjectileDS/Bambino2Rings/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileDS/Bambino2Rings/ThetaCorr");
  
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileDS/SegaDets/Summaries");
  */

  outFile->mkdir("Coincidence/ProjectileUS");
  outFile->mkdir("Coincidence/ProjectileUS/Rings");
  /*
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings");
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/Bambino2Rings/ThetaCorr");
  
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/ProjectileUS/SegaDets/Summaries");
  */
  
  outFile->mkdir("Coincidence/Recoil");
  outFile->mkdir("Coincidence/Recoil/Doppler");
  outFile->mkdir("Coincidence/Recoil/Doppler/Rings");
  outFile->mkdir("Coincidence/Recoil/Recon");
  outFile->mkdir("Coincidence/Recoil/Recon/Rings");
  /*
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings");
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/Bambino2Rings/ThetaCorr");

  outFile->mkdir("Coincidence/Recoil/SegaDets");
  outFile->mkdir("Coincidence/Recoil/SegaDets/PhiCorr");
  outFile->mkdir("Coincidence/Recoil/SegaDets/ThetaCorr");
  outFile->mkdir("Coincidence/Recoil/SegaDets/Summaries");
  */

  outFile->cd("Bambino2");

  bSum->Write();
  
  rPid0->Write();
  sPid0->Write();
  secD0->Write();
  
  rPid1->Write();
  sPid1->Write();
  secD1->Write();

  sPid1_p->Write();
  sPid1_r->Write();

  rPid1m1->Write();
  sPid1m1->Write();
  rPid1m2->Write();
  sPid1m2->Write();
  
  outFile->cd("SeGA");

  coreEnergy->Write();
  coreSum->Write();
  segEnergy->Write();
  segSum->Write();

  outFile->cd("Coincidence/ProjectileDS");

  sPidDS->Write();
  rPidDS->Write();

  pCoreEnergyDS->Write();
  pCoreSumDS->Write();

  outFile->cd("Coincidence/ProjectileDS/Doppler");
  
  pDopEnergyDS->Write();
  pDopSumDS->Write();

  pDopvPartDS->Write();
  pDopvPartNS2DS->Write();

  pThCorDS->Write();
  pThCrtDS->Write();

  pPhCorDS->Write();
  pPhCrtDS->Write();

  outFile->cd("Coincidence/ProjectileDS/Doppler/Rings");

  for(int i=0;i<24;i++) {
    pRingCoreEnDS.at(i)->Write();
    pRingDopEnDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/Recon");

  pReconEnergyDS->Write();
  pReconSumDS->Write();

  pReconvPartDS->Write();
  pReconvPartNS2DS->Write();

  pReconThCorDS->Write();
  pReconThCrtDS->Write();

  pReconPhCorDS->Write();
  pReconPhCrtDS->Write();

  outFile->cd("Coincidence/ProjectileDS/Recon/Rings");

  for(int i=0;i<24;i++) {
    pRingRecEnDS.at(i)->Write();
  }

  /*

  outFile->cd("Coincidence/ProjectileDS/Bambino2Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    pRingPhCorDS.at(i)->Write();
    pRingPhCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/Bambino2Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    pRingThCorDS.at(i)->Write();
    pRingThCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets");
  for(int i=0;i<16;i++) {
    pDetDopEnDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    pDetDopSumDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    pDetPhCorDS.at(i)->Write();
    pDetPhCrtDS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileDS/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    pDetThCorDS.at(i)->Write();
    pDetThCrtDS.at(i)->Write();
  }
  */

  outFile->cd("Coincidence/ProjectileUS");

  sPidUS->Write();
  rPidUS->Write();

  pCoreEnergyUS->Write();
  pCoreSumUS->Write();
  
  pDopEnergyUS->Write();
  pDopSumUS->Write();

  pDopvPartUS->Write();

  pThCorUS->Write();
  pThCrtUS->Write();

  pPhCorUS->Write();
  pPhCrtUS->Write();

  outFile->cd("Coincidence/ProjectileUS/Rings");

  for(int i=0;i<24;i++) {
    pRingCoreEnUS.at(i)->Write();
    pRingDopEnUS.at(i)->Write();
  }

  /*
  outFile->cd("Coincidence/ProjectileUS/Bambino2Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    pRingPhCorUS.at(i)->Write();
    pRingPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/Bambino2Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    pRingThCorUS.at(i)->Write();
    pRingThCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets");
  for(int i=0;i<16;i++) {
    pDetDopEnUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    pDetDopSumUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    pDetPhCorUS.at(i)->Write();
    pDetPhCrtUS.at(i)->Write();
  }

  outFile->cd("Coincidence/ProjectileUS/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    pDetThCorUS.at(i)->Write();
    pDetThCrtUS.at(i)->Write();
  }
  */
  
  outFile->cd("Coincidence/Recoil");

  sPidRec->Write();
  rPidRec->Write();

  rCoreEnergy->Write();
  rCoreSum->Write();

  outFile->cd("Coincidence/Recoil/Doppler");
  
  rDopEnergy->Write();
  rDopSum->Write();

  rDopvPart->Write();

  rThCor->Write();
  rThCrt->Write();

  rPhCor->Write();
  rPhCrt->Write();

  outFile->cd("Coincidence/Recoil/Doppler/Rings");

  for(int i=0;i<24;i++) {
    rRingCoreEn.at(i)->Write();
    rRingDopEn.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/Recon");

  rReconEnergy->Write();
  rReconSum->Write();

  rReconvPart->Write();

  rReconThCor->Write();
  rReconThCrt->Write();

  rReconPhCor->Write();
  rReconPhCrt->Write();

  outFile->cd("Coincidence/Recoil/Recon/Rings");
  for(int i=0;i<24;i++) {
    rRingRecEn.at(i)->Write();

    rRingRecThCor.at(i)->Write();
    rRingRecThCrt.at(i)->Write();
  }

  /* 
  outFile->cd("Coincidence/Recoil/Bambino2Rings/PhiCorr");
  for(int i=0;i<24;i++) {
    rRingPhCor.at(i)->Write();
    rRingPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/Bambino2Rings/ThetaCorr");
  for(int i=0;i<24;i++) {
    rRingThCor.at(i)->Write();
    rRingThCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets");
  for(int i=0;i<16;i++) {
    rDetDopEn.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/Summaries");
  for(int i=0;i<16;i++) {
    rDetDopSum.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/PhiCorr");
  for(int i=0;i<16;i++) {
    rDetPhCor.at(i)->Write();
    rDetPhCrt.at(i)->Write();
  }

  outFile->cd("Coincidence/Recoil/SegaDets/ThetaCorr");
  for(int i=0;i<16;i++) {
    rDetThCor.at(i)->Write();
    rDetThCrt.at(i)->Write();
  }
  */

  outFile->Close();

  std::cout << "Done!" << std::endl;

  return 0;
}
