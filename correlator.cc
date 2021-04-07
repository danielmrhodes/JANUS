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

////These should match the parameters defined in the simulation input////
//These correspond to the input file Examples/Macros/full.mac

const int beamZ = 48;
const double beam_mass = 98626.9; // MeV/c^2

//You should reduce this value by the energy loss in the target
const double beam_en = 265.0; // MeV
  
const int targZ = 22;
const double targ_mass = 44652.0; // MeV/c^2

//Silicon detector z-offsets (downstream and upstream)
double DS_Offset = 2.6; // cm 
double US_Offset = 3.4; // cm

//Beam spot position
const double beam_X = 0.0; // cm 
const double beam_Y = 0.0; // cm 

double SeGA_Offset = 3.1; // cm
/////////////////////////////////////////////////////////////////////////

////SeGA Resolution////
double Sigma(double en) {
  return 1.03753 + en*0.000274797;
}
///////////////////////

////Kinematics////
double Theta_CM_FP(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.) { //From projectile

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
 
  return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + TMath::Pi();

}

double Theta_CM_FR(double ThetaLAB, double Ep, bool sol2=false, double Ex=0.) { //From recoil

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
  
  return -std::asin(tau*std::sin(-ThetaLAB)) - ThetaLAB;

}

double Theta_LAB_Max(double Ep, double Ex=0.0) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  if(tau < 1.0) {
    return TMath::Pi();
  }
  
  return std::asin(1.0/tau);
  
}

double Theta_LAB(double thetaCM, double Ep, double Ex=0.) {

  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(thetaCM)/(std::cos(thetaCM) + tau);

  if(tanTheta > 0) {
    return std::atan(tanTheta);
  }
  
  return std::atan(tanTheta) + TMath::Pi();
  
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
    zoff = DS_Offset;
  }
  else {
    zoff = US_Offset;
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

  TVector3 pos(1.0,1.0,1.0);
  pos.SetPerp((outerRadius + innerRadius)/2.0);
  pos.SetPhi((quad+0.5)*2*PI/4.0);
  pos.SetZ((length/8.0)*(2.0*slice - 7.0));

  double rd = 12.975;
  double phid = (det-1)*(2.0*PI/8.0) + PI/8.0 ;
  double zd = length + 2*0.05 + 0.6;
  if(det > 8) {
    zd*=-1;
  }
  
  TVector3 origin(rd*TMath::Cos(phid),rd*TMath::Sin(phid),zd+SeGA_Offset);

  return origin+pos;

}
/////////////////

////////////Build data format////////////
struct BAM2 {

  //realistic info
  int det, ring, sector;
  double rEn, sEn;

  //perfect info
  TVector3 rPos, sPos;
  bool rP, rR, sP, sR;
  
};

struct SEGA {

  //realistic info
  int det, nsegs, segs[32];
  double cEn, sEn[32];
  
  //perfect info
  double x[32], y[32], z[32];
  bool fep, pfep;

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

    bam2[nBa].det = ringHit.det; //either hit works
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
/////////////////////////////////////////

//Unpack raw data into correlated data
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
    double x = dat.sData[i].x;
    double y = dat.sData[i].y;
    double z = dat.sData[i].z;
    
    bool FEP = dat.sData[i].fep;
    bool PFEP = dat.sData[i].pfep;
    
    if(!exists.at(detect-1)) {
      data.sega[nS].det = detect;

      if((bool)segment) {
        data.sega[nS].nsegs = 1;
	data.sega[nS].segs[0] = segment;
	data.sega[nS].sEn[0] = energy;

	data.sega[nS].x[0] = x;
        data.sega[nS].y[0] = y;
        data.sega[nS].z[0] = z;
      }
      else {
	data.sega[nS].nsegs = 0;
	data.sega[nS].cEn = energy;
	data.sega[nS].fep = FEP;
	data.sega[nS].pfep = PFEP;
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

	data.sega[index].x[Nsegs] = x;
        data.sega[index].y[Nsegs] = y;
        data.sega[index].z[Nsegs] = z;
	
	data.sega[index].nsegs++;
      }
      else {
	data.sega[index].cEn = energy;
	data.sega[index].fep = FEP;
	data.sega[index].pfep = PFEP;
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
  if(!strcmp(input_filename,output_filename)) {
    std::cout << "Give your input and output files different names" << std::endl;
    return 1;
  }
  
  //Bambino2 singles
  TH2* bSum = new TH2D("Summary","Janus Summary",120,1,121,500,0,500);
  TH2* pSum = new TH2D("pSummary","Janus Projectile Summary",120,1,121,500,0,500);
  TH2* rSum = new TH2D("rSummary","Janus Recoil Summary",120,1,121,500,0,500);

  TH2* pThvPhDS = new TH2D("pThvPhDS","Projectile #phi-#theta surface",1000,0,90,1000,-200,200);
  TH2* pThvPhDS1 = new TH2D("pThvPhDS1","Projectile #phi-#theta surface 1",1000,0,90,1000,-200,200);
  TH2* pThvPhDS2 = new TH2D("pThvPhDS2","Projectile #phi-#theta surface 2",1000,0,90,1000,-200,200);

  TH2* rThvPh = new TH2D("rThvPh","Recoil #phi-#theta surface",1000,0,90,1000,-200,200);
  
  TH1* secD0 = new TH1D("Sectors_Det0","US Sectors",32,1,33);
  TH1* secD1 = new TH1D("Sectors_Det1","DS Sectors",32,1,33);

  TH1* pSecDS = new TH1D("pSecDS","DS Projectile Sectors",32,1,33);
  TH1* rSec = new TH1D("rSec","Recoil Sectors",32,1,33);
  
  TH2* rPid0 = new TH2D("RingPID_Det0","US RingEn PID",24,1,25,500,0,500);
  TH2* rPid1 = new TH2D("RingPID_Det1","DS RingEn PID",24,1,25,500,0,500);

  TH2* sPid0 = new TH2D("SecPID_Det0","US SectorEn PID",24,1,25,500,0,500);
  TH2* sPid1 = new TH2D("SecPID_Det1","DS SectorEn PID",24,1,25,500,0,500);

  TH2* sPid1_p = new TH2D("SecPID_Det1_proj","DS SectorEn Projectile PID",24,1,25,500,0,500);
  TH2* sPid1_r = new TH2D("SecPID_Det1_rec","DS SectorEn Recoil PID",24,1,25,500,0,500);
  
  //SeGA singles
  TH1* coreEnergy = new TH1D("Core_Energy","SeGA Core Energy",3000,0,3000);
  TH1* segEnergy = new TH1D("Seg_Energy","SeGA Segment Energy",3000,0,3000);

  TH2* coreSum = new TH2D("Core_Summary","Core Energy Summary",16,1,17,3000,0,3000);
  TH2* segSum = new TH2D("Seg_Summary","Segment Energy Summary",512,1,513,3000,0,3000);

  TH1* coreEn_Fep = new TH1D("FEP","SeGA FEP",3000,0,3000);
  TH1* coreEn_NotFep = new TH1D("nFEP","SeGA Not FEP",3000,0,3000);

  TH2* sPosTP = new TH2D("sPosPT","SeGA Phi-Theta Surface",360,0,180,760,-10,370);

  //Coincidences
  //Projectile DS
  TH2* sPidDS = new TH2D("SecPID_DS","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidDS = new TH2D("RingPID_DS","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* pCoreEnergyDS = new TH1D("Core_EnergyDS","SeGA Core Energy",12000,0,4000);
  TH2* pCoreSumDS = new TH2D("Core_SummaryDS","Core Energy Summary",16,1,17,3000,0,3000);

  TH1* pDopEnergyDS = new TH1D("Dop_EnergyDS","Doppler Energy",12000,0,4000);
  TH2* pDopSumDS = new TH2D("Dop_SummaryDS","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH1* pCoreEnergyDS_fep = new TH1D("Core_EnergyDS_fep","SeGA Core Energy FEP",12000,0,4000);
  TH1* pCoreEnergyDS_pfep = new TH1D("Core_EnergyDS_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  TH1* pCoreEnergyDS_rfep = new TH1D("Core_EnergyDS_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);
  TH1* pCoreEnergyDS_nfep = new TH1D("Core_EnergyDS_nfep","SeGA Core Energy Not FEP",12000,0,4000);

  TH1* pDopEnergyDS_fep = new TH1D("Dop_EnergyDS_fep","Doppler Energy FEP",12000,0,4000);
  TH1* pDopEnergyDS_pfep = new TH1D("Dop_EnergyDS_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  TH1* pDopEnergyDS_rfep = new TH1D("Dop_EnergyDS_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  TH1* pDopEnergyDS_nfep = new TH1D("Dop_EnergyDS_nfep","Doppler Energy Not FEP",12000,0,4000);

  TH2* pDopvPartDS = new TH2D("DopEn_v_PartEn_DS","Doppler Energy vs Particle Energy",
			      3000,0,3000,500,0,500);

  TH2* pThCorDS = new TH2D("Theta_CorrDS","Theta Correlation",3000,0,3000,90,0,180);
  TH2* pThCrtDS = new TH2D("Theta_CrctDS","Theta Correction",6000,0,3000,90,0,180);

  double thing1 = 65.*180./32.0;
  TH2* pPhCorDS = new TH2D("Phi_CorrDS","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pPhCrtDS = new TH2D("Phi_CrctDS","Phi Correction",6000,0,3000,32,0,thing1);

  TH1* pReconEnergyDS = new TH1D("Recon_EnergyDS","Recon Energy",12000,0,4000);
  TH2* pReconSumDS = new TH2D("Recon_SummaryDS","Recon Energy Summary",16,1,17,6000,0,3000);

  TH1* pReconEnergyDS_fep = new TH1D("Recon_EnergyDS_fep","Recon Energy FEP",12000,0,4000);
  TH1* pReconEnergyDS_pfep = new TH1D("Recon_EnergyDS_pfep","Recon Energy Projectile FEP",12000,0,4000);
  TH1* pReconEnergyDS_rfep = new TH1D("Recon_EnergyDS_rfep","Recon Energy Recoil FEP",12000,0,4000);
  TH1* pReconEnergyDS_nfep = new TH1D("Recon_EnergyDS_nfep","Recon Energy Not FEP",12000,0,4000);

  TH2* pReconvPartDS = new TH2D("ReconEn_v_partEn_DS","Recon Energy vs Particle Energy",
				3000,0,3000,500,0,500);

  TH2* pReconThCorDS = new TH2D("ReconTheta_CorrDS","Recon Theta Correlation",3000,0,3000,90,0,180);
  TH2* pReconThCrtDS = new TH2D("ReconTheta_CrctDS","Recon Theta Correction",6000,0,3000,90,0,180);

  TH2* pReconPhCorDS = new TH2D("ReconPhi_CorrDS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pReconPhCrtDS = new TH2D("ReconPhi_CrctDS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Projectile US
  TH2* sPidUS = new TH2D("SecPID_US","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidUS = new TH2D("RingPID_US","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* pCoreEnergyUS = new TH1D("Core_EnergyUS","SeGA Core Energy",12000,0,4000);
  TH2* pCoreSumUS = new TH2D("Core_SummaryUS","Core Energy Summary",16,1,17,3000,0,3000);
  
  TH1* pDopEnergyUS = new TH1D("Dop_EnergyUS","Doppler Energy",12000,0,4000);
  TH2* pDopSumUS = new TH2D("Dop_SummaryUS","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH1* pCoreEnergyUS_fep = new TH1D("Core_EnergyUS_fep","SeGA Core Energy FEP",12000,0,4000);
  TH1* pCoreEnergyUS_pfep = new TH1D("Core_EnergyUS_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  TH1* pCoreEnergyUS_rfep = new TH1D("Core_EnergyUS_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);
  TH1* pCoreEnergyUS_nfep = new TH1D("Core_EnergyUS_nfep","SeGA Core Energy Not FEP",12000,0,4000);

  TH1* pDopEnergyUS_fep = new TH1D("Dop_EnergyUS_fep","Doppler Energy FEP",12000,0,4000);
  TH1* pDopEnergyUS_pfep = new TH1D("Dop_EnergyUS_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  TH1* pDopEnergyUS_rfep = new TH1D("Dop_EnergyUS_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  TH1* pDopEnergyUS_nfep = new TH1D("Dop_EnergyUS_nfep","Doppler Energy Not FEP",12000,0,4000);

  TH2* pDopvPartUS = new TH2D("Dop_PartEn_US","Doppler Energy vs Particle Energy",3000,0,3000,500,0,500);

  TH2* pThCorUS = new TH2D("Theta_CorrUS","Theta Correlation",3000,0,3000,90,0,180);
  TH2* pThCrtUS = new TH2D("Theta_CrctUS","Theta Correction",6000,0,3000,90,0,180);
  
  TH2* pPhCorUS = new TH2D("Phi_CorrUS","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pPhCrtUS = new TH2D("Phi_CrctUS","Phi Correction",6000,0,3000,32,0,thing1);

  TH1* pReconEnergyUS = new TH1D("Recon_EnergyUS","Recon Energy",12000,0,4000);
  TH2* pReconSumUS = new TH2D("Recon_SummaryUS","Recon Energy Summary",16,1,17,6000,0,3000);

  TH1* pReconEnergyUS_fep = new TH1D("Recon_EnergyUS_fep","Recon Energy FEP",12000,0,4000);
  TH1* pReconEnergyUS_pfep = new TH1D("Recon_EnergyUS_pfep","Recon Energy Projectile FEP",12000,0,4000);
  TH1* pReconEnergyUS_rfep = new TH1D("Recon_EnergyUS_rfep","Recon Energy Recoil FEP",12000,0,4000);
  TH1* pReconEnergyUS_nfep = new TH1D("Recon_EnergyUS_nfep","Recon Energy Not FEP",12000,0,4000);

  TH2* pReconvPartUS = new TH2D("ReconEn_v_partEn_US","Recon Energy vs Particle Energy",
				3000,0,3000,500,0,500);

  TH2* pReconThCorUS = new TH2D("ReconTheta_CorrUS","Recon Theta Correlation",3000,0,3000,90,0,180);
  TH2* pReconThCrtUS = new TH2D("ReconTheta_CrctUS","Recon Theta Correction",6000,0,3000,90,0,180);

  TH2* pReconPhCorUS = new TH2D("ReconPhi_CorrUS","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* pReconPhCrtUS = new TH2D("ReconPhi_CrctUS","Recon Phi Correction",6000,0,3000,32,0,thing1);

  //Recoil
  TH2* sPidRec = new TH2D("SecPID_Rec","DS SectorEn PID",24,1,25,500,0,500);
  TH2* rPidRec = new TH2D("RingPID_Rec","DS RingEn PID",24,1,25,500,0,500);
  
  TH1* rCoreEnergy = new TH1D("Core_EnergyRec","SeGA Core Energy",12000,0,4000);
  TH2* rCoreSum = new TH2D("Core_SummaryRec","Core Energy Summary",16,1,17,3000,0,3000);
  
  TH1* rDopEnergy = new TH1D("Dop_EnergyRec","Doppler Energy",12000,0,4000);
  TH2* rDopSum = new TH2D("Dop_SummaryRec","Doppler Energy Summary",16,1,17,6000,0,3000);

  TH1* rCoreEnergy_fep = new TH1D("Core_EnergyRec_fep","SeGA Core Energy FEP",12000,0,4000);
  TH1* rCoreEnergy_pfep = new TH1D("Core_EnergyRec_pfep","SeGA Core Energy Projectile FEP",12000,0,4000);
  TH1* rCoreEnergy_rfep = new TH1D("Core_EnergyRec_rfep","SeGA Core Energy Recoil FEP",12000,0,4000);
  TH1* rCoreEnergy_nfep = new TH1D("Core_EnergRec_nfep","SeGA Core Energy Not FEP",12000,0,4000);

  TH1* rDopEnergy_fep = new TH1D("Dop_EnergyRec_fep","Doppler Energy FEP",12000,0,4000);
  TH1* rDopEnergy_pfep = new TH1D("Dop_EnergyRec_pfep","Doppler Energy Projectile FEP",12000,0,4000);
  TH1* rDopEnergy_rfep = new TH1D("Dop_EnergyRec_rfep","Doppler Energy Recoil FEP",12000,0,4000);
  TH1* rDopEnergy_nfep = new TH1D("Dop_EnergyRec_nfep","Doppler Energy Not FEP",12000,0,4000);

  TH2* rDopvPart = new TH2D("DopEv_v_PartEn_Rec","Doppler Energy vs Particle Energy",3000,0,3000,500,0,500);

  TH2* rThCor = new TH2D("Theta_CorrRec","Theta Correlation",3000,0,3000,90,0,180);
  TH2* rThCrt = new TH2D("Theta_CrctRec","Theta Correction",6000,0,3000,90,0,180);;
  TH2* rPhCor = new TH2D("Phi_CorrRec","Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* rPhCrt = new TH2D("Phi_CrctRec","Phi Correction",6000,0,3000,32,0,thing1);

  TH1* rReconEnergy = new TH1D("Recon_EnergyRec","Recon Energy",12000,0,4000);
  TH2* rReconSum = new TH2D("Recon_SummaryRec","Recon Energy Summary",16,1,17,6000,0,3000);

  TH1* rReconEnergy_fep = new TH1D("Recon_EnergyRec_fep","Recon Energy FEP",12000,0,4000);
  TH1* rReconEnergy_pfep = new TH1D("Recon_EnergyRec_pfep","Recon Energy Projectile FEP",12000,0,4000);
  TH1* rReconEnergy_rfep = new TH1D("Recon_EnergyRec_rfep","Recon Energy Recoil FEP",12000,0,4000);
  TH1* rReconEnergy_nfep = new TH1D("Recon_EnergyRec_nfep","Recon Energy Not FEP",12000,0,4000);

  TH2* rReconvPart = new TH2D("ReconEn_v_PartEn_Rec","Recon Energy vs Particle Energy",3000,0,3000,500,0,500);

  TH2* rReconThCor = new TH2D("ReconTheta_CorrRec","Recon Theta Correlation",3000,0,3000,90,0,180);
  TH2* rReconThCrt = new TH2D("ReconTheta_CrctRec","Recon Theta Correction",6000,0,3000,90,0,180);

  TH2* rReconPhCor = new TH2D("ReconPhi_CorrRec","Recon Phi Correlation",3000,0,3000,32,0,thing1);
  TH2* rReconPhCrt = new TH2D("ReconPhi_CrctRec","Recon Phi Correction",6000,0,3000,32,0,thing1);

  std::cout << "Correlating and histograming data..." << std::endl;
  FILE* input_file = fopen(input_filename,"rb");
  
  const TVector3 incBeam = TVector3(0.0,0.0,1.0);
  const double Sol2_En = KE_LAB(Theta_CM_FP(Theta_LAB_Max(beam_en),beam_en),beam_en);
  const double r2d = TMath::RadToDeg();

  TRandom* rand = new TRandom(50747227);
   
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
      TVector3 pos = data.bam2[i].sPos;

      segPos.SetX(segPos.X() - beam_X);
      segPos.SetY(segPos.Y() - beam_Y);
      pos.SetX(pos.X() - beam_X);
      pos.SetY(pos.Y() - beam_Y);
      
      if(!det) { //Upstream

	bSum->Fill(sec,sec_en);
	bSum->Fill(ring+32,ring_en);
	
	rPid0->Fill(ring,ring_en);
	sPid0->Fill(ring,sec_en);
	secD0->Fill(sec);

	if(data.bam2[i].rP && data.bam2[i].sP) {
	  pSum->Fill(sec,sec_en);
	  pSum->Fill(ring+32,ring_en);
	}
	if(data.bam2[i].rR && data.bam2[i].sR) {
	  rSum->Fill(sec,sec_en);
	  rSum->Fill(ring+32,ring_en);
	}
      }
      else { //Downstream
	
	rPid1->Fill(ring,ring_en);
	sPid1->Fill(ring,sec_en);
	secD1->Fill(sec);

	bSum->Fill(sec+64,sec_en);
	bSum->Fill(ring+96,ring_en);

	if(data.bam2[i].rP && data.bam2[i].sP) { //Projectile
	  sPid1_p->Fill(ring,sec_en);

	  pSum->Fill(sec+64,sec_en);
	  pSum->Fill(ring+96,ring_en);

	  pSecDS->Fill(sec);
	  pThvPhDS->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  
	  if(ring%2) {
	    if(sec%2) {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }
	  else {
	    if(sec%2) {
	      pThvPhDS2->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	    else {
	      pThvPhDS1->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	    }
	  }  
 
	} //End projectile gate

	if(data.bam2[i].rR && data.bam2[i].sR) { //Recoil
	  sPid1_r->Fill(ring,sec_en);

	  rSum->Fill(sec+64,sec_en);
	  rSum->Fill(ring+96,ring_en);

	  rSec->Fill(sec);
	  rThvPh->Fill(pos.Theta()*r2d,pos.Phi()*r2d);
	  
	} //End recoil gate

      } //End downstream gate
      
    } //End Bambino2 singles

    //SeGA singles
    for(int i=0;i<data.nSe;i++) {

      int det = data.sega[i].det;
      double en = data.sega[i].cEn;
      double core_en = rand->Gaus(en,Sigma(en));
      
      coreEnergy->Fill(core_en);
      coreSum->Fill(det,core_en);

      for(int j=0;j<data.sega[i].nsegs;j++) {

	double seg_en = data.sega[i].sEn[j];
	int num = data.sega[i].segs[j] + 32*(det-1);

	double exactX = data.sega[i].x[j];
	double exactY = data.sega[i].y[j];
	double exactZ = data.sega[i].z[j];
	TVector3 exact_pos(exactX,exactY,exactZ);

	double theta = exact_pos.Theta()*r2d;
	double phi = exact_pos.Phi();
	if(phi < 0) {
	  phi += TMath::TwoPi();
	}
	phi *= r2d;	
	
	segEnergy->Fill(seg_en);
	segSum->Fill(num,seg_en);
	sPosTP->Fill(theta,phi);
	
      }

      if(data.sega[i].fep) {
	coreEn_Fep->Fill(core_en);
      }
      else {
	coreEn_NotFep->Fill(core_en);
      }
      
    } //End SeGA singles

    //Coincidences
    if(data.nBa > 0 && data.nSe > 0) {
      
      for(int i=0;i<data.nBa;i++) {

	int bDet = data.bam2[i].det;
	int ring = data.bam2[i].ring;
	int sector = data.bam2[i].sector;

	double ring_en = data.bam2[i].rEn;
        double sec_en = data.bam2[i].sEn;

	TVector3 bPos = GetPos(bDet,ring,sector);
	bPos.SetX(bPos.X() - beam_X);
	bPos.SetY(bPos.Y() - beam_Y);
	
	if(bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile DS gate

	  sPidDS->Fill(ring,sec_en);
	  rPidDS->Fill(ring,ring_en);

	  bool sol2 = false;
	  if(sec_en < Sol2_En) {
	    sol2 = true;
	  }
	  
	  double thetaCM = Theta_CM_FP(bPos.Theta(),beam_en,sol2);
	  
	  double energy = KE_LAB(thetaCM,beam_en);
	  double gam = (energy/beam_mass) + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  double recon_energy = Recoil_KE_LAB(thetaCM,beam_en);
	  double recon_gam = (recon_energy)/targ_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));
	  
	  TVector3 rPos(0,0,1); //Reconstructed position of the recoil
	  rPos.SetTheta(Recoil_Theta_LAB(thetaCM,beam_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());
	  
	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double en = data.sega[i].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[i].fep;
	    bool PFEP = data.sega[i].pfep;
	  
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
	    
	    pCoreEnergyDS->Fill(coreEn);
	    pCoreSumDS->Fill(det,coreEn);

	    if(FEP) {
	      pCoreEnergyDS_fep->Fill(coreEn);
	      pDopEnergyDS_fep->Fill(dopEn);
	      pReconEnergyDS_fep->Fill(recon_en);

	      if(PFEP) {
		pCoreEnergyDS_pfep->Fill(coreEn);
		pDopEnergyDS_pfep->Fill(dopEn);
		pReconEnergyDS_pfep->Fill(recon_en);
	      }
	      else {
		pCoreEnergyDS_rfep->Fill(coreEn);
		pDopEnergyDS_rfep->Fill(dopEn);
		pReconEnergyDS_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      pCoreEnergyDS_nfep->Fill(coreEn);
	      pDopEnergyDS_nfep->Fill(dopEn);
	      pReconEnergyDS_nfep->Fill(recon_en);
	    }
	    
	    pDopEnergyDS->Fill(dopEn);
	    pDopSumDS->Fill(det,dopEn);

	    pDopvPartDS->Fill(dopEn,sec_en);
	    
	    pThCorDS->Fill(coreEn,theta*r2d);
	    pThCrtDS->Fill(dopEn,theta*r2d);

	    pPhCorDS->Fill(coreEn,planeAng*r2d);
	    pPhCrtDS->Fill(dopEn,planeAng*r2d);

	    pReconEnergyDS->Fill(recon_en);
	    pReconSumDS->Fill(det,recon_en);

	    pReconvPartDS->Fill(recon_en,sec_en);

	    pReconThCorDS->Fill(coreEn,recon_theta*r2d);
	    pReconThCrtDS->Fill(recon_en,recon_theta*r2d);

	    pReconPhCorDS->Fill(coreEn,reconAng*r2d);
	    pReconPhCrtDS->Fill(recon_en,reconAng*r2d);
	    
	  } //End SeGA loop
	} //End DS projectile gate

	else if(!bDet && data.bam2[i].rP && data.bam2[i].sP) { //Projectile US gate

	  sPidUS->Fill(ring,sec_en);
	  rPidUS->Fill(ring,ring_en);
	  
	  double thetaCM = Theta_CM_FP(bPos.Theta(),beam_en,false);
	  
	  double energy = KE_LAB(thetaCM,beam_en);
	  double gam = (energy)/beam_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  double recon_energy = Recoil_KE_LAB(thetaCM,beam_en);
	  double recon_gam = (recon_energy)/targ_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));

	  TVector3 rPos(0,0,1); //Reconstructed position of the recoil
	  rPos.SetTheta(Recoil_Theta_LAB(thetaCM,beam_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());
	  
	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double en = data.sega[i].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[i].fep;
	    bool PFEP = data.sega[i].pfep;
	  
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

	    pCoreEnergyUS->Fill(coreEn);
	    pCoreSumUS->Fill(det,coreEn);
	    
	    pDopEnergyUS->Fill(dopEn);
	    pDopSumUS->Fill(det,dopEn);

	    if(FEP) {
	      pCoreEnergyUS_fep->Fill(coreEn);
	      pDopEnergyUS_fep->Fill(dopEn);
	      pReconEnergyUS_fep->Fill(recon_en);

	      if(PFEP) {
		pCoreEnergyUS_pfep->Fill(coreEn);
		pDopEnergyUS_pfep->Fill(dopEn);
		pReconEnergyUS_pfep->Fill(recon_en);
	      }
	      else {
		pCoreEnergyUS_rfep->Fill(coreEn);
		pDopEnergyUS_rfep->Fill(dopEn);
		pReconEnergyUS_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      pCoreEnergyUS_nfep->Fill(coreEn);
	      pDopEnergyUS_nfep->Fill(dopEn);
	      pReconEnergyUS_nfep->Fill(recon_en);
	    }
	    
	    pDopvPartUS->Fill(dopEn,sec_en);

	    pThCorUS->Fill(coreEn,theta*r2d);
	    pThCrtUS->Fill(dopEn,theta*r2d);

	    pPhCorUS->Fill(coreEn,planeAng*r2d);
	    pPhCrtUS->Fill(dopEn,planeAng*r2d);

	    pReconEnergyUS->Fill(recon_en);
	    pReconSumUS->Fill(det,recon_en);

	    pReconvPartUS->Fill(recon_en,sec_en);

	    pReconThCorUS->Fill(coreEn,recon_theta*r2d);
	    pReconThCrtUS->Fill(recon_en,recon_theta*r2d);

	    pReconPhCorUS->Fill(coreEn,reconAng*r2d);
	    pReconPhCrtUS->Fill(recon_en,reconAng*r2d);
	  
	  } //End SeGA loop
	} ////End US projectile gate

	else if(data.bam2[i].rR && data.bam2[i].sR) { //Recoil gate

	  sPidRec->Fill(ring,sec_en);
	  rPidRec->Fill(ring,ring_en);

	  double thetaCM = Theta_CM_FR(bPos.Theta(),beam_en);
	  
	  double energy = Recoil_KE_LAB(thetaCM,beam_en);
	  double gam = (energy)/targ_mass + 1.0;
	  double beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));

	  double recon_energy = KE_LAB(thetaCM,beam_en);
	  double recon_gam = (recon_energy)/beam_mass + 1.0;
	  double recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));

	  TVector3 rPos(0,0,1); //Reconstructed position of the projectile
	  rPos.SetTheta(Theta_LAB(thetaCM,beam_en));
	  rPos.SetPhi(bPos.Phi() - TMath::Pi());

	  for(int i=0;i<data.nSe;i++) {

	    int det = data.sega[i].det;
	    int seg = data.sega[i].MainSeg();
	    double en = data.sega[i].cEn;
	    double coreEn = rand->Gaus(en,Sigma(en));
	    bool FEP = data.sega[i].fep;
	    bool PFEP = data.sega[i].pfep;
	  
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

	    if(FEP) {
	      rCoreEnergy_fep->Fill(coreEn);
	      rDopEnergy_fep->Fill(dopEn);
	      rReconEnergy_fep->Fill(recon_en);

	      if(PFEP) {
		rCoreEnergy_pfep->Fill(coreEn);
		rDopEnergy_pfep->Fill(dopEn);
		rReconEnergy_pfep->Fill(recon_en);
	      }
	      else {
		rCoreEnergy_rfep->Fill(coreEn);
		rDopEnergy_rfep->Fill(dopEn);
		rReconEnergy_rfep->Fill(recon_en);
	      }
	    }
	    else {
	      rCoreEnergy_nfep->Fill(coreEn);
	      rDopEnergy_nfep->Fill(dopEn);
	      rReconEnergy_nfep->Fill(recon_en);
	    }

	    rDopvPart->Fill(dopEn,sec_en);

	    rThCor->Fill(coreEn,theta*r2d);
	    rThCrt->Fill(dopEn,theta*r2d);
	    
	    rPhCor->Fill(coreEn,planeAng*r2d);
	    rPhCrt->Fill(dopEn,planeAng*r2d);

	    rReconEnergy->Fill(recon_en);
	    rReconSum->Fill(det,recon_en);

	    rReconvPart->Fill(recon_en,sec_en);

	    rReconThCor->Fill(coreEn,recon_theta*r2d);
	    rReconThCrt->Fill(recon_en,recon_theta*r2d);

	    rReconPhCor->Fill(coreEn,reconAng*r2d);
	    rReconPhCrt->Fill(recon_en,reconAng*r2d);
	    
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
  outFile->mkdir("Coincidence/ProjectileDS/Recon");
  
  outFile->mkdir("Coincidence/ProjectileUS");
  outFile->mkdir("Coincidence/ProjectileUS/Doppler");
  outFile->mkdir("Coincidence/ProjectileUS/Recon");
  
  outFile->mkdir("Coincidence/Recoil");
  outFile->mkdir("Coincidence/Recoil/Doppler");
  outFile->mkdir("Coincidence/Recoil/Recon");

  outFile->cd("Bambino2");

  bSum->Write();
  pSum->Write();
  rSum->Write();
  
  pThvPhDS->Write();
  pThvPhDS1->Write();
  pThvPhDS2->Write();

  rThvPh->Write();
     
  rPid0->Write();
  sPid0->Write();
  secD0->Write();
  
  rPid1->Write();
  sPid1->Write();
  secD1->Write();

  pSecDS->Write();
  rSec->Write();

  sPid1_p->Write();
  sPid1_r->Write();
  
  outFile->cd("SeGA");

  coreEnergy->Write();
  coreSum->Write();
  segEnergy->Write();
  segSum->Write();

  coreEn_Fep->Write();
  coreEn_NotFep->Write();

  sPosTP->Write();

  outFile->cd("Coincidence/ProjectileDS");

  sPidDS->Write();
  rPidDS->Write();

  pCoreEnergyDS->Write();
  pCoreSumDS->Write();

  pCoreEnergyDS_fep->Write();
  pCoreEnergyDS_pfep->Write();
  pCoreEnergyDS_rfep->Write();
  pCoreEnergyDS_nfep->Write();

  outFile->cd("Coincidence/ProjectileDS/Doppler");
  
  pDopEnergyDS->Write();
  pDopSumDS->Write();

  pDopEnergyDS_fep->Write();
  pDopEnergyDS_pfep->Write();
  pDopEnergyDS_rfep->Write();
  pDopEnergyDS_nfep->Write();

  pDopvPartDS->Write();

  pThCorDS->Write();
  pThCrtDS->Write();

  pPhCorDS->Write();
  pPhCrtDS->Write();

  outFile->cd("Coincidence/ProjectileDS/Recon");

  pReconEnergyDS->Write();
  pReconSumDS->Write();

  pReconEnergyDS_fep->Write();
  pReconEnergyDS_pfep->Write();
  pReconEnergyDS_rfep->Write();
  pReconEnergyDS_nfep->Write();

  pReconvPartDS->Write();

  pReconThCorDS->Write();
  pReconThCrtDS->Write();

  pReconPhCorDS->Write();
  pReconPhCrtDS->Write();

  outFile->cd("Coincidence/ProjectileUS");

  sPidUS->Write();
  rPidUS->Write();

  pCoreEnergyUS->Write();
  pCoreSumUS->Write();

  pCoreEnergyUS_fep->Write();
  pCoreEnergyUS_pfep->Write();
  pCoreEnergyUS_rfep->Write();
  pCoreEnergyUS_nfep->Write();

  outFile->cd("Coincidence/ProjectileUS/Doppler");
  
  pDopEnergyUS->Write();
  pDopSumUS->Write();

  pDopEnergyUS_fep->Write();
  pDopEnergyUS_pfep->Write();
  pDopEnergyUS_rfep->Write();
  pDopEnergyUS_nfep->Write();

  pDopvPartUS->Write();

  pThCorUS->Write();
  pThCrtUS->Write();

  pPhCorUS->Write();
  pPhCrtUS->Write();

  outFile->cd("Coincidence/ProjectileUS/Recon");

  pReconEnergyUS->Write();
  pReconSumUS->Write();

  pReconEnergyUS_fep->Write();
  pReconEnergyUS_pfep->Write();
  pReconEnergyUS_rfep->Write();
  pReconEnergyUS_nfep->Write();

  pReconvPartUS->Write();

  pReconThCorUS->Write();
  pReconThCrtUS->Write();

  pReconPhCorUS->Write();
  pReconPhCrtUS->Write();
  
  outFile->cd("Coincidence/Recoil");

  sPidRec->Write();
  rPidRec->Write();

  rCoreEnergy->Write();
  rCoreSum->Write();

  rCoreEnergy_fep->Write();
  rCoreEnergy_pfep->Write();
  rCoreEnergy_rfep->Write();
  rCoreEnergy_nfep->Write();

  outFile->cd("Coincidence/Recoil/Doppler");
  
  rDopEnergy->Write();
  rDopSum->Write();

  rDopEnergy_fep->Write();
  rDopEnergy_pfep->Write();
  rDopEnergy_rfep->Write();
  rDopEnergy_nfep->Write();

  rDopvPart->Write();

  rThCor->Write();
  rThCrt->Write();

  rPhCor->Write();
  rPhCrt->Write();

  outFile->cd("Coincidence/Recoil/Recon");

  rReconEnergy->Write();
  rReconSum->Write();

  rReconEnergy_fep->Write();
  rReconEnergy_pfep->Write();
  rReconEnergy_rfep->Write();
  rReconEnergy_nfep->Write();

  rReconvPart->Write();

  rReconThCor->Write();
  rReconThCrt->Write();

  rReconPhCor->Write();
  rReconPhCrt->Write();

  outFile->Close();
  std::cout << "Done!" << std::endl;

  return 0;
}
