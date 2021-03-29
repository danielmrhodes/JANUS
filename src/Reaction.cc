#include "Reaction.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "Detector_Construction.hh"

Reaction::Reaction() : onlyP(false), onlyR(false), rDSpUS(false) {

  messenger = new Reaction_Messenger(this);

  recThresh = 0.0*MeV;
  
}

Reaction::Reaction(G4int bZ, G4int bA, G4double bM, G4int tZ, G4int tA, G4double tM)
  : beamZ(bZ), beamA(bA), beam_mass(bM), targZ(tZ), targA(tA), targ_mass(tM),
    onlyP(false), onlyR(false), rDSpUS(false) {

  messenger = new Reaction_Messenger(this);

  recThresh = 0.0*MeV;

}

Reaction::~Reaction() {

  delete messenger;
  
}

void Reaction::ConstructRutherfordCM(G4double Ep, G4double Ex) {

  G4cout << "\nBuilding Rutherford scattering distribution!\n" << G4endl;
  
  if(good_LAB_thetas.size()) {
    
    if(good_LAB_thetas.size()%2 == 0) {
      TruncatedRutherfordCM(Ep,Ex);
    }
    else {
      
      G4cout << "Desired LAB angle ranges were improperly defined! Resorting to default angle ranges!" << G4endl;
      
      Bambino2Thetas();
      TruncatedRutherfordCM(Ep,Ex);
    }
    
  }
  else {
    FullRutherfordCM(Ep,Ex);
  }

  return;
}

void Reaction::FullRutherfordCM(G4double Ep, G4double Ex) {

  const G4int nBins = 1000;
  G4double probDist[nBins-1];

  for(int i=0;i<nBins-1;i++) {

    G4double thetaCM = (pi/(double)nBins)*(i+1);
    probDist[i] = 2.0*pi*std::sin(thetaCM)*RutherfordCM(thetaCM,Ep,Ex);

  }

  AngleGenerator = new CLHEP::RandGeneral(probDist,nBins-1);

  return;
  
}

void Reaction::TruncatedRutherfordCM(G4double Ep, G4double Ex) {
  
  const G4int nBins = 1000;
  G4double probDist[nBins-1];

  for(int i=0;i<nBins-1;i++) {

    G4double thetaCM = (pi/(double)nBins)*(i+1);

    if(KeepThetaCM(thetaCM,Ep,Ex)) {
      probDist[i] = 2.0*pi*std::sin(thetaCM)*RutherfordCM(thetaCM,Ep,Ex);
    }
    else {
      probDist[i] = 0;
    }

  }
  
  AngleGenerator = new CLHEP::RandGeneral(probDist,nBins-1);

  return;
  
}

bool Reaction::KeepThetaCM(G4double thetaCM, G4double Ep, G4double Ex) {

  bool keep = false;

  if(onlyP) { //Only consider projectiles
    
    for(unsigned int i=0;i<good_LAB_thetas.size();i+=2) {
  
      if(
	 //Check if projectile will scatter into desired  range
         ((Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
          (Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)))
        ) {
      
        keep = true;
      } 
    }
    
  }
  else if(onlyR) { //Only consider recoils

    for(unsigned int i=0;i<good_LAB_thetas.size();i+=2) {
  
      if(
	 //Check if recoil will will scatter into desired  range
         ((Recoil_Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
          (Recoil_Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)) &&
	  (Recoil_KE_LAB(thetaCM,Ep,Ex) > recThresh))
        ) {
      
        keep = true;
      }
    }

  }
  else if(rDSpUS) { //Recoil downstream, projectile upstream

    for(unsigned int i=0;i<good_LAB_thetas.size();i+=2) {
  
      if(!i && //DS thetas only for recoil
	 //Check if recoil will will scatter into desired  range
         ((Recoil_Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
          (Recoil_Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)) &&
	  (Recoil_KE_LAB(thetaCM,Ep,Ex) > recThresh))
        ) {
      
        keep = true;
      }
      else if(i && //US thetas only for projectile
	      //Check if projectile will scatter into desired  range
	      ((Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
	       (Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)))
	     ) {
      
        keep = true;
      }
    }
    
  }
  else { //Consider both
    
    for(unsigned int i=0;i<good_LAB_thetas.size();i+=2) {
  
      if(
	 //Check if projectile will scatter into desired  range
         ((Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
          (Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)))
       
         ||
	 
         //Check if recoil will will scatter into desired  range
         ((Recoil_Theta_LAB(thetaCM,Ep,Ex) > good_LAB_thetas.at(i)) &&
          (Recoil_Theta_LAB(thetaCM,Ep,Ex) < good_LAB_thetas.at(i+1)) &&
	  (Recoil_KE_LAB(thetaCM,Ep,Ex) > recThresh))
       
        ) {	
      
        keep = true;
      }
    }
    
  }

  return keep;
}

void Reaction::Bambino2Thetas() {

  good_LAB_thetas.clear();

  //G4RunManager* runManager = G4RunManager::GetRunManager();
  Detector_Construction* con = (Detector_Construction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  //Bambino2 offsets
  G4double uso = con->GetUS_Offset();
  G4double dso = con->GetDS_Offset();

  //Bambino2 thickness
  G4double th = 300.0*um;

  G4ThreeVector v(0,0,1);

  //Downstream inner edge
  v.setRhoPhiZ(1.1*cm,0.0*rad,dso+(th/2.0));
  //v.setRhoPhiZ(1.1*cm - 2.2*mm,0.0*rad,dso+(th/2.0));
  good_LAB_thetas.push_back(v.theta());

  //Downstream outer edge
  v.setRhoPhiZ(3.5*cm,0.0*rad,dso-(th/2.0));
  //v.setRhoPhiZ(3.5*cm + 2.2*mm,0.0*rad,dso-(th/2.0));
  good_LAB_thetas.push_back(v.theta());

  //Upstream outer edge
  v.setRhoPhiZ(3.5*cm,0.0*rad,-uso+(th/2.0));
  good_LAB_thetas.push_back(v.theta());
  
  //Upstream inner edge
  v.setRhoPhiZ(1.1*cm,0.0*rad,-uso-(th/2.0));
  good_LAB_thetas.push_back(v.theta());

  return;
}

void Reaction::UpstreamThetas() {

  Bambino2Thetas();
  
  good_LAB_thetas.erase(good_LAB_thetas.begin());
  good_LAB_thetas.erase(good_LAB_thetas.begin());
  
  return;
}

void Reaction::DownstreamThetas() {

  Bambino2Thetas();
  
  good_LAB_thetas.erase(good_LAB_thetas.begin()+2);
  good_LAB_thetas.erase(good_LAB_thetas.begin()+2); 
  
  return;
}

void Reaction::AddThetaLAB(G4double T) {
  
  good_LAB_thetas.push_back(T);

  return;
}

void Reaction::SetOnlyP() {

  onlyR = false;
  onlyP = true;

  return;
}

void Reaction::SetOnlyR() {

  onlyP = false;
  onlyR = true;

  return;
}

/*
void Reaction::AddThetaLABRange(G4double Tmin, G4double Tmax) {

  good_LAB_thetas.push_back(Tmin);
  good_LAB_thetas.push_back(Tmax);

  return;
}
*/

G4double Reaction::SampleRutherfordCM() {

  return AngleGenerator->shoot()*pi;
  
}

//Everything below here come the GOSIA manual, chapter 5
/*
G4double Reaction::Theta_LAB_Max(G4double Ep, G4double Ex) {

  G4double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  if(tau < 1.0) {
    return pi;
  }
  else {
    return std::asin(1.0/tau);
  }
  
}

G4double Reaction::Recoil_Theta_LAB_Max(G4double Ep, G4double Ex) {

  G4double tau_t = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  return std::asin(1.0/tau_t); 
}
*/

G4double Reaction::RutherfordCM(G4double thetaCM, G4double Ep, G4double Ex) {

  G4double Ecm = Ep/(1.0 + ((G4double)beamA/(G4double)targA));
  G4double Esym2 = std::pow(Ecm,1.5)*std::sqrt(Ecm-Ex);
  G4double denom = Esym2*std::pow(std::sin(thetaCM/2.0),4);
  G4double factor = 1.29596*MeV*MeV*(millibarn/sr);
  
  return factor*std::pow(beamZ*targZ,2)/denom;
  
}
  
G4double Reaction::Theta_LAB(G4double thetaCM, G4double Ep, G4double Ex) {

  G4double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  G4double tanTheta = std::sin(thetaCM)/(std::cos(thetaCM) + tau);

  if(tanTheta > 0) {
    return std::atan(tanTheta);
  }
  else {
    return std::atan(tanTheta) + pi;
  }
}

G4double Reaction::Recoil_Theta_LAB(G4double thetaCM, G4double Ep, G4double Ex) {

  G4double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  G4double tanTheta = std::sin(pi - thetaCM)/(std::cos(pi - thetaCM) + tau);
  
  return std::atan(tanTheta);
  
}

G4double Reaction::KE_LAB(G4double thetaCM, G4double Ep, G4double Ex) {

  G4double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  G4double term1 = std::pow(targ_mass/(beam_mass + targ_mass),2);
  G4double term2 = 1 + tau*tau + 2*tau*std::cos(thetaCM);
  G4double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

G4double Reaction::Recoil_KE_LAB(G4double thetaCM, G4double Ep, G4double Ex) {

  G4double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  G4double term1 = beam_mass*targ_mass/std::pow(beam_mass + targ_mass,2);
  G4double term2 = 1 + tau*tau + 2*tau*std::cos(pi - thetaCM);
  G4double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

/*
G4double Reaction::Theta_CM_FP(G4double ThetaLAB, G4double Ep, G4bool sol2, G4double Ex) {

  G4double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += pi;
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }
  else {
    return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + pi;
  }
  
}

G4double Reaction::Theta_CM_FR(G4double ThetaLAB, G4double Ep, G4bool sol2, G4double Ex) {

  G4double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(std::sin(ThetaLAB) > 1.0/tau) {
    ThetaLAB = std::asin(1.0/tau);

    if(ThetaLAB < 0) {
      ThetaLAB += pi;
    }

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(!sol2) {
    return pi - (std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB);
  }
  else {
    return -std::asin(tau*std::sin(-ThetaLAB)) - ThetaLAB;
  }
  
}
*/
