#include "Primary_Generator.hh"

#include "G4GenericIon.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

#include "G4RunManager.hh"
#include "Detector_Construction.hh"
#include "G4EmCalculator.hh"

Primary_Generator::Primary_Generator() {

  messenger = new Primary_Generator_Messenger(this);
  gun = new G4ParticleGun(1);
  
  reac = new Reaction();
  source = new Gamma_Source();
  excite = new Excitation();

  projGS = NULL;
  recoilGS = NULL;

  dedx = 0.0*(MeV/mm);
  width = 0.0*mm;

  beam_X = 0.0*mm;
  beam_Y = 0.0*mm;
  beam_AX = 0.0*rad;
  beam_AY = 0.0*rad;
  beam_En = 0.0*MeV;
  
  sigma_X = 0.0*mm;
  sigma_Y = 0.0*mm;
  sigma_AX = 0.0*rad;
  sigma_AY = 0.0*rad;
  sigma_En = 0.0*MeV;

  deltaE = 0.0*MeV;
  
}

Primary_Generator::~Primary_Generator() {
  
  delete gun;
  delete messenger;
  delete reac;
  delete source;
  delete excite;
  
}

void Primary_Generator::GeneratePrimaries(G4Event* evt) {
  
  switch(mode) {
    
    case MODE::Scattering: {
 
      GenerateScatteringPrimaries(evt);
      /*
      G4ThreeVector bdir = G4RandomDirection();
      G4ThreeVector rdir = G4RandomDirection();
      //rdir.setTheta(bdir.theta());
      rdir.setPhi(bdir.phi()-pi);

      G4ThreeVector pos = G4ThreeVector(0.3*mm,0.63*mm,0.0);
      
      //Beam vertex
      gun->SetParticleDefinition(projGS);
      gun->SetParticleEnergy(200.0*MeV);
      gun->SetParticleMomentumDirection(bdir);
      gun->SetParticlePosition(pos);
      gun->GeneratePrimaryVertex(evt);
      
      //Recoil vertex
      gun->SetParticleDefinition(recoilGS);
      gun->SetParticleEnergy(100.0*MeV);
      gun->SetParticleMomentumDirection(rdir);
      gun->SetParticlePosition(pos);
      gun->GeneratePrimaryVertex(evt);
      */
      break;
    }

    case MODE::Source: {
      
      GenerateSourcePrimaries(evt);

      break;
    }

    case MODE::Full: {
      
      GenerateFullPrimaries(evt);
      
      break;
    }

  }

  return;
  
}

void Primary_Generator::GenerateScatteringPrimaries(G4Event* evt) {

  //Choose thetaCM 
  G4double th = reac->SampleRutherfordCM();
  
  //Randomize incoming beam energy using energy distribution
  G4double en = G4RandGauss::shoot(beam_En,sigma_En);

  //Randomly choose reaction depth in target
  G4double depth = G4RandFlat::shoot(width);

  //Energy loss
  en -= dedx*depth;
  
  //Reaction position
  //Randomize X and Y
  G4ThreeVector pos = G4ThreeVector(G4RandGauss::shoot(beam_X,sigma_X),
				    G4RandGauss::shoot(beam_Y,sigma_Y),
				    -(width/2.0) + depth); //reaction position
  //Outgoing vectors
  G4ThreeVector bdir = G4ThreeVector(0,0,1); //projectile direction
  bdir.setTheta(reac->Theta_LAB(th,en,deltaE)); //theta from kinematics
  bdir.setPhi(G4RandFlat::shoot(-pi,pi)); //randomly choose phi

  G4ThreeVector rdir = G4ThreeVector(0,0,1); //recoil direction
  rdir.setTheta(reac->Recoil_Theta_LAB(th,en,deltaE)); //theta from kinematics
  rdir.setPhi(bdir.phi()-pi); //Particles emerge back-to-back in LAB frame
  
  //Randomize direction using angle distributions
  G4double ax = G4RandGauss::shoot(beam_AX,sigma_AX);
  G4double ay = G4RandGauss::shoot(beam_AY,sigma_AY);
  
  bdir.rotateX(ax);
  bdir.rotateY(ay);
  rdir.rotateX(ax);
  rdir.rotateY(ay);
  
  //Beam vertex
  gun->SetParticleDefinition(projGS);
  gun->SetParticleEnergy(reac->KE_LAB(th,en,deltaE));
  gun->SetParticlePosition(pos);
  gun->SetParticleMomentumDirection(bdir);
  gun->GeneratePrimaryVertex(evt);
  
  //Recoil vertex
  gun->SetParticleDefinition(recoilGS);
  gun->SetParticleEnergy(reac->Recoil_KE_LAB(th,en,deltaE));
  gun->SetParticlePosition(pos);
  gun->SetParticleMomentumDirection(rdir);
  gun->GeneratePrimaryVertex(evt);

}

void Primary_Generator::GenerateSourcePrimaries(G4Event* evt) {

  //Simple isotropic gamma
  if(source->GetEnergy() > 0.0*MeV) {

    gun->SetParticleDefinition(G4Gamma::Definition());
    gun->SetParticlePosition(G4ThreeVector());
    gun->SetParticleEnergy(source->GetEnergy());
    gun->SetParticleMomentumDirection(G4RandomDirection());
    gun->GeneratePrimaryVertex(evt);

    return;
  }

  //Remove polarization
  source->Unpolarize();

  //Choose excited state
  G4int state_index = source->ChooseState();

  //Make vertex
  gun->SetParticleDefinition(source->GetDefinition(state_index));
  gun->SetParticleEnergy(0.0*MeV);
  gun->SetParticlePosition(G4ThreeVector());
  gun->SetParticleMomentumDirection(G4ThreeVector());
  gun->GeneratePrimaryVertex(evt);
  
  return;
}

void Primary_Generator::GenerateFullPrimaries(G4Event* evt) {

  //Remove all polarization
  excite->Unpolarize();
  
  //Choose thetaCM 
  G4double th = reac->SampleRutherfordCM();

  //Choose excited states
  G4int pI = excite->ChooseProjectileState(th);
  G4int rI = excite->ChooseRecoilState(th);

  //DelatE for inelastic scattering
  G4double ex = 0.0*MeV;
  ex += excite->GetProjectileExcitation(pI);
  ex += excite->GetRecoilExcitation(rI);

  //Randomize incoming beam energy using energy distribution
  G4double en = G4RandGauss::shoot(beam_En,sigma_En);

  //Randomly choose reaction depth in target
  G4double depth = G4RandFlat::shoot(width);

  //Energy loss
  en -= dedx*depth;

  //Reaction position
  //Randomize X and Y
  G4ThreeVector pos = G4ThreeVector(G4RandGauss::shoot(beam_X,sigma_X),
				    G4RandGauss::shoot(beam_Y,sigma_Y),
				    -(width/2.0) + depth); //reaction position
  
  //Outgoing vectors
  G4ThreeVector bdir = G4ThreeVector(0,0,1); //projectile direction
  bdir.setTheta(reac->Theta_LAB(th,en,ex)); //theta from kinematics
  bdir.setPhi(G4RandFlat::shoot(-pi,pi)); //randomly choose phi
  
  G4ThreeVector rdir = G4ThreeVector(0,0,1); //recoil direction
  rdir.setTheta(reac->Recoil_Theta_LAB(th,en,ex)); //theta from kinematics
  rdir.setPhi(bdir.phi()-pi); //Particles emerge back-to-back in LAB frame

  //Polarize excited states
  excite->Polarize(pI,rI,th,bdir.phi());
  
  //Randomize direction using angle distributions
  G4double ax = G4RandGauss::shoot(beam_AX,sigma_AX);
  G4double ay = G4RandGauss::shoot(beam_AY,sigma_AY);
  
  bdir.rotateX(ax);
  bdir.rotateY(ay);
  rdir.rotateX(ax);
  rdir.rotateY(ay);
  
  //Beam vertex
  gun->SetParticleDefinition(excite->GetProjectileDefinition(pI));
  gun->SetParticleEnergy(reac->KE_LAB(th,en,ex));
  gun->SetParticlePosition(pos);
  gun->SetParticleMomentumDirection(bdir);
  gun->GeneratePrimaryVertex(evt);
  
  //Recoil vertex
  gun->SetParticleDefinition(excite->GetRecoilDefinition(rI));
  gun->SetParticleEnergy(reac->Recoil_KE_LAB(th,en,ex));
  gun->SetParticlePosition(pos);
  gun->SetParticleMomentumDirection(rdir);
  gun->GeneratePrimaryVertex(evt);
  
}

void Primary_Generator::Update() {

  //Don't change the order of these functions
  switch(mode) {

    case MODE::Scattering: {
      
      G4IonTable* table = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
      
      projGS = table->GetIon(reac->GetBeamZ(),reac->GetBeamA(),0.0*MeV);
      recoilGS = table->GetIon(reac->GetRecoilZ(),reac->GetRecoilA(),0.0*MeV);
      projGS->SetPDGStable(true);
      recoilGS->SetPDGStable(true);

      UpdateReaction();
      
      break;
    }

    case MODE::Source: {
      
      source->BuildLevelScheme();
      
      break;
    }

    case MODE::Full: {

      G4int bZ = reac->GetBeamZ();
      G4int bA = reac->GetBeamA();
      G4int rZ = reac->GetRecoilZ();
      G4int rA = reac->GetRecoilA();
      
      excite->BuildLevelSchemes(bZ,bA,rZ,rA);
      excite->BuildProbabilities();
      
      projGS = excite->GetProjectileDefinition(0);
      recoilGS = excite->GetRecoilDefinition(0);
      UpdateReaction();

      G4double bEn = beam_En - 0.5*dedx*width;
      excite->BuildStatisticalTensors(bZ,bA,bEn,rZ,rA);
      
      break;
    }
  }
  
  return;
  
}

void Primary_Generator::UpdateReaction() {

  reac->SetBeamMass(projGS->GetPDGMass());
  reac->SetRecoilMass(recoilGS->GetPDGMass());
  reac->ConstructRutherfordCM(beam_En,deltaE);

  Detector_Construction* con =
    (Detector_Construction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

  if(con->GetTargetMaterial()) {
    G4EmCalculator* calc = new G4EmCalculator();
    dedx = calc->ComputeTotalDEDX(beam_En,projGS,con->GetTargetMaterial());
    width = con->GetTargetThickness();
    
  }
  
  return;
}

void Primary_Generator::SetMode(G4String md) {

  if(md == "Scattering") {
    mode = MODE::Scattering;
  }
  else if(md == "Source") {
    mode = MODE::Source;
  }
  else if(md == "Full") {
    mode = MODE::Full;
  }
  
  return;
  
}

