#include "Physics_List.hh"

#include "G4GenericIon.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4PhysicsListHelper.hh"
#include "G4ionIonisation.hh"

#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eBremsstrahlung.hh"
//#include "G4CoulombScattering.hh"
//#include "G4eCoulombScatteringModel.hh"
//#include "G4WentzelVIModel.hh"
//#include "G4UrbanMscModel96.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eplusAnnihilation.hh"

Physics_List::Physics_List() {}

Physics_List::~Physics_List() {}


void Physics_List::ConstructProcess() {

  AddTransportation();
  
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //Heavy ion interactions
  ph->RegisterProcess(new G4ionIonisation(),G4GenericIon::Definition());
  //ph->RegisterProcess(new G4hMultipleScattering("ionmsc"),G4GenericIon::Definition());

  //Gamma interactions
  ph->RegisterProcess(new G4PhotoElectricEffect(),G4Gamma::Definition());
  ph->RegisterProcess(new G4ComptonScattering(),G4Gamma::Definition());
  ph->RegisterProcess(new G4GammaConversion(),G4Gamma::Definition());

  //Scattering models (e+ and e-)
  G4eMultipleScattering* ms = new G4eMultipleScattering();
  //G4UrbanMscModel95* msm1 = new G4UrbanMscModel95();
  //G4WentzelVIModel* msm2 = new G4WentzelVIModel();
  //msm1->SetHighEnergyLimit(100*MeV);
  //msm2->SetLowEnergyLimit(100*MeV);
  //ms->AddEmModel(0,msm1);
  //ms->AddEmModel(0,msm2);

  //G4eCoulombScatteringModel* csm = new G4eCoulombScatteringModel(); 
  //G4CoulombScattering* cs = new G4CoulombScattering();
  //cs->SetEmModel(ssm, 1); 
  //cs->SetMinKinEnergy(100*MeV);
  //csm->SetLowEnergyLimit(100*MeV);
  //csm->SetActivationLowEnergyLimit(100*MeV);

  //Electron
  ph->RegisterProcess(new G4eIonisation(),G4Electron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Electron::Definition());
  ph->RegisterProcess(ms,G4Electron::Definition());
  //ph->RegisterProcess(cs,G4Electron::Definition());

  //Positron
  ph->RegisterProcess(new G4eIonisation(),G4Positron::Definition());
  ph->RegisterProcess(new G4eplusAnnihilation(),G4Positron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Positron::Definition());
  ph->RegisterProcess(ms,G4Positron::Definition());
  //ph->RegisterProcess(cs,G4Positron::Definition());
  
}

void Physics_List::ConstructParticle() {

  G4Gamma::Definition();
  G4Proton::Definition();
  G4Electron::Definition();
  G4Positron::Definition();
  
  G4GenericIon::Definition();
  
}

