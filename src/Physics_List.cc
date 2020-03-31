#include "Physics_List.hh"

#include "G4GenericIon.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4PhysicsListHelper.hh"
#include "G4StepLimiter.hh"

#include "G4ionIonisation.hh"

#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eBremsstrahlung.hh"

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
  ph->RegisterProcess(new G4StepLimiter(),G4GenericIon::Definition());
  
  //Gamma interactions
  ph->RegisterProcess(new G4PhotoElectricEffect(),G4Gamma::Definition());
  ph->RegisterProcess(new G4ComptonScattering(),G4Gamma::Definition());
  ph->RegisterProcess(new G4GammaConversion(),G4Gamma::Definition());

  //Electron
  ph->RegisterProcess(new G4eIonisation(),G4Electron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Electron::Definition());
  ph->RegisterProcess(new G4eMultipleScattering(),G4Electron::Definition());

  //Positron
  ph->RegisterProcess(new G4eIonisation(),G4Positron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Positron::Definition());
  ph->RegisterProcess(new G4eMultipleScattering(),G4Positron::Definition());
  ph->RegisterProcess(new G4eplusAnnihilation(),G4Positron::Definition());
  
}

void Physics_List::ConstructParticle() {

  G4Gamma::Definition();
  G4Proton::Definition();
  G4Electron::Definition();
  G4Positron::Definition();
  
  G4GenericIon::Definition();
  
}

