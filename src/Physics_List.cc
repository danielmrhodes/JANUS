#include "Physics_List.hh"

#include "G4SystemOfUnits.hh"
//#include "G4EmParameters.hh"

//#include "G4LossTableManager.hh"
//#include "G4UAtomicDeexcitation.hh"

#include "G4GenericIon.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4PhysicsListHelper.hh"
#include "G4StepLimiter.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
//#include "G4NuclearStopping.hh"
//#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
//#include "G4ePairProduction.hh"

/*
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4Generator2BS.hh"
#include "G4Generator2BN.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"
*/

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

Physics_List::Physics_List() {

  /*
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(9);
  param->SetMinEnergy(100*eV);
  param->SetMaxEnergy(10*TeV);
  param->SetLowestElectronEnergy(100*eV);
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetUseMottCorrection(true);          // use Mott-correction for e-/e+ msc gs
  param->SetMscStepLimitType(fUseSafetyPlus); // error-free stepping for e-/e+ msc gs
  param->SetMscSkin(3);                       // error-free stepping for e-/e+ msc gs
  param->SetMscRangeFactor(0.2);              // error-free stepping for e-/e+ msc gs
  //param->SetMuHadLateralDisplacement(true);
  //param->SetLatDisplacementBeyondSafety(true);
  //param->SetFluo(true);
  //param->SetAugerCascade(true);
  //SetPhysicsType(bElectromagnetic);
  */
  
}

Physics_List::~Physics_List() {}


void Physics_List::ConstructProcess() {

  AddTransportation();
  
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //Generic ion
  /*
  ph->RegisterProcess(new G4ionIonisation(),G4GenericIon::Definition());
  ph->RegisterProcess(new G4StepLimiter(),G4GenericIon::Definition());
  //ph->RegisterProcess(new G4hMultipleScattering("ionmsc"),G4GenericIon::Definition());
  */
  
  //Gamma
  /*
  ph->RegisterProcess(new G4PhotoElectricEffect(),G4Gamma::Definition());
  ph->RegisterProcess(new G4ComptonScattering(),G4Gamma::Definition());
  ph->RegisterProcess(new G4GammaConversion(),G4Gamma::Definition());
  ph->RegisterProcess(new G4RayleighScattering(),G4Gamma::Definition());
  */

  //Electron
  ph->RegisterProcess(new G4eIonisation(),G4Electron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Electron::Definition());
  ph->RegisterProcess(new G4eMultipleScattering(),G4Electron::Definition());

  //Positron
  ph->RegisterProcess(new G4eIonisation(),G4Positron::Definition());
  ph->RegisterProcess(new G4eBremsstrahlung(),G4Positron::Definition());
  ph->RegisterProcess(new G4eMultipleScattering(),G4Positron::Definition());
  ph->RegisterProcess(new G4eplusAnnihilation(),G4Positron::Definition());
  
  //Generic Ion
  //G4NuclearStopping* inuc = new G4NuclearStopping();
  
  G4ionIonisation* ionIoni = new G4ionIonisation();
  //ionIoni->SetEmModel(new G4IonParametrisedLossModel());
  //ionIoni->SetStepFunction(0.1, 1*um);

  //ph->RegisterProcess(new G4hMultipleScattering("ionmsc"),G4GenericIon::Definition());
  ph->RegisterProcess(ionIoni,G4GenericIon::Definition());
  //ph->RegisterProcess(inuc,G4GenericIon::Definition());
  ph->RegisterProcess(new G4StepLimiter(),G4GenericIon::Definition());

  //inuc->SetMaxKinEnergy(MeV);
  
  //Gamma
  G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
  pe->SetEmModel(new G4LivermorePhotoElectricModel());
  ph->RegisterProcess(pe,G4Gamma::Definition());
  //G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();

  // Compton scattering
  G4ComptonScattering* cs = new G4ComptonScattering;
  cs->SetEmModel(new G4KleinNishinaModel());
  
  G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
  theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
  cs->AddEmModel(0,theLowEPComptonModel);
  ph->RegisterProcess(cs,G4Gamma::Definition());

  // Gamma conversion
  G4GammaConversion* gc = new G4GammaConversion();
  G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
  
  thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
  gc->SetEmModel(thePenelopeGCModel);
  ph->RegisterProcess(gc,G4Gamma::Definition());

  // Rayleigh scattering
  ph->RegisterProcess(new G4RayleighScattering(),G4Gamma::Definition());

  /*
  //Electron and Positron
  G4double highEnergyLimit = 100*MeV; // energy limits for e+- scattering models
  G4double penEnergyLimit = 1*MeV; // energy limits for e+- ionisation models
  
  // multiple scattering
  G4eMultipleScattering* msc = new G4eMultipleScattering;
  // e-/e+ msc gs with Mott-correction
  // (Mott-correction is set through G4EmParameters)
  G4GoudsmitSaundersonMscModel* msc1 = new G4GoudsmitSaundersonMscModel();
  G4WentzelVIModel* msc2 = new G4WentzelVIModel();

  msc1->SetHighEnergyLimit(highEnergyLimit);
  msc2->SetLowEnergyLimit(highEnergyLimit);
  msc->SetEmModel(msc1);
  msc->SetEmModel(msc2);

  G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
  G4CoulombScattering* ss = new G4CoulombScattering();

  ss->SetEmModel(ssm); 
  ss->SetMinKinEnergy(highEnergyLimit);
  ssm->SetLowEnergyLimit(highEnergyLimit);
  ssm->SetActivationLowEnergyLimit(highEnergyLimit);

  // ionisation
  G4eIonisation* eIoni = new G4eIonisation();
  eIoni->SetStepFunction(0.2, 10*um);

  G4PenelopeIonisationModel* pen = new G4PenelopeIonisationModel();
  pen->SetHighEnergyLimit(penEnergyLimit);
  eIoni->AddEmModel(0,pen,new G4UniversalFluctuation());

  // bremsstrahlung
  G4eBremsstrahlung* brem = new G4eBremsstrahlung();
  G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
  G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();

  br1->SetAngularDistribution(new G4Generator2BS());
  br2->SetAngularDistribution(new G4Generator2BS());
  brem->SetEmModel(br1);
  brem->SetEmModel(br2);
  br2->SetLowEnergyLimit(GeV);

  // pair production
  G4ePairProduction* ee = new G4ePairProduction();

  // register processes
  ph->RegisterProcess(msc,G4Electron::Definition());
  ph->RegisterProcess(eIoni,G4Electron::Definition());
  ph->RegisterProcess(brem,G4Electron::Definition());
  ph->RegisterProcess(ee,G4Electron::Definition());
  ph->RegisterProcess(ss,G4Electron::Definition());

  ph->RegisterProcess(msc,G4Positron::Definition());
  ph->RegisterProcess(eIoni,G4Positron::Definition());
  ph->RegisterProcess(brem,G4Positron::Definition());
  ph->RegisterProcess(ee,G4Positron::Definition());
  ph->RegisterProcess(ss,G4Positron::Definition());
  ph->RegisterProcess(new G4eplusAnnihilation(),G4Positron::Definition());

  //G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  //G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  */

  return;
}

void Physics_List::ConstructParticle() {

  G4Gamma::Definition();
  G4Proton::Definition();
  G4Electron::Definition();
  G4Positron::Definition();
  
  G4GenericIon::Definition();

  return;
  
}

