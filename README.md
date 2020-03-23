JANUS
===========================================================================================
A GEANT4 simulation of low-energy Coulomb Excitation experiments with the SeGA-JANUS setup.
===========================================================================================

Requirements
------------------
- JANUS requires GEANT4 v10.04.02. It has not been tested on any other version.
- In order to unpack, correlate, and histogram the simulation output, a ROOT installation is required.

Running JANUS
-----------------
JANUS takes an input macro file and writes the output to a data file. To run the simulation and subsequently histogram the simulation output, do

- JANUS input.mac
- correlator output.dat hist_file.root

The ROOT file hist_file.root now contains many histograms and can be opened with any standard ROOT installation.

Functionality
-----------------
The JANUS simultation has three modes: Source, Scattering, and Full. Each mode has different functionality and input requirements. 

- Source: Simulates a simple isotropic gamma-ray source. No (massive) particles involved.
- Scattering: Simulates two-body scattering events. No gamma-rays involved. 
- Full: Simulates Coulomb Excitation events with user-definable level schemes and scattering-angle dependent excitation probabilities.

Macro Files
-----------------
The three simulation modes require different commands in their macro files. However all share a common structure 

/Mode mode\
*(optional geometry commands)*\
/Geometry/Update\
*(mode specific commands)*\
/run/beamOn nEvents

The /Mode command must come first, and the parameter (mode) can be Source, Scattering, or Full. The /Geometry/Update command is mandatory. 

Geometry Commands
-----------------
The /Geometry commands are common across all modes. With the exception of /Geometry/Update, all /Geometry commands are optional.

| Command | Description |
| --- | --- |
| /Geometry/Bambino2/UpstreamOffset *double unit* | Set (positive) z-offset of upstream silicon detector (Default: 3 cm) |
| /Geometry/Bambino2/DownstreamOffset *double unit* | Set (positive) z-offset of downstream silicon detector (Default: 3 cm) |
| /Geometry/SeGA/Offset *double unit* | Set z-offset of SeGA (Default: 0 cm) |
- /Geometry/Target/Z *int*
  - Set Z of target nucleus (Default: 82)
- /Geometry/Target/N *int*
  - Set N of target nucleus (Default: 126)
- /Geometry/Target/Density *double unit*
  - Set (volume) density of target material (Default: 11.382 g/cm3)
- /Geometry/Target/Mass *double unit*
  - Set mass of target material (Default: 207.97665 g/mole)
- /Geometry/Target/Thickness *double unit*
  - Set (linear) thickness of target (Default: 882 nm)
- /Geometry/Target/Radius *double unit*
  - Set radius of target (Default: 0.5 cm)
- /Geometry/Target/StandardTarget *string*
  - Set parameters for a standard target: 208Pb, 48Ti, or 196Pt
- /Geometry/Update
  - Update the simulation with your desired geometry

Note that the target does **NOT** define the recoiling nucleus for the kinematics or excittion, it only defines "bulk" material properties of the target.

Source Mode Commands
-----------------
There is only one /Source command, and it is mandatory.

| Command | Description |
| --- | --- |
| /Source/Energy *double unit* | Set energy of source gamma-rays |

Scattering Mode Commands
-----------------
The Scattering mode commands are divided into two categories: /Beam and /Reaction. The /Beam commands define the properties of the incoming beam. The /Reaction commands control the kinematics. Each have mandatory and optional commands for a Scattering mode simulation.

### Mandatory /Reaction commands
| Command | Description |
| --- | --- |
| /Reaction/ProjectileZ *int* | Set Z of projectile nucleus |
| /Reaction/ProjectileA *int* | Set A of projectile nucleus |
| /Reaction/RecoilZ *int* | Set Z of recoil nucleus |
| /Reaction/RecoilA *int* | Set A of recoil nucleus |

### Mandatory /Beam commands
- /Beam/Energy *double unit*
  - Set energy of incoming beam

### Optional /Reaction commands
- /Reaction/SendToJanus
  - Only sample parts of the Rutherford scattering distribution which will result in a particle entering a silicon detector.
- /Reaction/SendToUpstreamJanus
  - Only sample parts of the Rutherford scattering distribution which will result in a particle entering the upstream silicon detector.
- /Reaction/SendToDownstreamJanus
  - Only sample parts of the Rutherford scattering distribution which will result in a particle entering the downstream silicon detector.
- /Reaction/AddThetaLAB *double unit*
  - Add an angle to desired LAB scattering angle ranges. This command must always be used two at a time, with the smaller angle coming first. Otherwise it doesn't work.
- /Reaction/OnlyProjectiles
  - Only consider the projectile when defining desired scattering angle ranges (above commands)
- /Reaction/OnlyRecoils
  - Only consider the recoil when defining the desired scattering angle ranges (above commands)
- /Reaction/DeltaE *double unit*
  - Set (positive) deltaE to simulate inelastic scattering (Default: 0 MeV)

The use of optional \Reaction commands is highly encouraged. Without these commands, the entire scattering angle range will be sampled according the Rutheford scattering distribution. This means roughly 10^-4 of your simulated events will result in a particle entring the silicon detectors.

The optional \Reaction commands can be used together, so you can fully customize what scattering angles you simulate. The /Reaction/SendToX commands will overwrite any scattering angles defined via the /Reaction/AddThetaLAB command. To use these commands together, always call the /Reaction/SendToX command first.

There are no "safety checks" for these commands. For example, never do

/Reaction/SendToUpstreamJanus\
/Reaction/OnlyRecoils

This is condition is never satisified and will set entire Rutherford distribtion to zero.

### Optional /Beam commands
- /Beam/SigmaEn *double unit*
  - Set Gaussian sigma of the kinetic energy distribution of the incoming beam (Default: 0 MeV)
- /Beam/PositionX *double unit*
  - Set X position of incoming beam spot (Default: 0 mm)
- /Beam/PositionY *double unit*
  - Set Y position of incoming beam spot (Default: 0 mm)
- /Beam/AngleX" *double unit*
  - Set angle about x-axis of incoming beam (Default: 0 deg)
- /Beam/AngleY" *double unit*
  - Set angle about y-axis of incoming beam (Default: 0 deg)
- /Beam/SigmaX *double unit* 
  - Set Gaussian sigma of x position distribution of the incoming beam (Default: 0 mm)
- /Beam/SigmaY *double unit*
  - Set Gaussian sigma of y position distribution of the incoming beam (Default: 0 mm)
- /Beam/SigmaAX *double unit*
  - Set Gaussian sigma of angle distribution about the x-axis (Default: 0 deg)
- /Beam/SigmaAY *double unit*
  - Set Gaussian sigma of angle distribution about the y-axis (Default: 0 deg)

Full Mode Commands
-----------------
For the Full CoulEx simulation, all Scattering mode commands still apply (except /Reaction/DeltaE). Additionally, /Excitation commands must be called which determine the level scheme and excitation pattern in both the projectile and recoiling nucleus. No /Excitation commands are strictly optional or mandatory.

- /Excitation/Projectile/LevelScheme *string*
  - Name of the file to be read-in which defines the projectile level scheme.
- /Excitation/Projectile/Probabilities *string*
  - Name of the file to be read-in which defines the angle-dependent excitation probabilities for the projectile.
- /Excitation/Projectile/PopulateState *int* 
  - Choose one state to populate in the projectile, irrespectie of scattering angle. This overrides /Excitation/Projectile/Probabilities. 
- /Excitation/Projectile/Simple
  - Declare a "simple" level scheme and excitation pattern for the projectile. This means there is only one excited state in the projectile, and it will always be populated. This overrides the above three commands.
- /Excitation/Projectile/SimpleEnergy *double unit*
  - Set the energy of the simple state in the projectile
- /Excitation/Projectile/SimpleLifetime *double unit*
  - Set the lifetime of the simple state in the projectile

To control the recoil nucleus level scheme and excitations, replace /Excitation/Projectile/ with 
/Excitation/Recoil/. All commands are identically the same.

Level Scheme File Format
-----------------
The level scheme files are simple text files with the following format.



