JANUS
===========================================================================================
A GEANT4 simulation of low-energy Coulomb Excitation experiments with the SeGA-JANUS setup.
===========================================================================================

Requirements
------------------
- JANUS requires GEANT4 v10.04.02. It has not been tested on any other version.
- In order to unpack, correlate, and histogram the simulation output, a ROOT instatllation is required.

Running JANUS
-----------------
JANUS takes an input macro file and writes the output to a data file. To run the simulation and subsequently histogram the simulation output, do

- JANUS input.mac
- correlator output.dat hist_file.root

The ROOT file hist_file.root now contains many histograms and can be opened with any standard ROOT installation.

Functionality
-----------------
The JANUS simultation has three modes: Source, Scattering, and Full. Each mode has different functionality and input requirements. 

- Source: Simulates a simple gamma-ray soure. No (massive) particles involved.
- Scattering: Simulates two-body scattering events. No gamma-rays involved. 
- Full: Simulates Coulomb Excitation events with user-definable level schemes and scattering-angle dependent excitation probabilities.

Macro Files and Commands
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
The /Geometry commands are common across all modes. With the exception of /Geometry/Update, all /Geometry commands are optional. Note that the target material does **NOT** define the recoiling nucleus for the kinematics or excittion, it only defines "bulk" material properties of the target.

- /Geometry/Bambino2/UpstreamOffset *double unit*
  - Set (positive) z-offset of upstream silicon detector (Default: 3 cm)
- /Geometry/Bambino2/DownstreamOffset *double unit*
  - Set (positive) z-offset of downstream silicon detector (Default: 3 cm)
- /Geometry/SeGA/Offset *double unit*
  - Set z-offset of SeGA array (Default: 0 cm)
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

Source Mode Commands
-----------------
There is only one /Source command, and it is mandatory.

- /Source/Energy *double unit*
  - Set energy of source gamma-rays

Scattering Mode Commands
-----------------
The Scattering mode commands are divided into two categories: /Beam and /Reaction. The /Beam commands define the properties of the incoming beam. The /Reaction commands control the kinematics. Each have mandatory and optional commands for a Scattering mode simulation.

##### Mandatory /Beam commands
- /Beam/Energy *double unit*
  - Set energy of incoming beam

##### Mandatory /Reaction commands
- /Reaction/ProjectileZ *int*
  - Set Z of projectile nucleus
- /Reaction/ProjectileA *int*
  - Set A of projectile nucleus



