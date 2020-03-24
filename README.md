JANUS
=================================================================================
A GEANT4 simulation of low-energy Coulomb Excitation experiments with SeGA-JANUS.
=================================================================================

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

The correlator is essentially a stand-alone program compiled with ROOT libraries. To correctly sort the simulated data, the correlator.cc file needs to edited and recompiled. Directly after the inlcude statements, 5 variables need to be changed to match the simulation input: beamZ, beamA, beam_mass, targZ, targA, and targ_mass. Energy in MeV and mass in MeV/c<sup>2</sup>. 

To recompile the correlator, simply

- ./make_correlator.sh

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
/Output/Filename output.dat\
/run/beamOn nEvents

The /Mode command must come first, and the parameter (mode) can be Source, Scattering, or Full. The /Geometry/Update command is mandatory. The /Output/Filename command sets the name of the output data file. Example macros for each mode are in the Examples/Macros folder.

Geometry Commands
-----------------
The /Geometry commands are common across all modes. With the exception of /Geometry/Update, all /Geometry commands are optional.

| Command | Description |
| --- | --- |
| /Geometry/Bambino2/UpstreamOffset *double unit* | Set (positive) z-offset of upstream silicon detector. (Default: 3 cm) |
| /Geometry/Bambino2/DownstreamOffset *double unit* | Set (positive) z-offset of downstream silicon detector. (Default: 3 cm) |
| /Geometry/SeGA/Offset *double unit* | Set z-offset of SeGA. (Default: 0 cm) |
| /Geometry/Target/Z *int* | Set Z of target nucleus. (Default: 82) |
| /Geometry/Target/N *int* | Set N of target nucleus. (Default: 126) |
| /Geometry/Target/Density *double unit* | Set (volume) density of target material. (Default: 11.382 g/cm3) |
| /Geometry/Target/Mass *double unit* | Set mass of target material. (Default: 207.97665 g/mole) |
| /Geometry/Target/Thickness *double unit* | Set (linear) thickness of target. (Default: 882 nm) |
| /Geometry/Target/Radius *double unit* | Set radius of target. (Default: 0.5 cm) |
| /Geometry/Target/StandardTarget *string* | Set parameters for a standard target: 208Pb, 48Ti, or 196Pt. |
| /Geometry/Update | Update the simulation with your desired geometry. |

Note that the target does **NOT** define the recoiling nucleus for the kinematics or excitation, it only defines "bulk" material properties of the target.

Source Mode Commands
-----------------
There is only one /Source command, and it is mandatory.

| Command | Description |
| --- | --- |
| /Source/Energy *double unit* | Set energy of source gamma-rays. |

Scattering Mode Commands
-----------------
The Scattering mode commands are divided into two categories: /Beam and /Reaction. The /Beam commands define the properties of the incoming beam. The /Reaction commands control the kinematics. Each have mandatory and optional commands.

### Mandatory /Reaction commands
| Command | Description |
| --- | --- |
| /Reaction/ProjectileZ *int* | Set Z of projectile nucleus. |
| /Reaction/ProjectileA *int* | Set A of projectile nucleus. |
| /Reaction/RecoilZ *int* | Set Z of recoil nucleus. |
| /Reaction/RecoilA *int* | Set A of recoil nucleus. |

### Mandatory /Beam commands
| Command | Description |
| --- | --- |
| /Beam/Energy *double unit* | Set energy of incoming beam |

### Optional /Reaction commands
| Command | Description |
| --- | --- |
| /Reaction/SendToJanus | Only sample parts of the Rutherford scattering distribution which will result in a particle entering a silicon detector. |
| /Reaction/SendToUpstreamJanus | Only sample parts of the Rutherford scattering distribution which will result in a particle entering the upstream silicon detector. |
| /Reaction/SendToDownstreamJanus | Only sample parts of the Rutherford scattering distribution which will result in a particle entering the downstream silicon detector. |
| /Reaction/AddThetaLAB *double unit* | Add an angle to desired LAB scattering angle ranges. This command must always be used two at a time, with the smaller angle coming first. Otherwise it doesn't work. |
| /Reaction/OnlyProjectiles | Only consider the projectile when defining desired scattering angle ranges (above commands). |
| /Reaction/OnlyRecoils | Only consider the recoil when defining the desired scattering angle ranges (above commands). |
| /Reaction/DeltaE *double unit* | Set (positive) deltaE to simulate inelastic scattering. (Default: 0 MeV) |

The use of optional /Reaction commands is highly encouraged. Without these commands, the entire scattering angle range will be sampled according the Rutheford scattering distribution. This means roughly 10<sup>-4</sup> of your simulated events will result in a particle entring the silicon detectors.

The optional /Reaction commands can be used together, so you can fully customize what scattering angles you simulate. The /Reaction/SendToX commands will overwrite any scattering angles defined via the /Reaction/AddThetaLAB command. To use these commands together, always call the /Reaction/SendToX command first.

There are no "safety checks" for these commands. For example, never do

/Reaction/SendToUpstreamJanus\
/Reaction/OnlyRecoils

This is condition is never satisified and will set entire Rutherford distribtion to zero.

### Optional /Beam commands
| Command | Description |
| --- | --- |
| /Beam/SigmaEn *double unit* | Set Gaussian sigma of the kinetic energy distribution of the incoming beam (Default: 0 MeV) |
| /Beam/PositionX *double unit* | Set X position of incoming beam spot. (Default: 0 mm) |
| /Beam/PositionY *double unit* | Set Y position of incoming beam spot. (Default: 0 mm) |
| /Beam/AngleX" *double unit* | Set angle about x-axis of incoming beam. (Default: 0 deg) |
| /Beam/AngleY" *double unit* | Set angle about y-axis of incoming beam. (Default: 0 deg) |
| /Beam/SigmaX *double unit* | Set Gaussian sigma of x position distribution of the incoming beam. (Default: 0 mm) |
| /Beam/SigmaY *double unit* | Set Gaussian sigma of y position distribution of the incoming beam. (Default: 0 mm) |
| /Beam/SigmaAX *double unit* | Set Gaussian sigma of angle distribution about the x-axis. (Default: 0 deg) |
| /Beam/SigmaAY *double unit* | Set Gaussian sigma of angle distribution about the y-axis. (Default: 0 deg) |

Full Mode Commands
-----------------
For the Full CoulEx simulation, all Scattering mode commands still apply (except /Reaction/DeltaE). Additionally, /Excitation commands must be called which determine the level scheme and excitation pattern in both the projectile and recoil nucleus. No /Excitation commands are strictly optional or mandatory.

| Command | Description |
| --- | --- |
| /Excitation/Projectile/LevelScheme *string* | Name of the file to be read-in which defines the projectile level scheme. |
| /Excitation/Projectile/Probabilities *string* | Name of the file to be read-in which defines the angle-dependent excitation probabilities for the projectile. |
| /Excitation/Projectile/PopulateState *int* | Choose one state to populate in the projectile, irrespectie of scattering angle. This overrides /Excitation/Projectile/Probabilities. |
| /Excitation/Projectile/Simple | Declare a "simple" level scheme and excitation pattern for the projectile. This means there is only one excited state in the projectile, and it will always be populated. This overrides the above three commands. |
| /Excitation/Projectile/SimpleEnergy *double unit* | Set the energy of the simple state in the projectile. |
| /Excitation/Projectile/SimpleLifetime *double unit* | Set the lifetime of the simple state in the projectile. |

To control the recoil nucleus level scheme and excitations, replace /Excitation/Projectile/ with 
/Excitation/Recoil/. All commands are identically the same.

Level Scheme File Format
-----------------
The level scheme files are text files which describe the excited states of a nucleus. They have the following format.

II<sub>1</sub> en<sub>1</sub> tau<sub>1</sub> nb<sub>1</sub> \
&nbsp;IF<sup>(1)</sup><sub>1</sub> P<sup>(1)</sup><sub>1</sub>\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...\
&nbsp;IF<sup>(1)</sup><sub>nb<sub>1</sub></sub> P<sup>(1)</sup><sub>nb<sub>1</sub></sub>\
II<sub>2</sub> en<sub>2</sub> tau<sub>2</sub> nb<sub>2</sub> \
&nbsp;IF<sup>(2)</sup><sub>1</sub> P<sup>(2)</sup><sub>1</sub>\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...\
&nbsp;IF<sup>(2)</sup><sub>nb<sub>2</sub></sub> P<sup>(2)</sup><sub>nb<sub>2</sub></sub>\
...\
II<sub>N</sub> en<sub>N</sub> tau<sub>N</sub> nb<sub>N</sub> \
&nbsp;IF<sup>(N)</sup><sub>1</sub> P<sup>(N)</sup><sub>1</sub>\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;...\
&nbsp;IF<sup>(N)</sup><sub>nb<sub>N</sub></sub> P<sup>(N)</sup><sub>nb<sub>N</sub></sub>
 
Here II<sub>i</sub> is the index of the i-th level, en<sub>i</sub> is its energy in MeV, tau<sub>i</sub> is its mean-lifetime in ps, and nb<sub>i</sub> is the number of gamma decays from this state. IF<sup>(i)</sup><sub>j</sub> is the index of the final state for the j-th gamma decay of the i-th state. P<sup>(i)</sup><sub>j</sub> is the probability of that decay.

The states must be declared in order, i.e. II<sub>1</sub> = 1, II<sub>2</sub> = 2 and so on. This technically makes the initial state index redundant. The ground state (index 0) is not included in the level scheme file. There is no limit on the number of excited states or decays from a state. An example level scheme file is in the Examples/LevelSchemes folder.

Probability File Format
-----------------
The probability files are text files which describe the scattering-angle dependent excitation probabilitis of the excited states (defined in the level scheme file). They have the following format.


theta<sub>1</sub> P<sub>0</sub>(theta<sub>1</sub>) P<sub>1</sub>(theta<sub>1</sub>) ... P<sub>N</sub>(theta<sub>1</sub>)\
theta<sub>2</sub> P<sub>0</sub>(theta<sub>2</sub>) P<sub>1</sub>(theta<sub>2</sub>) ... P<sub>N</sub>(theta<sub>2</sub>)\
...\
theta<sub>K</sub> P<sub>0</sub>(theta<sub>K</sub>) P<sub>1</sub>(theta<sub>K</sub>) ... P<sub>N</sub>(theta<sub>K</sub>)

Here theta<sub>k</sub> is the center-of-mass frame scattering angle in radians. P<sub>i</sub>(theta<sub>k</sub>) is the excitation probability of i-th state for the CoM scattering angle theta<sub>k</sub>. The scattering angles must be entered smallest to largest, and there is no limit on the number of theta spline points. The state indices are defined by the level scheme file. Note that the ground state probabilities (index 0) must be included here. An example probabilities file is in the Examples/Probabilities folder.
