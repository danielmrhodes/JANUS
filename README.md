JANUS
=================================================================================
A GEANT4 simulation of low-energy Coulomb excitation experiments with SeGA-JANUS.
=================================================================================

Requirements
------------------
- JANUS has been tested on both GEANT4 v10.04.02 and v10.07.00. It *should* compile on any version of GEANT4 after v10.04.02.
- In order to use the file correlator.cc to sort and histogram the simulated data, a ROOT installation is required.

Running JANUS
-----------------
JANUS takes an input macro file and writes the output to a data file. To run the simulation and subsequently histogram the simulation output, do

- JANUS input.mac
- correlator output.dat hist_file.root

The ROOT file hist_file.root now contains many histograms and can be opened with any standard ROOT installation.

The correlator is a small program compiled with ROOT libraries. To correctly sort the simulated data, the correlator.cc file needs to be edited and recompiled. Directly after the inlcude statements there are several variables which need to be changed to match the simulation input.

To recompile the correlator, simply

- ./make_correlator.sh

Functionality
-----------------
The JANUS simultation has three modes: Source, Scattering, and Full. Each mode has different functionality and input requirements. 

- Source: Simulates either a simple isotropic gamma-ray, or a user-definable gamma-ray cascade. No massive particles involved.
- Scattering: Simulates two-body scattering events. No gamma-rays involved. 
- Full: Simulates Coulomb Excitation events with user-definable level schemes, scattering-angle dependent excitation probabilities, and scattering-angle dependent alignment of the excited states..

Macro Files
-----------------
The three simulation modes require different commands in their macro files. However all share a common structure: 

<pre>
/Mode mode\
*(optional geometry commands)*\
/Geometry/Update\
*(mode specific commands)*\
/Output/Filename output.dat\
/run/beamOn nEvents
</pre>

The /Mode command must come first, and the parameter (mode) can be Source, Scattering, or Full. The /Geometry/Update command is mandatory. The /Output/Filename command sets the name of the output data file. Example macros for each mode are in the Examples/Macros folder.

Geometry Commands
-----------------
The /Geometry commands are common across all modes. With the exception of /Geometry/Update, all /Geometry commands are optional.

| Command | Description |
| --- | --- |
| /Geometry/Bambino2/Construct | Include the silicon detectors in the simulation |
| /Geometry/Target/Construct | Include the target in the simulation |
| /Geometry/Bambino2/UpstreamOffset *double unit* | Set (positive) z-offset of upstream silicon detector. (Default: 3 cm) |
| /Geometry/Bambino2/DownstreamOffset *double unit* | Set (positive) z-offset of downstream silicon detector. (Default: 3 cm) |
| /Geometry/SeGA/Offset *double unit* | Set z-offset of SeGA. (Default: 0 cm) |
| /Geometry/Target/StandardTarget *string* | Set parameters for a standard target: 208Pb, 48Ti, or 196Pt. |
| /Geometry/Target/Z *int* | Set Z of target nucleus. (Default: 82) |
| /Geometry/Target/N *int* | Set N of target nucleus. (Default: 126) |
| /Geometry/Target/Density *double unit* | Set (volume) density of target material. (Default: 11.382 g/cm<sup>3</sup>) |
| /Geometry/Target/Mass *double unit* | Set mass of target material. (Default: 207.97665 g/mole) |
| /Geometry/Target/Thickness *double unit* | Set (linear) thickness of target. (Default: 882 nm) |
| /Geometry/Target/Radius *double unit* | Set radius of target. (Default: 0.5 cm) |
| /Geometry/Update | Update the simulation with your desired geometry. |

Note that the target does **NOT** define the recoiling nucleus for the kinematics or excitation, it only defines "bulk" material properties of the target. If you do not call the /Construct commands for the silicon detectors or the target, they will not be in the simulation.

Source Mode Commands
-----------------
There are two /Source commands, one of which must be called. They are mutually exclusive.

| Command | Description |
| --- | --- |
| /Source/Energy *double unit* | Simulate a single isotropic gamma-ray of this energy |
| /Source/LevelScheme *string* | Simulate a gamma-ray cascade defined by this file |

See the Source Level Scheme File Format section below for details on the level scheme file. 

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
| /Reaction/RecDS_ProjUS | Send *only* the recoiling nucleus to the downstream silicon detector while also sending the projectile nucleus to the upstream detector. This should not be used with any of the above optional /Reaction commands. |
| /Reaction/DeltaE *double unit* | Set (positive) deltaE to simulate inelastic scattering. (Default: 0 MeV) |
 /Reaction/RecoilThreshold *double unit* | Set recoil energy detection threshold. Removes very low energy recoils, which scatter into the silicon detector, that occur at small CM angles during inelastic scattering. (Default: 0 MeV) |

The use of optional /Reaction commands is strongly recommended. Without these commands, the entire scattering angle range [0,pi] will be sampled according the Rutheford scattering distribution. This means roughly 10<sup>-4</sup> of your simulated events will result in a particle entering the silicon detectors.

The optional /Reaction commands can be used together, so you can fully customize what scattering angles you simulate. The /Reaction/SendToX commands will overwrite any scattering angles defined via the /Reaction/AddThetaLAB command. To use these commands together, always call the /Reaction/SendToX command first.

There are no "safety checks" for these commands. For example, never do

/Reaction/SendToUpstreamJanus\
/Reaction/OnlyRecoils

This condition is never satisified and will set the entire Rutherford distribution to zero.

### Optional /Beam commands
| Command | Description |
| --- | --- |
| /Beam/SigmaEn *double unit* | Set Gaussian sigma of the kinetic energy distribution of the incoming beam (Default: 0 MeV) |
| /Beam/PositionX *double unit* | Set X position of incoming beam spot. (Default: 0 mm) |
| /Beam/PositionY *double unit* | Set Y position of incoming beam spot. (Default: 0 mm) |
| /Beam/AngleX *double unit* | Set angle about x-axis of incoming beam. (Default: 0 deg) |
| /Beam/AngleY *double unit* | Set angle about y-axis of incoming beam. (Default: 0 deg) |
| /Beam/SigmaX *double unit* | Set Gaussian sigma of x position distribution of the incoming beam. (Default: 0 mm) |
| /Beam/SigmaY *double unit* | Set Gaussian sigma of y position distribution of the incoming beam. (Default: 0 mm) |
| /Beam/SigmaAX *double unit* | Set Gaussian sigma of angle distribution about the x-axis. (Default: 0 deg) |
| /Beam/SigmaAY *double unit* | Set Gaussian sigma of angle distribution about the y-axis. (Default: 0 deg) |

Full Mode Commands
-----------------
For the Full CoulEx simulation, all Scattering mode commands still apply. Additionally, /Excitation commands can be called which determine the level scheme, excitation pattern, and nuclear alignment in both the projectile and recoil nucleus. All /Excitation commands are optional. Not calling any /Excitation commands will reduce the Full simulation to a Scattering simulation.

| Command | Description |
| --- | --- |
| /Excitation/Projectile/LevelScheme *string* | Name of the file to be read-in which defines the projectile level scheme. |
| /Excitation/Projectile/Probabilities *string* | Name of the file to be read-in which defines the angle-dependent excitation probabilities for the projectile. |
| /Excitation/Projectile/StatisticalTensors *string* | Name of the file to be read-in which contains the statistical tensors for the projectile to determine the angle-dependent nuclear alignment. |
| /Excitation/Projectile/OnlyConsiderState *int* | Turn off gamma decays from all states except this one. This requires a level scheme file and a probabilities file. |
| /Excitation/Projectile/PopulateState *int* | Choose one state to populate in the projectile, irrespective of scattering angle. This only requires a level scheme file (the probabilities file will simply be ignored). |
| /Excitation/Projectile/GroundStateSpin *double* | Set the spin of the projectile ground state. The spin must be integer or half-integer. (Default: 0.0) |
| /DeorientationEffect/Projectile/CalculateGk *bool* | Control whether the deorientation effect coefficients G<sub>k</sub> will be calculated for the projectile. These attenuate the nuclear alignment induced after CoulEx, and are only used if the statistical tensors file is provided. (Default: true) |
| /DeorientationEffect/Projectile/AverageJ *double* | Set the average atomic spin in the projectile for the deorentation effect two-state model (Default: 3.0) |
| /DeorientationEffect/Projectile/Gamma *double* | Set the FWHM of the frequency distribution (ps^-1 ) in the projectile for the deorentation effect two-state model (Default: 0.02) |
| /DeorientationEffect/Projectile/Lambda *double* | Set the transition rate (ps^-1 ) between static and fluctuating states in the projectile for the deorentation effect two-state model (Default: 0.0345) |
| /DeorientationEffect/Projectile/TauC *double* | Set the correlation time (ps) of the fluctating state in the projectile for the deorentation effect two-state model (Default: 3.5) |
| /DeorientationEffect/Projectile/GFactor *double* | Set the gyromagnetic ratio (g-factor) of the projectile for the deorentation effect two-state model (Default: Z/A) |
| /DeorientationEffect/Projectile/FieldCoefficient *double* | Set the hyperfine field coefficient (10^8 T) in the projectile for the deorentation effect two-state model (Default: 6*10<sup>-6</sup>) |
| /DeorientationEffect/Projectile/FieldExponent *double* | Set the hyperfine field exponent in the projectile for the deorentation effect two-state model (Default: 0.6) |

The recoiling target nucleus can be controlled with identical commands. Simply replace /Projectile/ with /Recoil/ in any of the above commands.

The statistical tensors [1] and deorientation effect coefficients G<sub>k</sub> are critical for reproducing the LAB frame gamma-ray spectra. See [2] for details on the two-state model, and the meaning of its parameters, which is used to describe the deorientation effect.

*If you input a level scheme, you must also input the probabilities or choose a state to populate. Otherwise the level scheme won't be used.*

*In a Full simulation, the /Reaction/DeltaE command only affects what CM angles will be sampled. The Q-value and LAB scattering anlges for each event are calculated based on which excited states get populated.*

Input Preparation
-----------------
The ROOT script MakeInput.C will make the probabilities file and statistical tensors file that can be given to JANUS. The script uses the Coulomb excitation code Cygnus [3] for the calculations. The Cygnus libraries must be loaded into the ROOT session before loading MakeInput.C, and the Cygnus nucleus file must already be created. See [3] for details.

The level scheme file, for either a Source or Full simulation, must be created manually. The file has a different format depending on the simulation mode; these are described below.   

Full Mode Level Scheme File Format
-----------------

<pre>
II<sub>1</sub> En<sub>1</sub> Sp<sub>1</sub> Tau<sub>1</sub> Nb<sub>1</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>1</sub></sub> P<sub>Nb<sub>1</sub></sub> L1<sub>Nb<sub>1</sub></sub> L2<sub>Nb<sub>1</sub></sub> DL<sub>Nb<sub>1</sub></sub> CC<sub>Nb<sub>1</sub></sub>
II<sub>2</sub> En<sub>2</sub> Sp<sub>2</sub> Tau<sub>2</sub> Nb<sub>2</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>2</sub></sub> P<sub>Nb<sub>2</sub></sub> L1<sub>Nb<sub>2</sub></sub> L2<sub>Nb<sub>2</sub></sub> DL<sub>Nb<sub>2</sub></sub> CC<sub>Nb<sub>2</sub></sub>
...
II<sub>N</sub> En<sub>N</sub> Sp<sub>N</sub> Tau<sub>N</sub> Nb<sub>N</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>N</sub></sub> P<sub>Nb<sub>N</sub></sub> L1<sub>Nb<sub>N</sub></sub> L2<sub>Nb<sub>N</sub></sub> DL<sub>Nb<sub>N</sub></sub> CC<sub>Nb<sub>N</sub></sub>
</pre>
 
Here II<sub>i</sub> is the index of state i, En<sub>i</sub> is its energy in keV, Sp<sub>i</sub> is its spin (J), Tau<sub>i</sub> is its mean-lifetime in ps, and Nb<sub>i</sub> is the number of gamma decays from this state. IF<sub>j</sub> is the index of the final state for gamma decay j of this state. P<sub>j</sub> is the probability of that decay (relative to the other gamma-ray decays), L1<sub>j</sub> is the higher multipolarity of the transition, L2<sub>j</sub> is the lower multipolarity of the transition, and DL<sub>j</sub> is the mixing ratio. CC<sub>j</sub> is the total conversion coefficient of this gamma-ray transition. 

The states must be declared in order, i.e. II<sub>1</sub> = 1, II<sub>2</sub> = 2 and so on. This technically makes the initial state index redundant. The ground state (index 0) is not included in the level scheme file. There is no limit on the number of excited states or decays from a state.

The spin of a state must be integer or half-integer. For the gamma transitions, L1 > L2. If DL = 0, the transition will be pure L1 and L2 is ignored. Example level scheme files for a Full simulation are in the Examples/LevelSchemes/Full folder.

Source Mode Level Scheme File Format
-----------------

<pre>
II<sub>1</sub> En<sub>1</sub> Sp<sub>1</sub> Tau<sub>1</sub> Pop<sub>1</sub> Nb<sub>1</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>1</sub></sub> P<sub>Nb<sub>1</sub></sub> L1<sub>Nb<sub>1</sub></sub> L2<sub>Nb<sub>1</sub></sub> DL<sub>Nb<sub>1</sub></sub> CC<sub>Nb<sub>1</sub></sub>
II<sub>2</sub> En<sub>2</sub> Sp<sub>2</sub> Tau<sub>2</sub> Pop<sub>2</sub> Nb<sub>2</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>2</sub></sub> P<sub>Nb<sub>2</sub></sub> L1<sub>Nb<sub>2</sub></sub> L2<sub>Nb<sub>2</sub></sub> DL<sub>Nb<sub>2</sub></sub> CC<sub>Nb<sub>2</sub></sub>
...
II<sub>N</sub> En<sub>N</sub> Sp<sub>N</sub> Tau<sub>N</sub> Pop<sub>2</sub> Nb<sub>N</sub>
 IF<sub>1</sub> P<sub>1</sub> L1<sub>1</sub> L2<sub>1</sub> DL<sub>1</sub> CC<sub>1</sub>
 ...
 IF<sub>Nb<sub>N</sub></sub> P<sub>Nb<sub>N</sub></sub> L1<sub>Nb<sub>N</sub></sub> L2<sub>Nb<sub>N</sub></sub> DL<sub>Nb<sub>N</sub></sub> CC<sub>Nb<sub>N</sub></sub>
</pre>

The Source level scheme file is the same the Full level schem file (above), but has one additional entry, Pop<sub>i</sub>, which comes before Nb<sub>i</sub>. This is the relative population of the state i. An example level scheme file (co60.lvl) for a Source simulation is in the Examples/LevelSchemes/Source folder.

Probability File Format
-----------------
The probability files are text files which describe the scattering-angle dependent excitation probabilities of the excited states. They have the following format.

<pre>
theta<sub>1</sub> P<sub>0</sub>(theta<sub>1</sub>) P<sub>1</sub>(theta<sub>1</sub>) ... P<sub>N</sub>(theta<sub>1</sub>)
theta<sub>2</sub> P<sub>0</sub>(theta<sub>2</sub>) P<sub>1</sub>(theta<sub>2</sub>) ... P<sub>N</sub>(theta<sub>2</sub>)
...
theta<sub>K</sub> P<sub>0</sub>(theta<sub>K</sub>) P<sub>1</sub>(theta<sub>K</sub>) ... P<sub>N</sub>(theta<sub>K</sub>)
</pre>

Here theta<sub>k</sub> is the center-of-mass frame scattering angle in radians. There is no limit on the number of theta spline points. P<sub>i</sub>(theta<sub>k</sub>) is the excitation probability of state i for the CoM scattering angle theta<sub>k</sub>. The scattering angles must be entered smallest to largest, and the probabilities must be entered in the order of the states 0 to N. The state indices are defined by the level scheme file, and all states must be included. Note that the ground state probabilities (index 0) must be included here. An example probabilities file is in the Examples/Probabilities folder.

Statistical Tensor File Format
-----------------
The statistical tensor files are text files which contatin the components of the statistical tensor for each excited state at different scattering angles. They have the following format.

<pre>
Theta [CM]: theta<sub>1</sub> rad
STATISTICAL TENSORS: LAB FRAME
INDEX:	    KA:	     KAPPA:	RHOC:
1	    0	     0		rho<sup>(1)</sup><sub>00</sub>(theta<sub>1</sub>)
1	    2	     0		rho<sup>(1)</sup><sub>20</sub>(theta<sub>1</sub>)
1	    2	     1		rho<sup>(1)</sup><sub>21</sub>(theta<sub>1</sub>)
1	    2	     2		rho<sup>(1)</sup><sub>22</sub>(theta<sub>1</sub>)
...
1	    k<sup>(1)</sup><sub>max</sub>    k<sup>(1)</sup><sub>max</sub>	rho<sup>(1)</sup><sub>k<sup>(1)</sup><sub>max</sub>k<sup>(1)</sup><sub>max</sub></sub>(theta<sub>1</sub>)
2	    0	     0		rho<sup>(2)</sup><sub>00</sub>(theta<sub>1</sub>)
2	    2	     0		rho<sup>(2)</sup><sub>20</sub>(theta<sub>1</sub>)
2	    2	     1		rho<sup>(2)</sup><sub>21</sub>(theta<sub>1</sub>)
2	    2	     2		rho<sup>(2)</sup><sub>22</sub>(theta<sub>1</sub>)
...
2	    k<sup>(2)</sup><sub>max</sub>    k<sup>(2)</sup><sub>max</sub>	rho<sup>(2)</sup><sub>k<sup>(2)</sup><sub>max</sub>k<sup>(2)</sup><sub>max</sub></sub>(theta<sub>1</sub>)
...
N	    0	     0		rho<sup>(N)</sup><sub>00</sub>(theta<sub>1</sub>)
N	    2	     0		rho<sup>(N)</sup><sub>20</sub>(theta<sub>1</sub>)
N	    2	     1		rho<sup>(N)</sup><sub>21</sub>(theta<sub>1</sub>)
N	    2	     2		rho<sup>(N)</sup><sub>22</sub>(theta<sub>1</sub>)
...
N	    k<sup>(N)</sup><sub>max</sub>    k<sup>(N)</sup><sub>max</sub>	rho<sup>(N)</sup><sub>k<sup>(N)</sup><sub>max</sub>k<sup>(N)</sup><sub>max</sub></sub>(theta<sub>1</sub>)

Theta [CM]: theta<sub>2</sub> rad
STATISTICAL TENSORS: LAB FRAME
INDEX:	    KA:	     KAPPA:	RHOC:
1	    0	     0		rho<sup>(1)</sup><sub>00</sub>(theta<sub>2</sub>)
1	    2	     0		rho<sup>(1)</sup><sub>20</sub>(theta<sub>2</sub>)
1	    2	     1		rho<sup>(1)</sup><sub>21</sub>(theta<sub>2</sub>)
1	    2	     2		rho<sup>(1)</sup><sub>22</sub>(theta<sub>2</sub>)
...
1	    k<sup>(1)</sup><sub>max</sub>    k<sup>(1)</sup><sub>max</sub>	rho<sup>(1)</sup><sub>k<sup>(1)</sup><sub>max</sub>k<sup>(1)</sup><sub>max</sub></sub>(theta<sub>2</sub>)
2	    0	     0		rho<sup>(2)</sup><sub>00</sub>(theta<sub>2</sub>)
2	    2	     0		rho<sup>(2)</sup><sub>20</sub>(theta<sub>2</sub>)
2	    2	     1		rho<sup>(2)</sup><sub>21</sub>(theta<sub>2</sub>)
2	    2	     2		rho<sup>(2)</sup><sub>22</sub>(theta<sub>2</sub>)
...
2	    k<sup>(2)</sup><sub>max</sub>    k<sup>(2)</sup><sub>max</sub>	rho<sup>(2)</sup><sub>k<sup>(2)</sup><sub>max</sub>k<sup>(2)</sup><sub>max</sub></sub>(theta<sub>2</sub>)
...
N	    0	     0		rho<sup>(N)</sup><sub>00</sub>(theta<sub>2</sub>)
N	    2	     0		rho<sup>(N)</sup><sub>20</sub>(theta<sub>2</sub>)
N	    2	     1		rho<sup>(N)</sup><sub>21</sub>(theta<sub>2</sub>)
N	    2	     2		rho<sup>(N)</sup><sub>22</sub>(theta<sub>2</sub>)
...
N	    k<sup>(N)</sup><sub>max</sub>    k<sup>(N)</sup><sub>max</sub>	rho<sup>(N)</sup><sub>k<sup>(N)</sup><sub>max</sub>k<sup>(N)</sup><sub>max</sub></sub>(theta<sub>2</sub>)
...

Theta [CM]: theta<sub>K</sub> rad
STATISTICAL TENSORS: LAB FRAME
INDEX:	    KA:	     KAPPA:	RHOC:
1	    0	     0		rho<sup>(1)</sup><sub>00</sub>(theta<sub>K</sub>)
1	    2	     0		rho<sup>(1)</sup><sub>20</sub>(theta<sub>K</sub>)
1	    2	     1		rho<sup>(1)</sup><sub>21</sub>(theta<sub>K</sub>)
1	    2	     2		rho<sup>(1)</sup><sub>22</sub>(theta<sub>K</sub>)
...
1	    k<sup>(1)</sup><sub>max</sub>    k<sup>(1)</sup><sub>max</sub>	rho<sup>(1)</sup><sub>k<sup>(1)</sup><sub>max</sub>k<sup>(1)</sup><sub>max</sub></sub>(theta<sub>K</sub>)
2	    0	     0		rho<sup>(2)</sup><sub>00</sub>(theta<sub>K</sub>)
2	    2	     0		rho<sup>(2)</sup><sub>20</sub>(theta<sub>K</sub>)
2	    2	     1		rho<sup>(2)</sup><sub>21</sub>(theta<sub>K</sub>)
2	    2	     2		rho<sup>(2)</sup><sub>22</sub>(theta<sub>K</sub>)
...
2	    k<sup>(2)</sup><sub>max</sub>    k<sup>(2)</sup><sub>max</sub>	rho<sup>(2)</sup><sub>k<sup>(2)</sup><sub>max</sub>k<sup>(2)</sup><sub>max</sub></sub>(theta<sub>K</sub>)
...
N	    0	     0		rho<sup>(N)</sup><sub>00</sub>(theta<sub>K</sub>)
N	    2	     0		rho<sup>(N)</sup><sub>20</sub>(theta<sub>K</sub>)
N	    2	     1		rho<sup>(N)</sup><sub>21</sub>(theta<sub>K</sub>)
N	    2	     2		rho<sup>(N)</sup><sub>22</sub>(theta<sub>K</sub>)
...
N	    k<sup>(N)</sup><sub>max</sub>    k<sup>(N)</sup><sub>max</sub>	rho<sup>(N)</sup><sub>k<sup>(N)</sup><sub>max</sub>k<sup>(N)</sup><sub>max</sub></sub>(theta<sub>K</sub>)
</pre>

Here theta<sub>k</sub> is the center-of-mass frame scattering angle in radians, and these must be entered smallest to largest. There is no limit on the number of theta spline points, and they do not need to be the same as in the probabilities file. The INDEX column specifies the state index, defined by the level scheme file, and these must entered in order from 1 to N. The KA and KAPPA columns specify the component of the statistical tensor, and these must be entered in odometer ordering as shown (lowest k first, then lowest kappa first). The RHOC column lists the value of the component. The largest k for state i (k<sup>(i)</sup><sub>max</sub>) is given by the lesser number of 2J<sup>(i)</sup> and 6. Note that the ground state is not included in this file.

The statistical tensor must be calculated in coordinate frame C as defined by [1]. This implies they are purely real and only have non-zero components for even k. Only positive kappa values should be inlcuded, with kappa running from 0 to k.

*Not all excited states must be included, but currently you cannot "skip" a state. If you have 5 excited states, including tensors for states 1, 2, and 3 is fine. Including tensors for states 1, 2, and 4 is not.*

*If the incoming particle trajectory is not parallel to the z-axis, which can be accomplished with optional /Beam commands, the statistical tensor will not be correct. Correcting this would require a rotation that is not yet implemented.*

References
-----------------
[1] K. Alder and A. Winter, *Electromagnetic Excitation, Theory of Coulomb Excitation with Heavy Ions*, North Holland, Amsterdam (1975).

[2] T. Czosnyka, D. Cline, and C.Y. Wu, *GOSIA User's Manual*, Bull. Am. Phys. Soc. **28**, 745 (1983). [http://www.pas.rochester.edu/~cline/Gosia/Gosia_Manual_20120510.pdf](http://www.pas.rochester.edu/~cline/Gosia/Gosia_Manual_20120510.pdf)

[3] J. Henderson, *Cygnus*, [https://github.com/jhenderson88/Cygnus](https://github.com/jhenderson88/Cygnus)
