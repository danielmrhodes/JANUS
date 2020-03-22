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
The JANUS simultation has three modes: source, scattering, and full. Each mode has different functionality and input requirements. 

- Source: Simulates a simple gamma-ray soure. No (massive) particles involved.
- Scattering: Simulates two-body scattering events. No gamma-rays involved. 
- Full: Simulates Coulomb Excitation events with user-definable level schemes and scattering-angle dependent excitation probabilities.
