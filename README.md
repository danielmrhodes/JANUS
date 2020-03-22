JANUS
================================================
A GEANT4 simulation of a JANUS CoulEx experiment
================================================

Requirements
------------------

- JANUS requires GEANT4 version 10.04.02. It has not been tested on any other version.
- In order to unpack, correlate, and histogram the simulation output, a ROOT instatllation is required.

Running JANUS
-----------------

- JANUS input.mac
- correlator simulation_output.dat histogram_file.root

Functionality
-----------------
The JANUS simultation has three modes: source, scattering, and full. Each mode has different functionality and input requirements. 

- Source: Simulates a simple gamma-ray soure. No (massive) particles involved.
- Scattering: Simulates two-body scattering events. No gamma-rays involved. 
- Full: Simulates Coulomb Excitation events with user-definable level schemes and scattering-angle dependent excitation probabilities. 
