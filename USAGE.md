Mini-Tutorial
=============

The following small tutorial has the aim to introduce the reader into the basic use of **Solvshift** package.
If your dream is modelling the solvatochromic response of molecule immerged in solvent environment of other
molecules **this package is just for you!** Note that this software is a fresh born baby and is still under 
current intensive work and the changes and new up-to-date versions are appearing from time to time.

For the time being, **Solvshift** is generally a tool for calculation of solvatochromic vibrational frequency 
shifts using **implicit** and **coarse-grained** models. 
The following computations can be performed at present:
********************
- **SolX-n** vibrational solvatochromic parameters for selected mode
- **Coarse-grained SolX-n** frequency shift calculations from molecular dinamics (MD) trajectory (currently *Amber*)
- **Kirkwood-Onsager-Buckingham-Cho** frequency shift calculations in *Rigid* and *Flexible* molecule approximations
- **Frequency shift decomposition** using Hybrid Variation-Perturbation Interaction Energy Decomposition (EDS)
- **Input files** for all of these routines can be created ( *Gaussian*, *Gamess* or *Coulomb*)

********************
In the future the solvatochromic transition dipole change calculation, anharmonic analysis in solvents as well as 
simulation of 1D and 2D infrared spectras are also planned to be added.

Basically, **Solvshift** offers two modes of action: 
1. [Input](https://github.com/globulion/slv/blob/master/Inputs.md) and 
2. [Calculation](https://github.com/globulion/slv/blob/master/Calculations.md) mode. 
The *Input* mode creates input files needed for finite difference calculations using various differentiation schemes. 
The default is 5-point central Stencil method. The *Calculation* mode offers calculation of SolX-n parameters
and saving them on the disk and calculation of frequency shifts. Although, the current version computes also 
solvatochromic polarizability tensor and corresponding frequency shifts, this feature is under current investigation 
and will be extended to distributed solvatochromic polarizabilities which is more accurate and correct treatment.

To see the version information type:
```
slv
```
The command
```
slv --help
```
will result in help box with a detailed description of Solvshift functionalities and usage.

