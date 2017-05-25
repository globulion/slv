Mini-Tutorial
=============

Bartosz Błasiak, Summer 2013, Updated: 25 April 2017

INTRODUCTION
------------

The following small tutorial has the aim to introduce the reader into the basic use of **Solvshift** package.
If your dream is modelling the solvatochromic response of molecule immerged in solvent environment of other
molecules **this package is just for you!** Note that this software is a relatively fresh born tool and is still under 
current intensive work. Therefore, it has not been published.

For the time being, **Solvshift** is generally a tool for calculation of the solvatochromic vibrational frequency 
shifts using mostly **coarse-grained** [models](https://github.com/globulion/slv/blob/master/README.md). 
In the future, the implementation of the solvatochromic interaction-induced transition dipole change,
simulation of 1D and 2D infrared spectra and general improvement of accuracy and performance 
are planned. See the [ToDo lists] for reference upon current development tasks.

INTERACTING WITH Solvshift
------------------------------

Solvshift can be used in many ways depending on what task is to be performed.


### Interacting through main executable.

In general, simple calculations and tasks can be specified by using the 
**Solvshift** main executable script, `slv`. 
To see the version information run:
```
slv
```
The command
```
slv --help
```
will result in help box with a detailed description of Solvshift functionalities and usage.
> Many of the options are now obsolete. This will be fixed in the future.

### Interacting through SLV tools.

Most of the basic steps in setting the calculations can be done by using utility scripts
that are stored in `util` directory. Below I list their names and function.

*******
**Table 1. SLV utility scripts.** The overall tool name is constructed from `slv_[Object Name]`.		

| Object Name                 | Purpose                                                     | Comment   |
| --------------------------- | ----------------------------------------------------------- | --------- |
|                             | **SolEDS Frequency Shifts**                                 |           |
|  `soleds-analyzer`          | SolEDS inputs and calculations                              | *Important* |
|  `make-soleds-xyz`          | Construct constrained cluster for SolEDS analysis           |           |
|                             | **SolEFP Frequency Shift***                                 |           |
|  `efp-xyz`                  | SolEFP computation on xyz structure                         | Useful    |
|  `efp-val`                  | SolEFP computations on set of xyz structures                |           |
|  `md-run`                   | Runs SolEFP/MD calculations with parallel implementation    | Currently not working well |
|  `md-run_nopp`              | Runs SolEFP/MD calculations without parallel mode           | *Important* |
|  `efpmd-anal`               | EFPMD trajectory analysis                                   | *Important* |
|  `efpmd-traj`               | EFPMD trajectory post-processing                            | *Important* |
|  `md-prepare`               | Generate input file for SolDMTP/MD calculations             | Useful    |
|                             | **Preparing FRG Files**                                     |           |
|  `make-frg`                 | Generates the EFP parameters for BSM solvent molecule       | *Important* |
|  `solefp-frg`               | Generates the SolEFP parameters for BSM solute molecule     | *Important* |
|  `der-dma`                  | Prepare the SolEFP parameters relating to DMTP              |           | 
|  `der-dpol`                 | Prepare the SolEFP parameters relating to polarizabilities  |           |
|  `der-wfn`                  | Prepare the SolEFP parameters relating to wavefunction      |           |
|  `calc-gijk-from-hessians`  | Compute cubic anharmonic constants from Hessian matrices    |           | 
|  `convert-gijj`             | Convert cubic anharmonic constants units                    |           |
|  `show-frg`                 | Display fragment info                                       |           |
|  `compare-frg`              | Compare two fragments                                       |           |
|  `frg_check`                | Check the fragment file                                     |           |
|  `gen-camm`                 | Generate CAMM distribution                                  | *Important* |
|                             | **Simulations of Spectra**                                  |           |
|  `calc-ftir`                | Compute FTIR spectrum                                       | *Important* |
|  `calc-tcf`                 | Calculate time correlation functions and response functions | *Important* |
|                             | **Other**                                                   |           |
|  `anal`                     | Calculates averages from SolEFP/MD output                   | Useful    |
|  `dat`                      | | |
|  `hist`                     | Plot histograms from SolEFP/MD output                       |           |
|  `check`                    | Compute derivatives of overlap integrals numerically and analytically | |
|  `dma-overall`              | Compute singe-centered multipoles from DMTP distribution    | Useful    |
|  `parse-itp`                | Parse the charges from GROMACS itp file                     |           |
|  `pick-clusters`            | Pick the solute-solvent clusters from MD trajectory file    | Useful    |

*******

Some of the above tools automatize work with SLV by the use of other tools. Therefore, not all 
of the tools are important for a standard use and most of them are for debugging or developing 
Solvshift. To see the usage of these tools run them without arguments or with `-h` option.
> The option parsing scheme is now implemented for only some of the tools whereas the rest
> rely on simple argument list from command line. This will be changed in the future to 
> unify the way of using the tools, including printing the usage information (with `-h` option).

### Interacting through Python level.

Solvshift and its modules can be accessed diretcly from Python level. 
This is particularly useful for writing custom analysis scripts and is the most 
advanced form of using SLV code. When starting writing own analysis script one can 
see the SLV tools as examples. Below, the most important classes and modules
are listed. The exemplary use of them is also presented in the tutorial [ADD].

*******
**Table 2. The most useful classes of Solvshift and its modules.** Classes are denoted by `class`
           whereas functions by `def`.

| Object Name                 | Module                                                      | Short Description   |
| --------------------------- | ----------------------------------------------------------- | --------- |
| `class Frag`                | `solvshift.slvpar` | EFP2 and SolEFP fragment |                        
| `class EFP`                 | `solvshift.solefp` | EFP2 and SolEFP calculation method |
| `class MDInput`             | `solvshift.md`     | SolEFP/MD input file | 
| `class UNITS`               | `libbbg.units`     | Units and standard converters | 
| `class DMA`                 | `libbbg.dma`       | Distributed multipole analysis object | 
| `class QMFile`              | `libbbg.utilities` | File with quantum chemistry calculation results |
| `class PotentialContourMap` | `libbbg.utilities` | Contour plots for electrostatic potential |
| `class RungeKutta`          | `libbbg.utilities` | Runge Kutta numerical integration method | 
| `class Allign`              | `libbbg.utilities` | Align the molecule along specified axis |
| `class Peak`                | `libbbg.utilities` | Spectral peak |
| `class Peak2DIR`            | `libbbg.utilities` | 2DIR peak |
| `class QMOscillator`        | `libbbg.utilities` | Quantum 1-dimensional oscillator in a specified potential |
| `class CUBFLE`              | `coulomb.cube`     | Gaussian CUBE file |
| `class QMMap`               | `coulomb.cube`     | 3D isocontour object of a scalar function                 |
| `class MULTIP`              | `coulomb.multip`   | Distributed multipole expansion method                    |
| `def get_elmtp`             | `libbbg.utilities` | Compute interaction energy between two DMTP objects |
*******






[ToDo lists]: https://github.com/globulion/slv/projects/1
