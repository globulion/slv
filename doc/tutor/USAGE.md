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

| Object Name                 | Purpose                                                |
| --------------------------- | ------------------------------------------------------ |
|  **Preparing FRG**          | |
|  `make-frg`                 | Generates the EFP parameters for BSM solvent molecule |
|  `solefp-frg`               | Generates the SolEFP parameters for BSM solute molecule |
|  `der-dma`                  | Prepare the SolEFP parameters relating to DMTP | 
|  `der-dpol`                 | Prepare the SolEFP parameters relating to polarizabilities |
|  `der-wfn`                  | Prepare the SolEFP parameters relating to wavefunction |
|  `calc-gijk-from-hessians`  | | 
|  `convert-gijj`             | |
|  `show-frg`                 | |
|  `compare-frg`              | | 
|  `frg_check`                | |
|  `gen-camm`                 | |
|  **SolEFP Frequency shift** | |
|  `efp-xyz`                  | |
|  `efp-val`                  | |
|  `md-run`                   | |
|  `md-run_nopp`              | |
|  `efpmd-anal`               | |
|  `efpmd-traj`               | |
|  `md-prepare`               | |
|  **SolEDS Frequency shift** | |
|  `soleds-analyzer`          | |
|  `make-soleds-xyz`          | |
|  **Spectrum Simulation**    | |
|  `calc-ftir`                | |
|  `calc-tcf`                 | |
|  **Other**                  | |
|  `anal`                     | |
|  `dat`                      | |
|  `hist`                     | |
|  `slv_check`                | |
|  `dma-overall`              | |
|  `parse-itp`                | |
|  `pick-clusters`            | |










[ToDo lists]: https://github.com/globulion/slv/projects/1
