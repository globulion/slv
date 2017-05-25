Vibrational Solvatochromism Based on Interaction Energy Decomposition
=====================================================================

Bartosz Błasiak, May 5, 2017  Updated: 

INTRODUCTION
------------

This section of the Mini-Tutorial describes the usage of SolEDS method
within Solvshift. SolEDS
aims in calculation of the vibrational frequency shifts in small systems
(up to few tens of atoms) which are in energy minimum. The method
requires constructing the, so called, **constrained cluster** which is a
model of a true energy minimum such that the internal coordinates within 
the solute molecule are identical to the ones in the gas phase. The constrained
cluster has to imitate the true energy minimum as close as possible.

Figure 1 shows the main principle of SolEDS. The gas phase solute structure
is depicted as 1, whereas the true energy minimum structure (fully energy optimized)
is depicted as 2. The constrained cluster is shown as number 3 in this figure.

******
**Figure 1. The principle of SolEDS.** 
<img src="soleds-scheme.png" alt="Drawing" style="width: 20px;"/>
******

### Constrained Cluster.

Solvshift has an automatic tool that creates the constrained cluster
based on the true energy minumum geometry and solute gas phase structure.
However, it is always necessary to check if the model is sufficiently
close in geometry to the true energy minimum. Even seemingly small 
deviations could lead to large errors when computing the frequency shifts.
Therefore, it is recommended to always check the model constructed automatically
by SLV, and, if necessary, upgrade it manually. 

>Remember that *the SolEDS 
>results depend highly on the quality of the constrained cluster!*

Below, I list the features of a "good" constrained cluster:
 1. The internal coordinates of solvent molecules with respect to the
    solute atoms involved in the normal mode of interest are almost identical
    to the true energy minimum
 2. All the other solvent molecules that are not in the closest contact
    with the solute normal mode atoms are placed around solute in almost exactly 
    same way as in the true energy minimum.

### Performing SolEDS calculation.

The SolEDS calculations are automatised in the SLV utility tool `slv_soleds-analyzer`.
First one needs to prepare GAMESS input files of the finite difference (FD) displaced structures
derived from the constrained cluster. These FD's are necessary to compute the first and second
derivatives of the interaction energy components numerically with respect to the solute's 
normal coordinates.

The basic syntax for generation of the constrained cluster and all the input files
is as follows:
```
slv_soleds-analyzer 1 [full opt xyz] [gas phase xyz] [anh] [mode] [line] [suplist] <n> <step-cart> <step-mode>
```
wherease the command line used to compute the SolEDS frequency shifts is
```
slv_soleds-analyzer 2 [anh] [eds dirs] [mode] [method]
```
The meanings of the arguments passed to `slv_soleds-analyzer` are explained in Table 1.

*******
**Table 1. Meaning of arguments passed to `slv_soleds-analyzer` tool.** 

| Argument        | Explanation                                                                          | 
| --------------- | ------------------------------------------------------------------------------------ | 
| `full opt xyz`  | xyz file with fully optimized solute-solvent cluster                                 |
| `gas phase xyz` | xyz file with solute in gas-phase (anharmonic file orientation)                      |
| `anh`           | g09 anharmonic file with CUBIC (remember to change the wrong signs!)                 |
| `mode`          | normal mode number (in Helico, normal numbers - not in Python convention)            |
| `line`          | line number in GAMESS input file after which solvent molecule is to be inserted      |
| `n`             | the number of solvent atoms (default 3 - one water molecule)                         |
| `step-cart`     | differentiation step in cartesian coordinates (default 0.006 Angs)                   |
| `step-mode`     | differentiation step in normal coordinates (default 15.0)                            |
| `eds dirs`      | For example: w20:sder/20 - parse w20 as fder_dir and w20/sder/20 as sder_dir         |
| `method`        | HF or MP2, depending on the type of EDS calculations.                                |

*******


Back to [Start](https://github.com/globulion/slv/tree/master/doc/tutor/README.md)

