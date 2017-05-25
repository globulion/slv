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
whereas the command line used to compute the SolEDS frequency shifts is
```
slv_soleds-analyzer 2 [anh] [eds dirs] [mode] [method]
```
The meanings of the arguments passed to `slv_soleds-analyzer` are explained in Table 1.

*******
**Table 1. Meaning of the arguments passed to the `slv_soleds-analyzer` tool.** 

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

#### Cubic anharmonic constant sign problem.

After performing step 1 (constrained cluster preparation and all the input files)
one has to manually inspect the signs of the cubic anharmonic constants in the 
Gaussian anharmonic analysis file **before the SolEDS calculations are performed!** 
This is necessary because the phase of the eigenvectors
from the Gaussian harmonic analysis is not always in agreement with the phase of eigenvectors
used in the anharmonic code of Gaussian. Therefore, certain cubic anharmonic constants 
may have signs which are incopmatible with the Wilson matrix printed (and read by SLV)
in the Gaussian anharmonic file. Therefore, Solvshift has another utility that creates 
separate input files for calculations of the cubic anharmonic constants. This utility
is automatically called when running soleds-analyzer. What one needs to do is to

Case study: MeNC-H2O cluster.
-----------------------------

Here I describe the procedure to obtain the SolEDS frequency shifts
for methylisocyanate interacting with one water molecule.

After running the soleds-analyzer command to create input files
```
slv_soleds-analyzer 1 menc-water.xyz menc.xyz menc-freq.log 4 17 1-4 3 0.006 15
```
Solvshift will ask for the GAMESS template input file `gamess.templ` if such file
is not present in the working directory.
The contents of this file is shown below:

******
**Scheme 1. GAMESS template file for SolEDS calculation.** The name of the file is `gamess.templ` 
            (name should should not be changed).
```
$system mwords=7 memddi=12 parall=.t. $end
$contrl scftyp=rhf runtyp=eds icharg=0 mult=1 units=angs
        maxit=100 exetyp=run ispher=-1
        mplevl=2 cctyp=none $end
$eds    mch(1)=0,0 mmul(1)=1,1 mnr(1)=7 ffeds=.f. $end
$scf    dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
        conv=1d-10 $end
$basis  gbasis=slvbas extfil=.t. $end
$data
--- CALCULATIONS OF SOLVATOCHROMIC FREQUENCY SHIFTS OF MENC-WATER FROM EDS ---
C1 0
@DATA
$end
```
******

After performing all the calculations (submitting each Gaussian and GAMESS input file)
one can compute the frequency shifts by the following command:
```
slv_soleds-analyzer 2 menc-freq.log ./:./sder/15 4 HF
```
which produces the following output:
```
E(EL,10)            5.97      3.02      8.98 
E(EL,MTP)           7.38      1.47      8.85
E(EL,PEN)          -1.42      1.55      0.13
E(EX,HL)           13.03      1.24     14.27
E(DEL)             -2.36     -0.41     -2.77
DE(HF)             16.63      3.84     20.48
```
The first and second columns are the mechanical and electric anharmonicity,
whereas the third column is their sum, respectively. In this example
the total Hartree-Fock vibrational frequency shift is +20.48 cm^-1.
Coulombic, exchange-repulsion and charge delocalization frequency shifts are
8.98, 14.27 and -2.77 cm-1, respectively.

*******
Back to [Start](https://github.com/globulion/slv/tree/master/doc/tutor/README.md)

