2. Calculation mode
===================

Calculation mode in *Solvshift* offers two types of calculations of frequency shifts using either 
1. goarse-grained or 2. *continuum* models. Because of the fact that the latter is of minor importance and 
its practical applicability is rather very poor we focus on **coarse-grained** calculations. 
Continuum model of solvatochromism is actually only an elegant test of solvatochromic operator
in very simple theoretical framework and should not be used in frequency shift calculations for
confronting the theory with experiment.

# 2.1 Coarse-grained calculations

The coarse-grained calculations of frequency shifts can be easily performed using *Solvshift*. 
Here it is assumed that you have created all the necessary input files and submitted these tasks
(so you have all the .log, .camm files, modified properly `slv.step` file and so on). 
The further prerequisites to be done in the following steps are:
  * solute gas-phase parameters (frequency analysis file, 
    SolX parameters, BSM moments and, in case of need, solute DMTP and dDMTP parameters)
  * files containing the information about the target system (solute and solvent positions etc.)
  
The great majority of these files is created automatically by *Solvshift* so you don't need to worry.
The creation of some of the files which have to be created by the user will be explained in the text below. 

The option which activates the Calculation mode is `-c` option. It has no additional arguments.
  
## 2.1.1 Preparation of parameters

The parameters can be calculated either in flow during calculations of frequency shifts
or as a separate run. The latter method is more convenient because it enables you to 
create your own library of parameters (they are saved to a file). When the parameters 
are created in flow during frequency shift evaluation they are not saved.

To create the set of parameters for SolCHELPG method you have to type the following command:

```
slv -c -g -a nma.anh --save --name nma-solchelpg-12.par
```
or shorter
```
slv -cg -a nma.anh -SN nma-solchelpg-12.par
```
The parameters of **uncontracted** SolCHELPG model were saved in `nma-solchelpg-12.par` file.
This file is written in *coulomb* format (see the link 
[here](https://github.com/globulion/clmb/blob/master/doc/clmb-format.md) for description). If
only `--save` option is used without specification of a name the file will be written to `slv.par`
file as a default name.

### Preparation of SolX parameters

Beneath the suitable commands are provided for creation of various SolX models for
a molecule under study:
* **SolCHELPG model**
```
slv -cg -a nma.anh -S
```
* **SolCAMM model**
```
slv -cgd -a nma.anh -S
```
* **SolMMM model**
```
slv -cgO -a nma.anh -S --transl coe
```

The last case relates to SolMMM centered at COE molecular point **(add reference!)**.
The other arguments of `--transl` can be `com` (center of mass) and others like weighted
center between two atoms. For further referene see the `--help` option.

## 2.1.2 Making contracted models

It is very easy to create an united atom (contracted) model: use `-u` option followed by a string
of numbers and dashes. For instance, our molecule is NMA (12 atoms) and we want to make two united atoms
on methyl groups. Assuming that the two carbon atoms have index 1 and 5 (first and fourth atom in the geometry
file) and the respective protons has the indexes 2,3,11 and 10,12,7, respectively, the suitable command is
as follows:

```
slv -cg -a nma.anh --save --name nma-solchelpg-12.par -u 1,2,3,11-5,10,12,7
```

The moments for protons in methyl groups are now zeroed out. We remark here, however, that
technically there is still 12 solvatochromic centers but 6 of them are **zero** and contribute nothing, so
as if they were absent.

## 2.1.3 Frequency shift calculations

To calculate the frequency shifts you have to provide the structural information
about solute-solvent system. You can prepare the necessary files directly from using
*Gaussian* log file in which population analysis were performed (the form of output
is important). The second way is to use Molecular Dynamics simulation trajectory files.
The latter way enables you to confront the results with experiment whereas the first 
route is rather for the testing of the model on the chosen level of theory. If there is 
no data about the quality of the parameters you should first do some *ab initio* or DFT
optimizations followed by vibrational analyses to correlate the exact frequency shifts
with the ones coming from chosen coarse-grained SolX-n models. Remember that SolX-n models
are purely electrostatic in first order.

### 2.1.3.1 Validation of the SolX-n models

To validate the parameters you have to make a set of various clusters of solute
with several solvent molecules (at least 30 clusters). For each system perform
frequency analysis with ChelpgG population analysis using *Gaussian*. The ChelpgG
run is necessary just for the technical purposes, mainly *slv* reads atomic coordinates
from *Gaussian* ChelpgG output. An examplary input for such model calculations is depicted
beneath:
```
#p DFT/6-311++G** iop(3/8=2) freq integral(grid=199974) scf(conver=10,xqc) pop(chelpg)
```
Of course the system must be in energy minimum.

Now, if you managed with the model frequency analysis for the solute-solvent cluster
you need to make special files: 1. `solute_<type>` and 2. `solvent_<type>`. These
files can be readily made by running the following command:
```
slv -V <N>,<type> <file>
```
For example, our solute is NMA molecule (12 atoms) and we want to create the respective
files for `A.log` and `B.log` frequency files. We simply run this:
```
slv -V 12,A A.log
slv -V 12,B B.log
```
If everything was done correctly you should see on the screen the following message:
```
 The files: 
      < solute_A > and < solvent_A > 
 have been saved in the actuall directory!
```
and the same for the second file. To calculate frequency shift basing on this files
you have to use `-f` option, which specifies frequency shift mode for model validation.
Now create the file with names of your systems, i.e. for our example the names are `A` and `B`
and the file content (let it be named `systems`) would be:
```
A
B
```
You have also know the directory with your target systems (i.e. the solute and olvent files).
Now you have them written in the actual directory but it is more convienient
to store them in special place, e.g. `target` directory. In the case of SolCAMM
and SolMMM models you need also **benchmark solvent molecule** electrostatic moments
which are multipole moments for equilibrium gas phase solvent molecule. This file
is in coulomb format and can be easily prepared using *Coulomb.py* package. 

Thus, the command running the frequency shift calculation looks like this:
```
slv -m 7 -cgfd -a nma.log -t system --read ../par/nma-solcamm6.par -b ../bsm/water.camm -D ./target
```
The new options used here are: 1. `--read` or `-R` reads the solvatochromic
parameters prepared previously, 2. `--bsm` or `-b` for BSM moments, 3. `--typ` or `-t`
specifies the file with system names and 4. `--target` or `-D` tells where are solute
and solvent files. If you are using SolChelpG parameters you don't need to use BSM as well
as you skip `-d` option.

If everything went correct you should see the final log with information on the shifts
on your screen. For example for one of the systems it looks like this:
```
 -------------------------------
          RMS analysis
 -------------------------------
  - solvent rms:   0.014602
  - solvent rms:   0.015199
  - solvent rms:   0.009917
  - solute  rms:   0.014725

 --------------------------------:--------------------------
 INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]
 --------------------------------:--------------------------
 Total               -45.28      :
 --------------------------------:--------------------------
    q-q              -21.42      :  1            -21.42
    q-D              -23.00      :  1+2          -44.42
    q-Q               -3.89      :  1+2+3        -54.42
    q-O                6.32      :  1+2+3+4      -47.74
    D-D               -6.11      :  1+2+3+4+5    -45.23
    D-Q                0.36      :
    D-O                0.58      :
    Q-Q                1.93      :
    Q-O               -5.09      :
    O-O                5.04      :
 --------------------------------:--------------------------
```
The first window tells you the RMS (in atomic units) from superimpositions of solute and solvent molecules
with gas-phase solute and BSM. The RMS values should be smaller than 0.05 au. If the RMS
is larger than 0.5 probably you have done incorrect files (check the atom ordering in those cases!).
The next window shows the frequency shift analysis with providing the asymptotic
convergence (right panel) and the magnitudes of each separate term (left panel).
If you don't use SolMMM model **remember to ignore** `1+2+3+4+5` **value because hexadecapoles are not implemented!**
This convergence value is valid only for SolMMM analysis when solute and solvent molecules
are **neutral** (terms with hexatecapoles are then zero). Thus for SolCAMM models
the most accurate prediction of frequency shift is `1+2+3+4` (up to R-4 terms included).

### 2.1.3.2 Calculations of frequency shifts from MD trajectories

In preparation. As for now you can calculate frequency shifts
from:
* *Gromacs* using `.xtc` trajectories and `.itp` file with charges
* *Amber* using `.mdcrd` trajectories and `.prmtop` file with charges

The command is as follows:
```
slv -m 7 -cfd -a nma.log --read nma-solcamm6.par --md-package amber -M traj.dcd topology.prmtop
```
The option `--md-package` or `U` provides the package and `-M` provides trajectory. Charges
are treated as last argument so it should be nothigh after their specification within a command line.
In this run don't use `-f` option because it is not validation of the model mode.

## 2.1.4 Corrections to the SolX-n frequency shifts

The correction terms to frequency shifts coming from solely distributed 
solvatochromic moments can be easily evaluated in Solvshift. However, you have to prepare
semi-manually the next two files containing the distributed multipole moments and their 
derivatives with respect
to the normal coordinate of interest. To do this first:
* calculate the moments.
If you were creating the solvatochromic parameters **you already have this file**. Otherwise
create it using *Coulomb.py* package.
Next 
* calculate the first derivatives of moments by typing:
```
slv -m 7 -cg -a nma.log --print > derivatives.txt
```
The derivatives and other very useful information are printed to `derivatives.txt` file.
Then select the first derivatives of your mode of interest and copy them **maintaining coulomb format**
to the new file, let say, `nma-dcamm12.par`. Now copy the structure of the solute molecule
from parameter file (`Structure` section) to the derivative file (under the derivatives).
 
It is important to note here, that your  electrostatic moments as well as their derivatives 
are fully distributed (no united atoms). However, suppose we want to use SolCAMM-6 model
so we should contract the derivatives as well as electrostatic moments in the same way.
Fortunately, it is very easy to contract the models as you wish. You just have to have LIBBBG
installed and everything will came very strightforward. Do it e.g. within python:
```
python
>>> from utilities import ParseDMA
>>> ualist = [ (3,5,6,7), (8,9,10,11) ]
>>> dma = ParseDMA('nma-dcamm12.par','coulomb')[0]
>>> dma.MakeUa(ualist)
>>> dma.write('nma-dcamm6.par')
```
I think that these commands are self-explanatory. You can make a script which contracts
for you any distributed multipole moment set.

When you have these two necessary additional files type the command like this:
```
slv -m 7 -cgfd -a nma.log -t system -R ../par/nma-solcamm6.par -b ../bsm/water.camm -D ./target \
    -Xz nma-camm-6.par -Z nma-dcamm6.par
```
The output is similar as after normal frequency shift calculation run now but the next two  windows
are present at the bottom:
```
 -------------------------------
          RMS analysis
 -------------------------------
  - solvent rms:   0.014602
  - solvent rms:   0.015199
  - solvent rms:   0.009917
  - solute  rms:   0.014725

 --------------------------------:--------------------------
 INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]
 --------------------------------:--------------------------
 Total               -45.28      :
 --------------------------------:--------------------------
    q-q              -21.42      :  1            -21.42
    q-D              -23.00      :  1+2          -44.42
    q-Q               -3.89      :  1+2+3        -54.42
    q-O                6.32      :  1+2+3+4      -47.74
    D-D               -6.11      :  1+2+3+4+5    -45.23
    D-Q                0.36      :
    D-O                0.58      :
    Q-Q                1.93      :
    Q-O               -5.09      :
    O-O                5.04      :
 --------------------------------:--------------------------
 CORRECTION TERMS [cm-1]         : CORRECTED SHIFTS [cm-1]  
 --------------------------------:--------------------------
           MA         EA         :  1            -21.42
 R-2      -5.48      -1.16       :  1+2          -51.06
 R-3      -6.80       3.79       :  1+2+3        -64.07
 R-4      -4.83       1.68       :  1+2+3+4      -60.55
                                 :  1+2+3+4+5       ???
 --------------------------------:--------------------------
```
The left panel contains the information about the corrections to 
mechanical (`MA`) and electonic (`EA`) anharmonicities, respectively,
providing each of R-n contributions. The right panel contains corrected
frequency shifts.

# 2.2 Continuum models of solvatochromism
