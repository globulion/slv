2. Calculation mode
===================

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
slv -cg -a nma.anh --save --name nma-solchelpg-12.par
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

### 2.1.1.1 Preparation of SolX parameters

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
slv -cgO -a nma.anh -S --transl coe [can be com and others, see the --help option)
```

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
