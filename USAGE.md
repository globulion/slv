Mini-Tutorial
=============

**Solvshift** is a tool for calculation of solvatochromic vibrational frequency shifts using **implicit** and **coarse-grained** models. 
The following computations can be performed at present:
********************
- **Kirkwood-Onsager-Buckingham-Cho** frequency shift calculations in *Rigid* and *Flexible* molecule approximations
- **SolX-n** vibrational solvatochromic parameters for selected mode
- **Coarse-grained SolX-n** frequency shift calculations
- **Input files** for all of these routines to create (Gaussian, Gamess or Coulomb)

********************
The *Input* mode creates input files needed for finite difference calculations using various differentiation schemes. 
The default is 5-point central Stencil method. 

To see the version information type:
```
slv
```
The command
```
slv --help
```
will result in help box with a detailed description of Solvshift functionalities and usage.

(1) Input file mode 
-------------------
## Finite field input files ##

To make cartesian displacements you need three files:
- **Geometry file** for solute molecule optimized in gas phase. 
  This file is typical `xyz` format file, e.g. for N-methylacetamide (NMA)
  this file can look like this:

```
 12

C    2.536781   -0.064900    1.271356 
N    1.220383   -0.160514    1.881594 
C    0.071609    0.197234    1.231828 
C   -1.209581    0.027796    2.030513 
O    0.061663    0.624870    0.088444 
H    3.183118    0.614997    1.833511 
H   -1.873205   -0.645760    1.484891 
H    1.149365   -0.505374    2.824985 
H   -1.053574   -0.363828    3.037984 
H   -1.708339    0.996557    2.097766 
H    3.016569   -1.045992    1.214026 
H    2.405729    0.324290    0.263094 


```
The first line contains the number of atoms, second line is blank and from the third line 
atomic coordinates are listed (units are *Angstroms*)

- **Anharmonic file** containing vibrational analysis of solute molecule 
  (having geometry identical as in the geometry file).
  The crucial variables stored in this file are as follows:
  - harmonic vibrational frequencies
  - reduced masses
  - cubic and quartic anharmonic constants
  - matrix of eigenvectors (**L** matrix transforming the derivtives from atomic Cartesian
    space to normal coordinate space and *vice versa* for transposed matrix)
  This file is technically a Gaussian `log` file for which 
  the following option is used in the anharmonic input:

```
  Freq(HPModes,Anharm,VibRot)
```
  which specifies the anharmonic analysis to be performed along with vibration-rotation
  coupling estimation. `HPModes` makes the matrix of eigenvectors printed in 5-digit high precision
  format. Beneath, examplary input specification for NMA molecule is depicted:
  
```
#p RB3LYP/6-311++G** Freq(HPModes,Anharm,VibRot,SelectAnharmonicModes)
#  scf(conver=10,xqc) iop(7/33=1) iop(3/8=2) GFInput NOSYMM 
#  integral(grid=199974) density=current
```

### **ATTENTION!!!** ###
  - In the current implementation, it is **important** to use **Cartesian** basis sets instead of 
    its spherical counterparts. It is due to the fact that calculation of distributed multipoles 
    (CAMM) by *Coulomb.py* routines requires such basis sets (there is no code for spherical Gaussians 
    as for now unfortunately). In Gaussian program the usage of Cartesian basis can be achieved either by 
    specifying `iop(3/8=2)` anywhere in task specification or `6D/10F` just after basis set querry.
  - It is customary not to use symmetry at all (`NOSYMM` keyword; otherwise the frequencies could have incorrect
    ordering for *slv* routines and **all computations could be probably WRONG if the molecule is symmetric!!!**)

--------
Having the geometry (e.g. `nma.xyz`) and anharmonic (e.g. `nma-anh.log`) files created one can make input files! 
The basic command for creating input files for Gaussian tasks and cartesian displacements of 0.005 Angstrom 
is as follows:

```
slv -s 0.005 -gi nma.xyz -a nma-anh.log
```
You should see on your screen among other informations something like this:
```
              There is NO template file!
 I am creating it on the actual directory. Please check it!

                 < gaussian.templ >

```
Now the new file `gaussian.templ` have been created. It looks like this:
```
%chk=@CHK
%mem=1800mb
%nproc=4
#p MP2/6-31G scf(conver=10,xqc) nosymm
# iop(7/33=1) GFInput freq(HPModes)
# integral(grid=199974) density=current
              
Finite Field computations

0 1
@DATA
```
The `@` delimiter specifies varying text pattern for various input files (`chk` file name or 
displaced coordinates due to finite different procedures). The rest of a template is not changing.
You can freely modify the content of template file accordingly to your needs.

If you want to create normal coordinate displacements for a selected mode create the directory for them, e.g.:
```
mkdir sder/20
cp nma.xyz ./sder/20
cd sder/20
```
`20` means the normal coordinate step. Now the command is as follows:
```
slv -m 7 -x 0 -s 20 -gi nma.xyz -a ../../nma-anh.log
```
Here, `-m` specifies normal mode number in **Python convention**. In the present example this normal mode
is thus 8th normal mode according to the list `Fundamental Bands` at the end of the anharmonic file.
The `-x` option is switched now to `0` value which indicates normal coordinates displacements.

#### `slv.step` File ####

As you perhaps noticed, after each input files creation some system file is born: `slv.step`. Its contents
is as follows:
```
DIFFERENTIATION STEP (Å):                   0.0050000000
NUMERICAL METHOD POINTS:                                    5
PACKAGE:                                                  gaussian
SECOND DERIVATIVE MODE: -1
SECOND DERIVATIVE WORKING DIRECTORY: ./sder
SECOND DERIVATIVE DIFFERENIATION STEP (DIMENSIONLESS): 10.0000000000
TRANSLATION:
```

