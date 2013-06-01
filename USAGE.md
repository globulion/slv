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
- **Geometry file** for solute molecule optimized in gas phase. This file is typical `xyz` format file, e.g. for N-methylacetamide
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
The first line contains the number of atoms, second line is blank and from the third line atomic coordinates are listed (units are *Angstroms*)

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
  coupling estimation. `HPModes` makes the matrix of eigenvectors printet in 5-digit hihg precision
  format. Beneath, examplary input specification for NMA molecule is depicted:
  
```
#p RB3LYP/6-311++G** Freq(HPModes,Anharm,VibRot,SelectAnharmonicModes)
  scf(conver=10,xqc) iop(7/33=1) iop(3/8=2) GFInput NOSYMM integral
 (grid=199974) density=current
```

  ### **ATTENTION** ###
  - In the current implementation, it is **important** to use **Cartesian** basis sets instead of 
    its spherical counterparts. This is achieved either by specifying `iop(3/8=2)` or `6D/10F`
    just after basis set querry.
  - It is customary not to use symmetry at all (otherwise the frequencies could have incorrect
    ordering for *slv* routines and **all computations will be WRONG!!!**)
  
