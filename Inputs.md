1. Input file mode 
==================

## Finite field input files ##

To make cartesian displacements you need two files:
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
  - matrix of eigenvectors ( **L** matrix transforming the derivtives from atomic Cartesian
    space to normal coordinate space and *vice versa* for transposed matrix).

  This file is technically a Gaussian `log` file for which 
  the following option is used in the anharmonic input:`Freq(HPModes,Anharm,VibRot)`
  which specifies the anharmonic analysis to be performed along with vibration-rotation
  coupling estimation. `HPModes` makes the matrix of eigenvectors printed in 5-digit high precision
  format. Beneath, examplary input specification for NMA molecule is depicted:
  
```
#p RB3LYP/6-311++G** Freq(HPModes,Anharm,VibRot,SelectAnharmonicModes)
#  scf(conver=10,xqc) iop(7/33=1) iop(3/8=2) GFInput NOSYMM 
#  integral(grid=199974) density=current
```

### **ATTENTION!!!** ###
  - It is **indispensable** not to use symmetry at all (`NOSYMM` keyword), otherwise the frequencies could 
    have incorrect ordering for *slv* routines and 
    **all computations could be probably WRONG if the molecule is symmetric!!!**
  - In the current implementation, it is **necessary** to use **Cartesian** basis sets instead of 
    its spherical counterparts. It is due to the fact that calculation of distributed multipoles 
    (CAMM) by *Coulomb.py* routines requires such basis sets (there is no code for spherical Gaussians 
    as for now unfortunately). In *Gaussian* program the usage of Cartesian basis can be achieved either by 
    specifying `iop(3/8=2)` anywhere in task specification or `6D/10F` just after basis set querry.

--------
Having the files with geometry (e.g. `nma.xyz`) and anharmonic analysis (e.g. `nma-anh.log`) created 
one can make input files! 
The basic command for creating input files for *Gaussian* tasks and cartesian displacements of 0.005 Angstrom 
is as follows:

```
slv -s 0.005 -g -i nma.xyz -a nma-anh.log
```
or shorter
```
slv -s 0.005 -gi nma.xyz -a nma-anh.log
```
The option `-g` indicates *Gaussian* tasks to be run. If it is ommited, the corresponding *Gamess* 
input mode is switched on as a default. Try to run one of the above commands. If everything is OK 
you should see on your screen among other informations something like this:
```
              There is NO template file!
 I am creating it on the actual directory. Please check it!

                 < gaussian.templ >

```
Note, that the new file `gaussian.templ` have been created which is necessary to create input file set. 
The content of this file looks like this:
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
The `@` delimiter specifies varying text pattern for various input files: 
- `@CHK` for `chk` file names (*Gaussian* checkpoint files; relevant if *Gaussian* input mode is choosen)
- `@DATA` for displaced coordinates due to finite different procedures. 
- `@FRAG` for specifying fragments (optional; see later)

The rest of a template is not changing.
You can freely modify the content of template file accordingly to your needs. In the case of *Gamess* input mode
the `gamess.templ` files is created instead.

If you want to create normal coordinate displacements for a selected mode **create the directory for them**, e.g.:
```
mkdir sder/20
cp nma.xyz ./sder/20
cd sder/20
```
`20` means the normal coordinate step (just a name of directory, it has no meaning for `slv`). 
Now the command is as follows:
```
slv -m 7 -x 0 -s 20 -gi nma.xyz -a ../../nma-anh.log
```
Here, `-m` specifies normal mode number in **Python convention** (N+1). In the present example this normal mode
is thus 8th normal mode according to the list `Fundamental Bands` at the end of the anharmonic file.
The `-x` option is switched now to `0` value which indicates normal coordinate displacements.

#### `slv.step` File ####

As you perhaps noticed, after each input file set creation some system file is born: `slv.step`. Its contents
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
The line `SECOND DERIVATIVE MODE:` needs some comment. If the mode is negative then `slv` does not 
take into account **electronic anharmonicity** for parameter calculations. Thus, you should  manually change this line
for your purpose to include this important contribution. For the presented example the line should look like this:
```
SECOND DERIVATIVE MODE: 7
```
Secondly, you should type actual second derivative directory, so in our case we would have to change the 
next line into:
```
SECOND DERIVATIVE WORKING DIRECTORY: ./sder/20/
```
Finally, the correct step for normal mode differentiation should be retyped too:
```
SECOND DERIVATIVE DIFFERENIATION STEP (DIMENSIONLESS):   20
```
#### Submiting the tasks ####

After you succeded in creation of inputs for cartesian displacements as well as normal coordinate displacements
and also you changed appropriately `slv.step` file you can submit all of the input files to some computer
cluster where the desired quantum chemistry program is installed. Because the number of input files is generally
quite huge (e.g. for NMA 145 inputs is required for first derivatives and 5  additional inputs are necessary for 
second derivatives wrt selected mode) the queueing system is very onvenient because you can submit all of the inputs
to the queue. 
