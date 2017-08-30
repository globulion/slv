Vibrational Solvatochromism from Molecular Dynamics
===================================================

Bartosz Błasiak, May 5, 2017  Updated: Aug 30, 2017

INTRODUCTION
------------

This section of the Mini-Tutorial describes how to use Solvshift
to analyze the molecular dynamics (MD) trajectory files. 

There are as for now two kinds of molecular systems containing
the IR probe that can be
analyzed by SLV with the aid of the MD trajectory as the source of
structures:

 - **Homogeneous systems**. This mostly refers to IR probes dissolved
   in solutions, such as MeSCN in water, chloroform or any mixture of solvents
   that produces homogeneous molecular environment around the IR reporter.

 - **Heterogeneous systems (proteins)**. In this case the IR probe is embedded
   within the polypeptide or protein framework.

HOW IT WORKS
------------

Application of SolEFP theory is always performed on a MD snapshot by reading 
the topology of the system and atomic coordinates from the trajectory. The solute IR probe
BSM and all the environmental BSM’s (and eventually electrostatic multipoles) are
deduced based on the target system topology and
superimposed according to the rigid molecule algorithm.
The solvatochromic solute-solvent interactions are quenched by using the plain cut-off scheme
with different cut-offs for Coulombic, polarizability and exchange-repulsion interactions.
At the present moment, only vibrational frequency shifts are computed.

> **Important!** Remember to *always use the post-processed MD trajectories* such that
> your IR reporter is **centered** in the periodic box in *all* MD snapshots as well
> as all the other molecules are put back into the primary cell! Otherwise,
> the results will be incorrect due to lack of sufficient layers of environment around the solute
> and migration of molecules out of primary cell
> in most of MD simulations.

The methods are implemented in the `solvshift.solefp.EFP` class, which also handles the EFP2 method.
The interface requires the information about the system’s topology and trajectory, and global settings of the computation
such as SolEFP cut-offs, solute’s normal mode etc.

## SolEFP/MD Topology conversion.

### Communication through MDInput instance

The system of interest is internally specified by MDInput object, implemented
in `solvshift.md.MDInput` class. This class reads a string from stdout or from other stream
(e.g., from a regular file) which has the following syntax:

for one fragment modelled by the
independent fragment molecule parameter object (IBM).
```sh
$Frag
  [name]     [(Sol)EFP parameter file]
  reorder    [reorder list]
  supl       [superimposition list]
  atoms      [atomic locations]            [number of EFP’s]
  ...
$endFrag
```
**Scheme 1. The syntax of SolEFP input file.** Input consists of several `$Frag...$endFrag` entries
specifying all the types of BSM in the system, their positions, superimposition and reordering lists.

First line within each `$Frag...$endFrag` section refers to molecule name and file name with parameters. 
The names of built-in BSM’s can be found by reading the documentation of `Frag` class:
```python
from solvshift.slvpar import Frag
help(Frag)
```
However, it is also possible to generate a new BSM fragment file
by using fully automatised SLV routines.

We can reorder atoms in the next line
by keyword `reorder [list]`. Also, we can provide 
superimposition list for that fragment. It can be
provided in this example either by `1 2 3 4 5`
or `1-5` (both work the same). `reorder` and `supl` are
optional. The last entry or entries refer to atomic positions in MD trajectory 
files and the number of
(Sol)EFP fragments. Note, that if frequency shift calculation is set, the first fragment
is assumed to be an IR probe molecule.

> SLV can also perform
> EFP2 interaction energy calculation and frequency shift calculation
> can be switched off.

Note also that the turn of entries within one fragment is
arbitrary, except for the first line where the fragment name and parameter file is provided.

For example, for N-methylacetamide (NMA) in water with 918 water molecules in a periodic cell
with water described by Tip4P model we could 
write the following input:

```sh
$Frag
    NMA        nma                                
    reorder    1 3 5 2 7 10 11 12 4 6 8 9
    supl       1-5
    atoms      1-12                           1
$endFrag                                          
                                                  
$Frag                                             
    water      water
    atoms      13-3684                        918
    supl       1-3
$endFrag
```

#### Passing MDInput file to the EFP class.

Once input string for `MDInput` is created correctly, one has to
initialize the SolEFP solver by invoking `EFP` instance. For example:

```python
# read the input
input = MDInput(solinp)
args, idx = input.get()

# initialize EFP object
efp = EFP(elect=True, pol=True, rep=True, disp=True, corr=True,
          ccut=35.0, pcut=16.0, ecut=13.0,
          freq=True, cunit=True, mode=4, ea=True)
```
reads the `MDInput` string (stored here as `solinp`), parses
data `args` and `idx` from `MDInput` object
and initializes SolEFP solver in `efp` object. Note that the SolEFP cut-off 
values are provided in Bohr. For more information on `EFP` class
read the documentation:
```python
from solvshift.solefp import EFP
help(EFP)
```

The next step is to set the data by
```python
# set the 
efp.set(frame[idx], *args)
efp.eval(0)
```
where `frame` is the `numpy.ndarray` object containing the atomic coordinates
of all atoms from MD topology in question. Note that SLV choses only those atoms,
which are within provided Coulomb cut-off by taking only a **slice** of `frame[idx]`.

The more elaborate example for use with one IR probe BSM for homogeneous topology
is given in the SLV tools (see `slv_md-run_nopp` utility script).


### Using topology conversion files

In the case when the molecular topology is very complex (like inside the protein)
writing input string for `MDInput` instance is tedious and impractical.
Therefore, the system of *SolEFP/MD topology conversion files* (STC) has been introduced within SLV.
These files contain general prescription on how to translate the complicated molecular
topology into the set of (Sol)EFP BSM’s. Example files are stored in `$SLV_DATA/dat/*.tc`.
However, User is urged to write his/her own STC’s.



*******
Back to [Start](https://github.com/globulion/slv/tree/master/doc/tutor/README.md)

