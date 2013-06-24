2. Calculation mode
===================

2.1 Coarse-grained calculations

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
  
2.1.1 Preparation of parameters

The parameters can be calculated either in flow during calculations of frequency shifts
or as a separate run. The latter method is more convenient because it enables you to 
create your own library of parameters (they are saved to a file). When the parameters 
are created in flow during frequency shift evaluation they are not saved.

To create the set of parameters for SolCHELPG method you have to type the following command:

```
slv -cg -a nma.anh --save --name nma-solchelpg-12.par
```

The parameters of **uncontracted** SolCHELPG model were saved in `nma-solchelpg-12.par` file.
