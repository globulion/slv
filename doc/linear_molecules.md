Tackling linear molecules with Solvhsift
========================================

As you perhaps noticed, there is no possibility to run the anharmonic analysis
in Gaussian for linear molecules such as thiocyanate anion. In order to create 
anharmonic file you have to create input files for calculation of Hessian, cubic
anharmonic constants as well as reduced masses and eigenvector **L** matrix.
There is a very easy tool in *Solchsift* to do this. The following steps are 
summarized below.

## 1. Input files creation

You need only the geometry file here.
