 SLV Installation guide
 
 Bartosz Błasiak, Fri 14 Nov 2014

 I.   INTRODUCTION
 
 This is the installation instruction to succesfully install SLV on UNIX systems. 
 SLV has a lot of dependencies which has to be installed as for now by hand.
 
 The list of the dependency packages is enclosed here:
                                                                                                                
      
 Table 1: SLV dependencies
 | --------------- | -------------------------------------- | ------------ |
 | pp              | Parallel Python                        | 1.6.4        |
 | numpy           | Numerical Python                       | 1.7.1        |
 | scipy           | Scientific Python Libraries            | 0.12.1       |
 | MDAnalysis      | Molecular Dynamics Analysis Tools      | 0.8.1        |
 | scitools        | Scientific Tools                       | 0.8          |
 | PyQuante        | PyQuante Modified*                     | 1.0.3-BBG    |
 | libbbg          | Quantum Chemistry Libraries package    | 1.0.4        |
 | coulomb         | Coulomb package                        | 1.0.8        |
 | --------------- | -------------------------------------- | ------------ |
 | solvshift       | Solvshift package                      | 1.0.2        |
 | --------------- | -------------------------------------- | ------------ |
 * version of original PyQuante suite modified by me with a couple of extensions.
 They are: molecular multipole 1-electron integrals, overlap and kinetic 1-electron
 integral derivatives wrt nuclear coordinates, adding 6-311++G**/6D basis set for H, N, C, O, F, S and Na,
 modification of orbital shell structure for d and f shells, extension of PyQuante.Molecule and
 PyQuante.BasisSet classes, hybrid double-molecular 2-electron J and K integrals
                                                                                                           
 The order of installation is according to the Table 1. Particular caution has to be kept
 when installing a correct version of the packages. As for now, the versions listed
 in the Table 1 are the 'tested' ones which do not cause the API incompatibilities. Therefore,
 installing a newer version of NumPy or SciPy than recommended here could cause some problems, though
 it is not necessary at all. 


 II.  INSTALLATION STEPS


 II.0 $PYTHONPATH variable

      Change $PYTHONPATH variable if you want to install these Python modules in a non-standard directory.  
      During the installation the scheme is generally like this:
                                                                                                           
          python setup.py <some options may be here> install --prefix=/your/installation/path
                                                                                                           
      Go to your $HOME/.bashrc file and add the following line:
                                                                                                           
          export PYTHONPATH=/your/installation/path/lib/pythonX.Y/site-packages:$PYTHONPATH
      
      or 
                                                                                                           
          export PYTHONPATH=/your/installation/path/pythonX.Y:$PYTHONPATH
                                                                                                           
      depending on the version of Python and maybe UNIX system.
                                                                                                           
      The searchable directories can be checked by:
      
      python
      >>> import sys
      >>> for i in sys.path: print i
      
      For this instruction the prefix will be set to $HOME/lib64 directory. 

      Note: sometimes after the 'install' step you might require logout/login step to update the changes in your system
            (in PyQuante and NumPy cases probably)

 II.1 Parallel Python

      python setup.py install --prefix=$HOME


 II.2 Numerical Python

      II.2.1 Linking with MKL Intel libraries

             1) create site.cfg file:                                                                                                            
                                                                                                                                                
                [mkl]                                                          
                library_dirs = /opt/intel/mkl/lib/intel64
                include_dirs = /opt/intel/mkl/include
                mkl_libs = mkl_rt
                                                                                                                                                
             2) Modify cc_exe in numpy/distutils/intelccompiler.py
                                                                                                                                                
                self.cc_exe = 'icc -O3 -g -fPIC -fp-model strict -fomit-frame-pointer -openmp -xhost' 
                
                In the case of ILP64 interface add also -DMKL_ILP64 to these options. Add other options if needed.
                                                                                                                                                
             3) Modify file numpy/distutils/fcompiler/intel.py
                                                                                                                                                
                change the 42 line:
                   pic_flags = ['-xhost', '-openmp', '-fp-model', 'strict', '-fPIC']
                                                                                                                                                
             4) install
                                                                                                                                                
                python setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install --prefix=$HOME
                                                                                                                                         
                                                                                                                                         
      II.2.2 Linking with ATLAS
                                                                                                                                         
             python setup.py install --prefix=$HOME


 II.3 Scientific Python with MKL **(not necessary for solvshift to function properly, so this step can be skipped)**

      II.3.1 Scipy 0.12.1

         Install by

         python setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install --prefix=$HOME

      II.3.2 Older Scipy versions

      1) modify file scipy/spatial/qhull/src/qhull_a.h (line 106)                                                                                                                                    
                                                                                                                                                                                                     
         template <typename T>                                                                                                                
         inline void qhullUnused(T &x) { (void)x; }
         #  define QHULL_UNUSED(x) qhullUnused(x);
         
         becomes:
         
         #define QHULL_UNUSED(x) (x)
         
         if curious, see the discussion here: http://scipy-user.10969.n7.nabble.com/Building-SciPy-on-Debian-with-Intel-compilers-td1888.html
                                                                                                                                                                                                     
      2) install
                                                                                                                                                                                                     
         python setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install --prefix=$HOME


 II.4 MDAnalysis

      python setup.py install --prefix=$HOME


 II.5 SciTools

      python setup.py install --prefix=$HOME


 II.6 PyQuante-Modified

      python setup.py install --prefix=$HOME


 II.7 LibBBG

      1) modify install (bash script)                                                                                              
                                                                                                                                   
         Set the appropriate directories for prefix. 
                                                                                                                                   
      2) install
                                                                                                                                   
         ./install -p $HOME
                                                                                                                                   
      3) check if it works. On epsilon it gives message similar to the following:
         
         python                                                                                                                          
         >>> import libbbg                                                                                                             
         openbabel not found in path, switching to PyQuante backend
         libint extension not found, switching to normal ERI computation
         Gtk-Message: Failed to load module "atk-bridge": libatk-bridge.so: cannot open shared object file: No such file or directory
                                                                                                                                      
         Second line means that libint is not installed. Third line is because incomplete installation of packages which should not
         disrupt slv in work. However, the $DISPLAY variable has to be set so ssh logging has to be done using -XY options!


 II.8 Coulomb

      1) If Python 2.6: modify the file coulomb/multip.py around line 191:                                                                          
                                                                                            
         hash code for python-2.7 and unhash code for python-2.7. It should look like this:
                                                                                            
             ### [2] bond moments
             # code for python-2.7
             #qB  = {bond:0                            for bond in self.bonds}
             #MB  = {bond:zeros( 3     ,dtype=float64) for bond in self.bonds}
             #QB  = {bond:zeros((3,3  ),dtype=float64) for bond in self.bonds}
             #OB  = {bond:zeros((3,3,3),dtype=float64) for bond in self.bonds}
             # code for python-2.6
             qB = {}
             MB = {}
             QB = {}
             OB = {}
             for bond in self.bonds:
                 qB.update({bond:0})
                 MB.update({bond:zeros( 3     ,dtype=float64)})
                 QB.update({bond:zeros((3,3  ),dtype=float64)})
                 OB.update({bond:zeros((3,3,3),dtype=float64)})
                                                                                            
         >>>if Python2.7: no changes are necessary
                                                                                            
      2) install
                                                                                            
         python setup.py install --prefix=$HOME


 II.9 SLV

      1) modify install (bash script)                                                                                              
                                                                                                                                   
         Set the appropriate directories for prefix. 

      2) install

         ./install -p $HOME

      3) set the $SLV_DAT variable in your $HOME/.bashrc file

      4) check if it works:
     
         python
         >>> import solvshift

 The End 



