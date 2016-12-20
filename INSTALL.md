SLV Installation guide
======================

Bartosz Błasiak, Fri 14 Nov 2014, Updated: 2 Nov 2016

INTRODUCTION
------------
 
This is the installation instruction to succesfully install the Solvshift package (SLV) on UNIX systems. 
SLV has a lot of dependencies which has to be installed as for now by hand. 
The list of the dependency packages and the supported SLV version 
is enclosed here:
                                                                                                                
      
**Table 1.** To install SLV, install the packages following the order given in this Table.

| Import name       | Package Full Name                      | Recommended version |                |
| ----------------- | -------------------------------------- | ------------------- | -------------- |
| `pp`              | [Parallel Python]                      | at least 1.6.4      |    Optional    |
| `numpy`           | [Numerical Python]                     | 1.7.1               | **_Required_** |
| `scipy`           | [Scientific Python Libraries]          | 0.12.1              |    Optional    |
| `MDAnalysis`      | [Molecular Dynamics Analysis Tools]    | 0.8.1               | **_Required_** |
| `scitools`        | [Scientific Tools]                     | 0.8                 | **_Required_** |
| `PyQuante`        | [PyQuante Modified]                    | 1.0.3-BBG           | **_Required_** |
| `libbbg`          | [Quantum Chemistry Libraries package]  | at least 1.0.4      | **_Required_** |
| `coulomb`         | [Coulomb package]                      | at least 1.0.8      | **_Required_** |
| `solvshift`       | [Solvshift package]                    | at least 1.0.2      | **_Required_** |
                                                                                                          
Please follow the order of installation according to the Table 1. Particular caution has to be kept 
when installing a correct version of the packages. As for now, the versions listed 
in the Table 1 are the _tested_ ones which do not cause the API incompatibilities. Therefore, 
installing a newer version of NumPy or SciPy than recommended here could cause some problems, though 
it is not necessary at all. 

-----------------------------------------------------------------------------------------------------


INSTALLATION STEPS
------------------

### 0. `PYTHONPATH` variable.

Change `PYTHONPATH` variable if you want to install these Python modules in a non-standard directory. 
During the installation the scheme is generally like this:
                                                                                                     
```bash
python setup.py <some options may be here> install --prefix=/your/installation/path
```
                                                                                                     
Go to your `$HOME/.bashrc` file and add the following line:

```bash                                                                                                     
export PYTHONPATH=/your/installation/path/lib/pythonX.Y/site-packages:$PYTHONPATH
```

or 

```bash
export PYTHONPATH=/your/installation/path/pythonX.Y:$PYTHONPATH
```
depending on the version of Python and maybe UNIX system. 
The searchable directories can be checked by:

```python
python
>>> import sys
>>> for i in sys.path: print i
```

For this instruction the prefix will be set to `$HOME/lib64` directory. 
> Note: sometimes after the `install` step you might require logout/login step to update the changes in your system 
> (in PyQuante and NumPy cases probably)


### I. Parallel Python - `pp` module.

```bash
python setup.py install --prefix=$HOME
```

### II. Numerical Python - `numpy` module.

1. Linking with MKL Intel libraries

  1. create `site.cfg` file. Its content is given below:

     ```cfg
     [mkl]                                                          
     library_dirs = /opt/intel/mkl/lib/intel64
     include_dirs = /opt/intel/mkl/include
     mkl_libs = mkl_rt
     ```
                                                                                                                                   
  2. Modify `self.cc_exe` in `numpy/distutils/intelccompiler.py`

     ```python                                                                                                                              
     self.cc_exe = 'icc -O3 -g -fPIC -fp-model strict -fomit-frame-pointer -openmp -xhost' 
     ```
   
     In the case of ILP64 interface add also `-DMKL_ILP64` to these options. **Add other options if necessary.**
                                                                                                                                   
  3. Modify file `numpy/distutils/fcompiler/intel.py`
                                                                                                                                   
     change the line no 42 into:
     ```python
     pic_flags = ['-xhost', '-openmp', '-fp-model', 'strict', '-fPIC']
     ```
                                                                                                                                   
  4. Install by
     
     ```bash                                                                                                                              
     python setup.py config --compiler=intelem build_clib --compiler=intelem build_ext --compiler=intelem install --prefix=$HOME
     ```
                                                                                                                                         
2. Default installation.
  
   ```bash                                                                                                                                       
   python setup.py install --prefix=$HOME
   ```


### III. Scientific Python with MKL **(not necessary for Solvshift to function properly, so this step can be skipped)**

1. Scipy 0.12.1

   Install by

   ```bash
   python setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install --prefix=$HOME
   ```

2. Older Scipy versions

  1. Modify file `scipy/spatial/qhull/src/qhull_a.h` (at line 106)
     ```c++
     template <typename T>
     inline void qhullUnused(T &x) { (void)x; }
     #  define QHULL_UNUSED(x) qhullUnused(x);
     ```
     becomes:
     ```c++
     #define QHULL_UNUSED(x) (x)
     ```
         
     If you want to know more about this step, see the discussion [here](http://scipy-user.10969.n7.nabble.com/Building-SciPy-on-Debian-with-Intel-compilers-td1888.html).
                                                                                                                                                                                                     
  2. Install by

     ```bash
     python setup.py config --compiler=intelem --fcompiler=intelem build_clib --compiler=intelem --fcompiler=intelem build_ext --compiler=intelem --fcompiler=intelem install --prefix=$HOME
     ```

### IV. MDAnalysis - `MDAnalysis` module.

Install by
```bash
python setup.py install --prefix=$HOME
```

### V. Scientific Tools - `scitools` package.

Install by
```bash
python setup.py install --prefix=$HOME
```

### VI. PyQuante-Modified package.

Install by
```bash
python setup.py install --prefix=$HOME
```

### VII. Chemistry LibBBG libraries - `libbbg` module.

1. Modify the `install` script by setting the appropriate directories for prefix. 
                                                                                                                             
2. Install by
   ```bash
   bash install -p $HOME
   ```
                                                                                                                             
3. Check if you can import `libbbg` module. On certain computer clusters you may see some warnings which are OK. For example, 
   this was one of the messages:
   
   ```python
   python                                                                                                                          
   >>> import libbbg                                                                                                             
   openbabel not found in path, switching to PyQuante backend
   libint extension not found, switching to normal ERI computation
   Gtk-Message: Failed to load module "atk-bridge": libatk-bridge.so: cannot open shared object file: No such file or directory
   ```
                                                                                                                             
   Second line means that the LibInt package is not installed. 
   Third line is because incomplete installation of packages which should not
   disrupt SLV in work. However, the `$DISPLAY` variable has to be set so logging through `ssh` has to be done using `-XY` options!
   Nevertheless, ss long as there is no `ImportError` raised by Python interpreter, everything should work fine.

### VIII. Coulomb - `coulomb` package.

1. **Only if you use Python 2.6**: Modify the file `coulomb/multip.py` around line 191.
   Hash the `code for python-2.7` and unhash the `code for python-2.7`. 
   The final code should look like this:

   ```python                                                                                      
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
   ```
   >> If you use Python2.7 no changes are necessary and you can go directly to step 2 (installation).
                                                                                      
2. Install by

   ```bash
   python setup.py install --prefix=$HOME
   ```

### IX. SLV - `solvshift` module.

1. Modify the `install` script by setting the appropriate directories for prefix.
   >> Remember to set the `SLV_LIB` variable properly! As a default it is set to `$PREFIX/lib/python-2.7/dist-packages`
   >> which may be a wrong destination on your computer.

2. Install by
   ```bash
   bash install -p $HOME
   ```

3. Set the `SLV_DATA` variable and update your `PYTHONPATH` variable
   **according to the message printed on the screen**. You can do these
   by modifying your `$HOME/.bashrc` file.

4. Check if you can import `solvshift` module.

   ```python
   python
   >>> import solvshift
   ```

   It should produce neither warnings nor errors.

----------------------------------------------------------------------------------------------------

[Parallel Python]: http://www.parallelpython.com/
[Numerical Python]: http://www.numpy.org/
[Scientific Python Libraries]: https://www.scipy.org/
[Molecular Dynamics Analysis Tools]: http://www.mdanalysis.org/
[Scientific Tools]: https://pypi.python.org/pypi/SciTools
[PyQuante Modified]: https://github.com/globulion/pyq-mod
[Quantum Chemistry Libraries package]: https://github.com/globulion/libbbg
[Coulomb package]: https://github.com/globulion/clmb
[Solvshift package]: https://github.com/globulion/slv
