#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------- #
#             EFFECTIVE FRAGMENT PARAMETER FORMAT MODULE                #
# --------------------------------------------------------------------- #

from numpy     import array, float64, zeros, newaxis, sqrt,  \
                      dot, asfortranarray, transpose, arange,\
                      concatenate as con
from units     import *
from dma       import DMA
from utilities import order, reorder, SVDSuperimposer as svd_sup, MakeMol
from efprot    import vecrot, vc1rot, tracls, rotdma, tracl1, rotdm1
from PyQuante.Ints import getbasis
import sys, copy, os, re, math
sys.stdout.flush()

__all__     = ['Frag',]
__version__ = '2.0.1'

class Frag(object,UNITS):
    """
 =============================================================================
                 Effective Fragment Parameter Format System                   
                                                                              
            Bartosz Błasiak                last revision: 9 Nov 2013          
 =============================================================================
                                                                              
 Represents Efective Fragment Potential (EFP) fragment. Designed to unify func
 tionality of the QM molecular object that contains information about all inte
 rmolecular interactions based on gas-phase molecular properties, mainly first
 order Coulomb electrostatics, second-order polarization forces, first-order c
 orrection due to overlap and Pauli repulsion, charge-transfer and dispersion.
 Supports the following functionalities for condensed phase physics modeling: 
 
   a) intermolecular interaction energies                                     
   b) nonlinear electrooptical properties                                     
   c) infra-red solvatochromism                                               
 
 It contains various memory sections necessary for evaluation of molecular agg
 regate property using coarse-grained models. The major groups of molecular pa
 rameters are:                                                                
 
   I  ) Basic molecular specification                                         
        a) atomic coordinates, numbers and masses                             
 
   II ) Electrostatic population parameters                                   
        a) distributed multipole moments                                      
        b) distributed dipole-dipole polarizabilities                         
 
   III) Nonlinear electrooptical properties                                   
        a) molecular charge and dipole moment                                 
        b) molecular polarizability                                           
        c) molecular first- and second-order hyperpolarizabilities            
 
   IV ) Frequency analysis data                                               
        a) harmonic frequencies, reduced masses and mass-weighted eigenvectors
        b) cubic anharmonic force constants                                   
 
   V  ) Wave function and QM energetics (CMO and LMO space)                   
        a) wave-function coefficients                                         
        b) LMO centroids                                                      
        c) Fock matrix                                                        
        d) derivatives of the above wrt electric field and normal coordinates 
 
 -----------------------------------------------------------------------------
 All of the above mentioned quantities are stored in ATOMIC UNITS!!!          
 -----------------------------------------------------------------------------
 
 Beneath, I list the further contents of this documentation box:              
    * How to use Frag class                                                   
    * Listing of memorials stored in the Frag instance                        
    * Information about storage dimensions of Frag memorials                  
 
 =============================================================================
  [1] USAGE (if '=': the argument is default)                     
 =============================================================================
 
  1) create the object (initialize)                                           
                                                                             
     frag = Frag()       initialize empty object or ...
     frag = Frag(name)   where name is the name of molecule or fragment writte
                         n in SLV standard fragment library $SLV_PATH/frg
                         The following are possible:
    ----------------------------------------------------------------------------                   
       name       molecule                EINT   FREQ   NLO      comment                           
    ----------------------------------------------------------------------------                   
     o water      H2O                     YES    NO     NO       with DMA-5                        
     o water2     H2O                     YES    NO     NO       with CAMM-3                       
     o dmso       DMSO                    YES    NO     NO                                         
     o meoh       Methanol                YES    NO     NO                                         
     o chcl3      Chloroform              YES    NO     NO                                         
     o dcm        Dichloromethane         YES    NO     NO  
     o na+        Sodium cation           YES    NO     NO
     o me-so3-    Methyl sulphonate (-1)  YES    NO     NO
     o so3--      SO3 (2-) anion          YES    NO     NO
     o meoac      Methyl acetate          YES    YES    NO       C=O stretch                       
     o nma        N-methyl acetamide      YES    YES    NO       amide I mode                      
     o nma-d7     NMA-D7                  YES    YES    NO       amide I' mode
                                                                (all H deuterated)
     o mescn      Methyl thiocyanide      YES    YES             C=N stretch
                                                                            
    ----------------------------------------------------------------------------                   
     frag = Frag(file)   where file is the Frag format file (*.frg) with parame
                         ters created by a user. See slv_der-* for fragment par
                         ameter generators
     frag(file)          after initialization of Frag instance read the *.frg file
                                                                              
  2) set the data to frag from other source than Frag format file (*.frg)     
                                                                              
     frag.set(mol=None, anh=None, frag=None, dma=None, chelpg=None, esp=None) 
                         mol  - Molecule    class object                      
                         anh  - FREQ        class object                      
                         frag - FragFactory class object                      
                         dma  - DMA         class object                      
                         and other:                                           
                         * elpg - array of ChelpG charges                     
                         * esp    - array of ESP    charges                   
                         * other keywors that match list of memorials (see [2])
                                                                              
  3) write the parameters to a file                                           
                                                                              
     frag.write(file='slv.frg')                                               
                                                                              
  4) actions on Frag object                                                   
                                                                              
     frag.translate(transl)     # translate all tensors by cartesian vector   
     frag.rotate(rot)           # rotate all tensors by unitary 3x3 matrix    
     frag.sup(xyz)              # superimpose (translate + rotate) to the pos
                                # indicated by xyz Cartesian coordinates
     frag.reorder(ord)          # reorder the turn of atoms by specifying ord
                                # vector, e.g.: ord=[3,1,2] is equivalent to
                                #        A-B-C ----> B-C-A
                                #        1 2 3
     f_copy = frag.copy()       # create a deep copy of frag objects
     print frag                 # print the status of frag object

  5) return objects from frag object

     par = frag.get()           # get parameter dictionary. The keys are the
                                # 'Shortcut' names in section [1] of this docu
                                # mentation, e.g.:
     vecl1 = par['vecl1']       # gets LMO wave function first derivatives wrt
                                # all normal coordinates
     name= frag.get_name()      # return the name of fragment 
     pos = frag.get_pos()       # get position Cartesian coordinates of atoms
     bfs = frag.get_bfs()       # get BasisSet object                        
     qad, oct = frag.get_traceless() # get traceless quadrupoles and octupoles
                                                                              
 =============================================================================
  [2] LISTING OF MEMORIAL NAMES                                               
      The *.frg file has two parts:                                           
      * preambule (molecule section)                                          
      * parameters sections                                                   
 =============================================================================
  GR : Section name                      Shortcut       Dimension             
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : molecule                    · · · mol                                  
     :                             ·      name                                
     :                             P      shortname                           
     :                             R      basis                               
     :                             E      nbasis                              
     :                             A      method                              
     :                             M      atoms                               
     :                             B      natoms                              
  I  :                             U      ndma                                
     :                             L      npol                                
     :                             E      nmos                                
     :                             ·      ncmos                               
     :                             · · ·  nmodes                              
     : Atomic numbers                    atno            natoms               
     : Atomic masses                     atms            natoms               
     : Atomic coordinates                pos             natoms, 3            
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : ESP charges                       esp             natoms               
     : ChelpG charges                    chlpg           natoms               
     : DMTP centers                      rdma            ndma  , 3            
     : DMTP charges                      dmac            ndma                 
     : DMTP dipoles                      dmad            ndma  , 3            
 II  : DMTP quadrupoles                  dmaq            ndma  , 6            
     : DMTP octupoles                    dmao            ndma  , 10           
     : DMTP charges fder                 dmac1           nmodes, ndma         
     : DMTP dipoles fder                 dmad1           nmodes, ndma ,3      
     : DMTP quadrupoles fder             dmaq1           nmodes, ndma ,6      
     : DMTP octupoles fder               dmao1           nmodes, ndma ,10     
     : DMTP charges sder                 dmac2           ndma                 
     : DMTP dipoles sder                 dmad2           ndma ,3              
     : DMTP quadrupoles sder             dmad2           ndma ,6              
     : DMTP octupoles sder               dmad2           ndma ,10             
     : Polarizable centers               rpol            npol  , 3            
     : Distributed polarizabilities      dpol            npol  , 9            
     : Distributed polarizabilities      dpol1           nmodes, npol, 9      
     : - first derivatives                                                    
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : Harmonic frequencies              freq            nmodes               
     : Reduced masses                    redmass         nmodes               
 III : Mass-weighted eigenvectors        lvec            nmodes, natoms, 3    
     : Cubic anharmonic constants        gijk            nmodes, nmodes, nmodes
 - - : - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     : LMO centroids                     lmoc            nmos  , 3            
     : LMO centroids                     lmoc1           nmodes, nmos  , 3    
     : - first derivatives                                                    
     : Fock matrix                       fock            nmos  , nmos         
     : Fock matrix                       fock1           nmodes, nmos  , nmos 
     : - first derivatives                                                    
     : Canonical Fock matrix             fckc            ncmos , ncmos        
  V  : Canonical Fock matrix             fckc1           nmodes, ncmos , ncmos
     : - first derivatives                                                    
     : AO->LMO matrix                    vecl            nmos  , nbasis       
     : AO->LMO matrix                    vecl1           nmodes, nmos  , nbasis
     : - first derivatives                                                    
     : AO->CMO matrix                    vecc            ncmos , nbasis       
     : AO->CMO matrix                    vecc1           nmodes, ncmos , nbasis
     : - first derivatives                                                    
 =============================================================================
  [3] TENSOR DIMENSIONS AND STORAGE CONVENTIONS                               
 =============================================================================
   All of the nonsymmetric tensors exept for DMA and DPOL analyses are in full 
   format. DMA and DPOL are stored in reduced formats. Other symmetric tensor 
   elements are stored in triangular matrix format (Fock matrices + derivs, cu
   bic anharmonic constants).

   Reduced format for electrostatic descriptors:
                                                                              
             CHG(i) .
                    1
             DIP(i) X   Y   Z
                    1   2   3
             QAD(i) XX  YY  ZZ  XY  XZ  YZ
                    1   2   3   4   5   6
             OCT(i) XXX YYY ZZZ XXY XXZ XYY YYZ XZZ YZZ XYZ
                    1   2   3   4   5   6   7   8   9   10
             POL(i) XX  XY  XZ  YX  YY  YZ  ZX  ZY  ZZ 
                    1   2   3   4   5   6   7   8   9
                                                                              
   * Turn of indices in the case of DMA distributions is identical as GAMESS  
   * Turn of indices in distributed polarizabilities is different as in GAMESS
 =============================================================================
 Notes:                                                                       
   1) -ALL- the data in the instance of Frag are stored in ATOMIC UNITS!!!    
                                                                              
       [ Angstrom, kcal/mole, cm-1, C, amu, ... ] ----> [a.u.]                
                                                                              
 References:                                                                  
  * add Gordon's papers, M.Cho's papers, R.W.Góra's papers and your papers    
  * add MDAnalysis reference                                                  
  * add LibBBG, Coulomb.py, PyQuante, Gaussian09 and GAMESS references        
 =============================================================================
  Globulion@                                                                  
"""
    ### built-in parameters
    params = ['water'  ,'water2' ,'nma'     ,'nma-d7'  ,'mescn'   ,'meoac'   ,
              'meoh'   ,'dmso'   ,'dcm'     ,'chcl3'   ,'na+'     ,'so3--'   ,
              'me-so3-',]
              
    def __init__(self,file=None):
        self.__file = file
        self._create()
        self._nlines = lambda n: n/5+bool(n%5)
        if file is not None: self.__call__(file)
        
    def __call__(self,file):
        """parse the parameters from a file"""
        if file in self.params:
           self.__file = os.environ['SLV_PATH'] + '/frg/%s/%s.frg' % (file,file)
           if not os.path.isfile(self.__file):
              message = " No such file: <%s>. There is no *.frg file for that molecule,\n" % self.__file
              message+= " the name is misspelled or check if the $SLV_PATH variable is set properly\n"
              message+= " (echo $SLV_PATH)" 
              raise IOError, message
        else: self.__file = file
        filef = open(self.__file,'r')
        text = filef.read()
        filef.close()
        templ = re.compile(r'\[',re.DOTALL)
        sections = re.split(templ,text)[1:]
        del text
        for section in sections:
            self._read(section)
        return
    
    def __repr__(self):
        """print the information about the status"""
        log = ' \n'
        log+= ' ================================= \n'
        log+= ' SLV SOLVATOCHROMIC EFP PARAMETERS \n'
        log+= ' ================================= \n\n'
        log+= ' MOLECULE SPECIFICATION \n'
        log+= ' name     %s    \n' %  self.__name
        log+= ' method   %s/%s \n' % (self.__method, self.__basis)
        log+= ' natoms   %s    \n' %  self.__natoms
        log+= ' nbasis   %s    \n' %  self.__nbasis
        log+= ' nmos     %s    \n' %  self.__nmos
        log+= ' nmodes   %s    \n' %  self.__nmodes
        log+= ' \n'
        return str(log)
    
    # public
    
    def set(self,mol=None,anh=None,frag=None,
            dma=None,**kwargs):
        """set the parameters providing appropriate objects"""
        # molecular structure
        if mol is not None:
           self.__name = mol.get_name()
           self.__pos  = mol.get_pos()
           #self.__atoms = mol.get_atoms() # dopisz do Molecule class!
           self.__natoms = len(self.__pos)
           self.__nmodes = 3 * self.__natoms - 6
           self.__nmos   = int(sum(mol.get_atno())/2)
           self.__method = mol.get_method()
           self.__basis  = mol.get_basis()
           self.__nbasis = len( mol.get_bfs() )
           self.__ncmos  = self.__nbasis
           self.__atno   = mol.get_atno()
           self.__atms   = mol.get_atms() * self.AmuToElectronMass
           self.__npol   = self.__nmos
        # electrostatic data
        if dma is not None:
           assert not dma.if_traceless(), "DMA object is ALREADY traceless!!!"
           self.__pos    = dma.get_pos()
           self.__rdma   = dma.get_origin()
           self.__ndma   = len(self.__pos)
           self.__dmac   = dma[0]
           self.__dmad   = dma[1]
           self.__dmaq   = dma[2]
           self.__dmao   = dma[3]
        # anharmonic file object (FREQ)
        if anh is not None:
           assert anh.if_w(), 'Anharmonic object is in wrong units! Supply anh.w() object'
           self.__redmass = anh.redmass
           self.__freq = anh.freq
           self.__lvec = self._tr_lvec(anh.L,self.__nmodes,self.__natoms)
           self.__gijk = anh.K3
        # EFP fragment parameters
        if frag is not None:
           self._parse_dict( frag.get() )
        # Other special data
        if kwargs: self.set_par(**kwargs)
        return
    
    def set_par(self,**kwargs):
        """set the particular parameters"""
        for key, arg in kwargs.items():
            if key == 'rdma' : self.__rdma  = arg
            if key == 'dmac' : self.__dmac  = arg
            if key == 'dmad' : self.__dmad  = arg
            if key == 'dmaq' : self.__dmaq  = arg
            if key == 'dmao' : self.__dmao  = arg
            if key == 'fock' : self.__fock  = arg
            if key == 'fock1': self.__fock1 = arg
            if key == 'vecl' : self.__vecl  = arg
            if key == 'vecl1': self.__vecl1 = arg
            if key == 'lmoc' : self.__lmoc  = arg
            if key == 'lmoc1': self.__lmoc1 = arg
            if key == 'chlpg': self.__chlpg = arg
            if key == 'esp'  : self.__esp   = arg
            if key == 'dmac1': self.__dmac1 = arg
            if key == 'dmad1': self.__dmad1 = arg
            if key == 'dmaq1': self.__dmaq1 = arg
            if key == 'dmao1': self.__dmao1 = arg
            if key == 'dmac2': self.__dmac2 = arg
            if key == 'dmad2': self.__dmad2 = arg
            if key == 'dmaq2': self.__dmaq2 = arg
            if key == 'dmao2': self.__dmao2 = arg
        return
            
    def get(self):
        """returns dictionary with parameters"""
        return self._make_dict()
    
    def get_name(self):
        """return the name of the fragment"""
        return self.__name
    
    def get_natoms(self):
        """return the number of atoms"""
        return self.__natoms
    
    def get_bfs(self):
        """return basis set object"""
        mol = MakeMol(self.__atno,self.__pos,name=self.__name)
        bfs = getbasis(mol,self.__basis)
        return bfs
    
    def get_pos(self):
        """return position of a fragment"""
        return self.__pos
    
    def get_traceless(self):
        """return traceless quadrupoles and octupoles"""
        dmaq, dmao = tracls(self.__dmaq, self.__dmao)
        return dmaq, dmao
    
    def get_traceless_1(self,ravel=False):
        """return traceless 1st derivatives wrt modes of quadrupoles and octupoles"""
        self.__dmaq1 = self.__dmaq1.ravel()
        self.__dmao1 = self.__dmao1.ravel()
        dmaq1, dmao1 = tracl1(self.__dmaq1, self.__dmao1, self.__nmodes, self.__ndma)
        if not ravel:
           self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes,self.__ndma,6)
           self.__dmao1 = self.__dmao1.reshape(self.__nmodes,self.__ndma,10)
        return dmaq1, dmao1
    
    def get_traceless_2(self):
        """return traceless  2nd derivatives wrt modes of quadrupoles and octupoles"""
        dmaq2, dmao2 = tracls(self.__dmaq2, self.__dmao2)
        return dmaq2, dmao2
    
    def write(self,file='slv.frg',par=None):
        """writes the parameters into a file"""
        f = open(file,'w')
        if par is not None:
           pass
        # basic molecular data
        if self.__name   is not None: self._write_preambule(f)
        if self.__pos    is not None: self._write_pos(f)
        if self.__origin is not None: self._write_origin(f)
        if self.__atno   is not None: self._write_atno(f)
        if self.__atms   is not None: self._write_atms(f)
        # population analysis
        if self.__esp   is not None: self._write_esp(f)
        if self.__chlpg is not None: self._write_chlpg(f)
        if self.__rdma  is not None: self._write_rdma(f)
        if self.__dmac  is not None: self._write_dmac(f)
        if self.__dmad  is not None: self._write_dmad(f)
        if self.__dmaq  is not None: self._write_dmaq(f)
        if self.__dmao  is not None: self._write_dmao(f)
        if self.__dmac1 is not None: self._write_dmac1(f)
        if self.__dmad1 is not None: self._write_dmad1(f)
        if self.__dmaq1 is not None: self._write_dmaq1(f)
        if self.__dmao1 is not None: self._write_dmao1(f)
        if self.__dmac2 is not None: self._write_dmac2(f)
        if self.__dmad2 is not None: self._write_dmad2(f)
        if self.__dmaq2 is not None: self._write_dmaq2(f)
        if self.__dmao2 is not None: self._write_dmao2(f)
        if self.__rpol  is not None: self._write_rpol(f)
        if self.__dpol  is not None: self._write_dpol(f)
        if self.__dpol1 is not None: self._write_dpol1(f)
        # frequency analysis
        if self.__redmass  is not None: self._write_redmass(f)
        if self.__freq     is not None: self._write_freq(f)
        if self.__lvec     is not None: self._write_lvec(f)
        if self.__gijk     is not None: self._write_gijk(f)
        # EFP parameters
        if self.__lmoc  is not None: self._write_lmoc(f)
        if self.__lmoc1 is not None: self._write_lmoc1(f)
        if self.__fock  is not None: self._write_fock(f)
        if self.__fock1 is not None: self._write_fock1(f)
        if self.__vecl  is not None: self._write_vecl(f)
        if self.__vecl1 is not None: self._write_vecl1(f)
        if self.__fckc  is not None: self._write_fckc(f)
        if self.__fckc1 is not None: self._write_fckc1(f)
        if self.__vecc  is not None: self._write_vecc(f)
        if self.__vecc1 is not None: self._write_vecc1(f)
        f.close()
        return
    
    def translate(self,transl):
        """translate tensors by <transl> cartesian displacement vector"""
        if self.__pos   is not None: self.__pos   += transl
        if self.__lmoc  is not None: self.__lmoc  += transl
        if self.__rpol  is not None: self.__rpol  += transl
        if self.__rdma  is not None: self.__rdma  += transl
        return
    
    def rotate(self,rot):
        """rotate the tensors by <rot> unitary matrix"""
        # transform the atomic position static and dynamic information
        if self.__pos   is not None:
           self.__pos    = dot(self.__pos , rot)
        if self.__lmoc  is not None:
           self.__lmoc   = dot(self.__lmoc, rot)
        if self.__lmoc1 is not None:
           self.__lmoc1  = dot(self.__lmoc1,rot)
        if self.__lvec  is not None:
           self.__lvec   = dot(self.__lvec, rot)
        if self.__rpol  is not None:
           self.__rpol   = dot(self.__rpol, rot)
        if self.__rdma  is not None:
           self.__rdma   = dot(self.__rdma, rot)
        # transform dipoles, quadrupoles and octupoles!
        if self.__dmac  is not None:
           self.__dmad, self.__dmaq, self.__dmao = \
           rotdma(self.__dmad,self.__dmaq,self.__dmao,rot)
        if self.__dmac1 is not None:
           self.__dmad1 = self.__dmad1.ravel()
           self.__dmaq1 = self.__dmaq1.ravel()
           self.__dmao1 = self.__dmao1.ravel()
           self.__dmad1, self.__dmaq1, self.__dmao1 = \
           rotdm1(self.__dmad1,self.__dmaq1,self.__dmao1,rot,
                  self.__nmodes,self.__ndma)
           self.__dmad1 = self.__dmad1.reshape(self.__nmodes,self.__ndma,3)
           self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes,self.__ndma,6)
           self.__dmao1 = self.__dmao1.reshape(self.__nmodes,self.__ndma,10)
        if self.__dmac2 is not None:
           self.__dmad2, self.__dmaq2, self.__dmao2 = \
           rotdma(self.__dmad2,self.__dmaq2,self.__dmao2,rot)
        # transform distributed polarizabilities!
        if self.__dpol   is not None:
           for i in xrange(self.__npol):
               self.__dpol[i] = dot(transpose(rot),dot(self.__dpol[i],rot))
        if self.__dpol1  is not None:
           for i in xrange(self.__nmodes):
               for j in xrange(self.__npol):
                   self.__dpol1[i,j]  = dot(transpose(rot),dot(self.__dpol1[i,j],rot))
        # transform the wave function!
        if self.__vecl  is not None:
           bfs = self.get_bfs()
           typs= bfs.get_bfst().sum(axis=1)
           self.__vecl = vecrot(self.__vecl,rot,typs)
           if self.__vecc  is not None:
              self.__vecc  = vecrot(self.__vecc,rot,typs)
           if self.__vecl1 is not None:
              self.__vecl1 = vc1rot(self.__vecl1,rot,typs)
           if self.__vecc1 is not None:
              self.__vecc1 = vc1rot(self.__vecc1,rot,typs)
        return
    
    def sup(self,xyz,suplist=None):
        """superimpose structures <str> and <self.__pos>. Return rms in A.U."""
        s = svd_sup()
        if suplist is None: s.set(xyz,self.__pos)
        else:               s.set(xyz[suplist],self.__pos[suplist])
        s.run()
        rms         = s.get_rms()
        rot, transl = s.get_rotran()
        # perform transformations
        #self.__pos       = s.get_transformed()
        self.__pos = dot(self.__pos , rot) + transl
        if self.__lmoc  is not None: self.__lmoc   = dot(self.__lmoc , rot) + transl
        if self.__rpol  is not None: self.__rpol   = dot(self.__rpol , rot) + transl
        if self.__rdma  is not None: self.__rdma   = dot(self.__rdma , rot) + transl
        if self.__lmoc1 is not None: self.__lmoc1  = dot(self.__lmoc1, rot)
        if self.__lvec  is not None: self.__lvec   = dot(self.__lvec , rot)
        # transform dipoles, quadrupoles and octupoles!
        if self.__dmac  is not None:
           self.__dmad, self.__dmaq, self.__dmao = \
           rotdma(self.__dmad,self.__dmaq,self.__dmao,rot)
        if self.__dmac1 is not None:
           self.__dmad1 = self.__dmad1.ravel()
           self.__dmaq1 = self.__dmaq1.ravel()
           self.__dmao1 = self.__dmao1.ravel()
           self.__dmad1, self.__dmaq1, self.__dmao1 = \
           rotdm1(self.__dmad1,self.__dmaq1,self.__dmao1,rot,
                  self.__nmodes,self.__ndma)
           self.__dmad1 = self.__dmad1.reshape(self.__nmodes,self.__ndma,3)
           self.__dmaq1 = self.__dmaq1.reshape(self.__nmodes,self.__ndma,6)
           self.__dmao1 = self.__dmao1.reshape(self.__nmodes,self.__ndma,10)
        if self.__dmac2 is not None:
           self.__dmad2, self.__dmaq2, self.__dmao2 = \
           rotdma(self.__dmad2,self.__dmaq2,self.__dmao2,rot)
        # transform distributed polarizabilities!
        if self.__dpol   is not None:
           for i in xrange(self.__npol):
               self.__dpol[i] = dot(transpose(rot),dot(self.__dpol[i],rot))
        if self.__dpol1  is not None:
           for i in xrange(self.__nmodes):
               for j in xrange(self.__npol):
                   self.__dpol1[i,j]  = dot(transpose(rot),dot(self.__dpol1[i,j],rot))
        # - wave function
        if self.__vecl  is not None:
           bfs = self.get_bfs()
           typs= bfs.get_bfst().sum(axis=1)
           self.__vecl = vecrot(self.__vecl,rot,typs)
           if self.__vecc  is not None:
              self.__vecc  = vecrot(self.__vecc,rot,typs)
           if self.__vecl1 is not None:
              self.__vecl1 = vc1rot(self.__vecl1,rot,typs)
           if self.__vecc1 is not None:
              self.__vecc1 = vc1rot(self.__vecc1,rot,typs)
        return rms
    
    def copy(self):
        """return a deepcopy of me!"""
        return copy.deepcopy(self)
    
    def reorder(self,ord,dma=True):
        """Reorder the turn of atoms in the parameter objects according to <ord>.
<dma> option is designed for the case when dma is not atomic and mid-bond-based (then False).
In other words - dma=True if DMA, CAMM and CBAMM are used. False if LMTP and others ..."""
        LIST1 = self.get_bfs().get_bfsl() + 1
        sim = zeros((self.__natoms,2),int)
        sim[:,0] = arange(1,self.__natoms+1,1)
        sim[:,1] = ord
        # reorder atomic positions
        self.__pos = reorder(self.__pos,sim)
        self.__atno= reorder(self.__atno,sim)
        self.__atms= reorder(self.__atms,sim)
        # reorder multipole moments
        if dma:
           if self.__rdma is not None:
              self.__rdma[:self.__natoms] = reorder(self.__rdma[:self.__natoms],sim)
           if self.__dmac is not None:
              self.__dmac[:self.__natoms] = reorder(self.__dmac[:self.__natoms],sim)
              self.__dmad[:self.__natoms] = reorder(self.__dmad[:self.__natoms],sim)
              self.__dmaq[:self.__natoms] = reorder(self.__dmaq[:self.__natoms],sim)
              self.__dmao[:self.__natoms] = reorder(self.__dmao[:self.__natoms],sim)
           #
           if self.__dmac1 is not None:
              self.__dmac1[:self.__natoms] = reorder(self.__dmac1[:self.__natoms],sim,axis=1)
              self.__dmad1[:self.__natoms] = reorder(self.__dmad1[:self.__natoms],sim,axis=1)
              self.__dmaq1[:self.__natoms] = reorder(self.__dmaq1[:self.__natoms],sim,axis=1)
              self.__dmao1[:self.__natoms] = reorder(self.__dmao1[:self.__natoms],sim,axis=1)
              #
              self.__dmac2[:self.__natoms] = reorder(self.__dmac2[:self.__natoms],sim,axis=0)
              self.__dmad2[:self.__natoms] = reorder(self.__dmad2[:self.__natoms],sim,axis=0)
              self.__dmaq2[:self.__natoms] = reorder(self.__dmaq2[:self.__natoms],sim,axis=0)
              self.__dmao2[:self.__natoms] = reorder(self.__dmao2[:self.__natoms],sim,axis=0)
        # reorder eigenvectors
        if self.__lvec is not None: 
           self.__lvec = reorder(self.__lvec,sim,axis=1)
        # reorder the wavefunction
        if self.__vecl is not None:
           self._reorder_wfn(sim,LIST1)
        return
    
    
    # protected
    
    def _create(self):
        """creates the namespace of parameter variables"""
        self.__name ,self.__natoms, self.__nbasis = None, None, None
        self.__nmos ,self.__nmodes, self.__basis  = None, None, None
        self.__atoms,self.__shortname             = None, None
        self.__pos  ,self.__origin, self.__nsites = None, None, None
        self.__atno ,self.__atms  , self.__mode   = None, None, None
        #
        self.__redmass, self.__freq, self.__lvec = None, None, None
        self.__gijk                              = None
        #
        self.__fock ,self.__lmoc ,self.__vecl = None, None, None
        self.__fock1,self.__lmoc1,self.__vecl1= None, None, None
        self.__ncmos,self.__vecc ,self.__vecc1= None, None, None
        self.__fckc ,self.__fckc1             = None, None
        #
        self.__esp  ,self.__chlpg,self.__npol = None, None, None
        self.__dpol ,self.__dpol1,self.__rpol = None, None, None
        #
        self.__rdma, self.__dmac, self.__dmad = None, None, None
        self.__dmaq, self.__dmao, self.__ndma = None, None, None
        self.__dmac1,self.__dmad1,self.__dmaq1= None, None, None
        self.__dmao1                          = None
        self.__dmac2,self.__dmad2,self.__dmaq2= None, None, None
        self.__dmao2                          = None
        #
        mol_names = ('name' ,'basis' ,'method','natoms'   ,'nbasis',
                     'nmos' ,'nmodes','atoms' ,'shortname','nsites',
                     'ncmos','npol'  ,'ndma'  ,'mode'     ,)
        sec_names = {  'mol': '[ molecule ]'                                        ,
                      'atno': '[ Atomic numbers ]'                                  ,
                      'atms': '[ Atomic masses ]'                                   ,
                       'pos': '[ Atomic coordinates ]'                              ,
                      'freq': '[ Harmonic frequencies ]'                            ,
                   'redmass': '[ Reduced masses ]'                                  ,
                      'lvec': '[ Mass-weighted eigenvectors ]'                      ,
                      'gijk': '[ Cubic anharmonic constants ]'                      ,
                    'origin': '[ DMTP origins ]'                                    ,
                       'esp': '[ ESP charges ]'                                     ,
                     'chlpg': '[ ChelpG charges ]'                                  ,
                      'rdma': '[ DMTP centers ]'                                    ,
                      'dmac': '[ DMTP charges ]'                                    ,
                      'dmad': '[ DMTP dipoles ]'                                    ,
                      'dmaq': '[ DMTP quadrupoles ]'                                ,
                      'dmao': '[ DMTP octupoles ]'                                  ,
                     'dmac1': '[ DMTP charges - first derivatives ]'                ,
                     'dmad1': '[ DMTP dipoles - first derivatives ]'                ,
                     'dmaq1': '[ DMTP quadrupoles - first derivatives ]'            ,
                     'dmao1': '[ DMTP octupoles - first derivatives ]'              ,
                     'dmac2': '[ DMTP charges - second derivatives wrt mode]'       ,
                     'dmad2': '[ DMTP dipoles - second derivatives wrt mode]'       ,
                     'dmaq2': '[ DMTP quadrupoles - second derivatives wrt mode]'   ,
                     'dmao2': '[ DMTP octupoles - second derivatives wrt mode]'     ,
                      'rpol': '[ Polarizable centers ]'                             ,
                      'dpol': '[ Distributed polarizabilities ]'                    ,
                     'dpol1': '[ Distributed polarizabilities - first derivatives ]',
                     'lmoc' : '[ LMO centroids ]'                                   ,
                     'lmoc1': '[ LMO centroids - first derivatives ]'               ,
                     'fock' : '[ Fock matrix ]'                                     ,
                     'fock1': '[ Fock matrix - first derivatives ]'                 ,
                     'fckc' : '[ Canonical Fock matrix ]'                           ,
                     'fckc1': '[ Canonical Fock matrix - first derivatives ]'       ,
                     'vecl' : '[ AO->LMO matrix ]'                                  ,
                     'vecl1': '[ AO->LMO matrix - first derivatives ]'              ,
                     'vecc' : '[ AO->CMO matrix ]'                                  ,
                     'vecc1': '[ AO->CMO matrix - first derivatives ]'              ,}
        self.__mol_names = mol_names
        self.__sec_names = sec_names
        return

    def _parse_dict(self,par):
        """save the memorials from par dictionary"""
        for key, val in par.items():
            if   key == 'fock' : self.__fock  = val
            elif key == 'fock1': self.__fock1 = val
            elif key == 'lmoc' : self.__lmoc  = val
            elif key == 'lmoc1': self.__lmoc1 = val
            elif key == 'vecl' : self.__vecl  = val
            elif key == 'vecl1': self.__vecl1 = val
            elif key == 'vecc' : self.__vecc  = val
            elif key == 'vecc1': self.__vecc1 = val
            elif key == 'fckc' : self.__fckc  = val
            elif key == 'fckc1': self.__fckc1 = val
        return

    def _make_dict(self):
        """returns dictionary with parameters"""
        par = {}
        # basic molecular data
        if self.__pos    is not None: par['pos'   ] = self.__pos
        if self.__origin is not None: par['origin'] = self.__origin
        if self.__atno   is not None: par['atno'  ] = self.__atno
        if self.__atms   is not None: par['atms'  ] = self.__atms
        if self.__ndma   is not None: par['ndma'  ] = self.__ndma
        if self.__npol   is not None: par['npol'  ] = self.__npol
        if self.__nmodes is not None: par['nmodes'] = self.__nmodes
        if self.__natoms is not None: par['natoms'] = self.__natoms
        if self.__name   is not None: par['name'  ] = self.__name
        if self.__basis  is not None: par['basis' ] = self.__basis
        if self.__nbasis is not None: par['nbasis'] = self.__nbasis
        if self.__nmos   is not None: par['nmos'  ] = self.__nmos
        if self.__mode   is not None: par['mode'  ] = self.__mode
        # frequency analysis
        if self.__redmass  is not None: par['redmass'] = self.__redmass
        if self.__freq     is not None: par['freq'   ] = self.__freq
        if self.__lvec     is not None: par['lvec'   ] = self.__lvec
        if self.__gijk     is not None: par['gijk'   ] = self.__gijk
        # EFP parameters
        if self.__lmoc  is not None: par['lmoc' ] = self.__lmoc
        if self.__lmoc1 is not None: par['lmoc1'] = self.__lmoc1
        if self.__fock  is not None: par['fock' ] = self.__fock
        if self.__fock1 is not None: par['fock1'] = self.__fock1
        if self.__vecl  is not None: par['vecl' ] = self.__vecl
        if self.__vecl1 is not None: par['vecl1'] = self.__vecl1
        if self.__vecc  is not None: par['vecc' ] = self.__vecc
        if self.__vecc1 is not None: par['vecc1'] = self.__vecc1
        if self.__fckc  is not None: par['fckc' ] = self.__fckc
        if self.__fckc1 is not None: par['fckc1'] = self.__fckc1
        # Population analysis
        if self.__esp   is not None: par['esp'  ] = self.__esp
        if self.__chlpg is not None: par['chlpg'] = self.__chlpg
        if self.__rdma  is not None: par['rdma' ] = self.__rdma
        if self.__dmac  is not None: par['dmac' ] = self.__dmac
        if self.__dmad  is not None: par['dmad' ] = self.__dmad
        if self.__dmaq  is not None: par['dmaq' ] = self.__dmaq
        if self.__dmao  is not None: par['dmao' ] = self.__dmao
        if self.__dmac1 is not None: par['dmac1'] = self.__dmac1
        if self.__dmad1 is not None: par['dmad1'] = self.__dmad1
        if self.__dmaq1 is not None: par['dmaq1'] = self.__dmaq1
        if self.__dmao1 is not None: par['dmao1'] = self.__dmao1
        if self.__dmac2 is not None: par['dmac2'] = self.__dmac2
        if self.__dmad2 is not None: par['dmad2'] = self.__dmad2
        if self.__dmaq2 is not None: par['dmaq2'] = self.__dmaq2
        if self.__dmao2 is not None: par['dmao2'] = self.__dmao2
        if self.__rpol  is not None: par['rpol' ] = self.__rpol
        if self.__dpol  is not None: par['dpol' ] = self.__dpol
        if self.__dpol1 is not None: par['dpol1'] = self.__dpol1
        #
        return par
        
    def _reorder_wfn(self,sim,LIST1):
        """reorders the turn of basis functions according to <ord> list of atomic interchange"""
        P1, P2 = [], []
        #
        T1 = transpose(self.__vecl,axes=(1,0))
        #
        b=[]; n=1; curr= LIST1[0]
        for i in range(len(LIST1)-1):
            if curr==LIST1[i+1]: n+=1
            else:
                b.append(n)
                n=1
            curr=LIST1[i+1]
        b.append(n)
        #
        c = [ sum(b[:x]) for x in range(len(b)) ]
        c.append(len(LIST1))
        #
        for i in range(self.__natoms):
            P1.append( T1[c[i]:c[i+1]] )
        P1 = array(P1)
        #
        N1 = zeros(self.__natoms,object)
        #
        for i,j in sim:
            N1[i-1] = P1[j-1]
        #
        self.__vecl = con(tuple(N1))
        self.__vecl = transpose(self.__vecl,axes=(1,0))
        ### DERIVATIVES
        if self.__vecl1 is not None:
           T2 = transpose(self.__vecl1,axes=(2,1,0))
           #
           for i in range(self.__natoms):
               P2.append( T2[c[i]:c[i+1]] )
           P2 = array(P2)
           #
           N2 = zeros(self.__natoms,object)
           for i,j in sim:
               N2[i-1] = P2[j-1]
           #
           self.__vecl1 = con(tuple(N2))
           self.__vecl1 = transpose(self.__vecl1,axes=(2,1,0))
        return
    
    def _tr_lvec(self,lvec,nmodes,natoms):
        """Transpose axes and reshape the L-matrix from FREQ instance"""
        return lvec.transpose().reshape(nmodes,natoms,3)
    
    # --------------------------------------------------------- #
    #            R E A D I N G    P R O C E D U R E S           #
    # --------------------------------------------------------- #
    
    def _read(self,section):
        """read the appropriate field from section and save it in SLVPAR instance"""
        section = section.split('\n')
        querry = section[0]
        for key, val in self.__sec_names.items():
            if val[2:] in querry: break
        ### basic molecular specification
        if key == 'mol':
            data = section[1:]
            for line in data:
                if '=' in line:
                    name, arg = line.split('=')
                    name = name.strip(); arg = arg.strip()
                    if name == 'name':
                        self.__name = arg
                    if name == 'shortname':
                        self.__shortname = arg
                    if name == 'atoms':
                        self.__atoms = arg.split(',')
                    if name == 'natoms':
                        self.__natoms = int(arg)
                    if name == 'nsites':
                        self.__nsites = int(arg)
                    if name == 'basis':
                        self.__method, self.__basis = arg.split('/')
                    if name == 'nbasis':
                        self.__nbasis = int(arg)
                    if name == 'nmos':
                        self.__nmos = int(arg)
                    if name == 'ncmos':
                        self.__ncmos = int(arg)
                    if name == 'nmodes':
                        self.__nmodes = int(arg)
                    if name == 'npol':
                        self.__npol = int(arg)
                    if name == 'ndma':
                        self.__ndma = int(arg)
                    if name == 'mode':
                        self.__mode = int(arg)
                
        ### more advanced information
        else:
            data = [] ; N = int(querry.split('N=')[1])
            n = self._nlines(N)
            for i in xrange(n): data += section[i+1].split()
            data = array(data,dtype=float64)
            # ------------------------------------ MOL --------------------------------------------
            # Atomic coordinates
            if   key == 'pos':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic coordinates ]!'
                 assert self.__natoms == N/3, merror
                 data = data.reshape(self.__natoms,3)
                 self.__pos = data
            # DMTP origins
            elif key == 'origin':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP origins ]!'
                 assert 1==1, merror
                 data = data.reshape(N/3,3)
                 self.__origin = data
            # Atomic numbers
            elif key == 'atno':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic numbers ]!'
                 assert self.__natoms == N, merror
                 self.__atno = array(data,int)
            # Atomic masses
            elif key == 'atms':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic masses ]!'
                 assert self.__natoms == N, merror
                 self.__atms = data
            # ------------------------------------ DMTP -------------------------------------------
            # ESP charges
            elif key == 'esp':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ ESP charges ]!'
                 assert self.__natoms == N, merror
                 self.__esp = data
            # ChelpG charges
            elif key == 'chlpg':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ ChelpG charges ]!'
                 assert self.__natoms == N, merror
                 self.__chlpg = data
            # DMA/CAMM
            # distributed center coordinates
            elif key == 'rdma':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP centers ]!'
                 assert self.__ndma == N/3, merror
                 data = data.reshape(self.__ndma,3)
                 self.__rdma = data
            # distributed charges
            elif key == 'dmac':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma == N, merror
                 self.__dmac = data
            # distributed dipoles
            elif key == 'dmad':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma == N/3, merror
                 data = data.reshape(self.__ndma,3)
                 self.__dmad = data
            # distributed quadrupoles
            elif key == 'dmaq':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma == N/6, merror
                 data = data.reshape(self.__ndma,6)
                 self.__dmaq = data
            # distributed octupoles
            elif key == 'dmao':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma == N/10, merror
                 data = data.reshape(self.__ndma,10)
                 self.__dmao = data
            # - DMTP derivatives -
            # first derivatives
            # distributed charge 1st derivatives
            elif key == 'dmac1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma*self.__nmodes == N, merror
                 self.__dmac1 = data.reshape(self.__nmodes,self.__ndma)
                 temp = sqrt(self.__redmass)[:,newaxis]
                 self.__dmac1 = temp * self.__dmac1
            # distributed dipole 1st derivatives
            elif key == 'dmad1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma*self.__nmodes == N/3, merror
                 data = data.reshape(self.__nmodes,self.__ndma,3)
                 temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                 self.__dmad1 = temp * data
            # distributed quadrupoles 1st derivatives
            elif key == 'dmaq1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma*self.__nmodes == N/6, merror
                 data = data.reshape(self.__nmodes,self.__ndma,6)
                 temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                 self.__dmaq1 = temp * data
            # distributed octupoles 1st derivatives
            elif key == 'dmao1':
                 merror = 'nmodes and ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma*self.__nmodes == N/10, merror
                 data = data.reshape(self.__nmodes,self.__ndma,10)
                 temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                 self.__dmao1 = temp * data
            # second derivatives
            # distributed charges 2nd derivatives
            elif key == 'dmac2':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP charges ]!'
                 assert self.__ndma == N, merror
                 temp = self.__redmass[self.__mode-1]
                 self.__dmac2 = temp * data
            # distributed dipoles 2nd derivatives
            elif key == 'dmad2':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP dipoles ]!'
                 assert self.__ndma == N/3, merror
                 data = data.reshape(self.__ndma,3)
                 temp = self.__redmass[self.__mode-1]
                 self.__dmad2 = temp * data
            # distributed quadrupoles 2nd derivatives
            elif key == 'dmaq2':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP quadrupoles ]!'
                 assert self.__ndma == N/6, merror
                 data = data.reshape(self.__ndma,6)
                 temp = self.__redmass[self.__mode-1]
                 self.__dmaq2 = temp * data
            # distributed octupoles 2nd derivatives
            elif key == 'dmao2':
                 merror = 'ndma in section [ molecule ] '
                 merror+= 'is not consistent with section [ DMTP octupoles ]!'
                 assert self.__ndma == N/10, merror
                 data = data.reshape(self.__ndma,10)
                 temp = self.__redmass[self.__mode-1]
                 self.__dmao2 = temp * data
            # CABMM
            elif key == 'cabmm':
                 pass
            # LMTP
            elif key == 'lmtp':
                 pass
            # ------------------------------------ DPOL -------------------------------------------
            elif key == 'rpol':
                 merror = None
                 self.__rpol = data.reshape(self.__npol,3)
            elif key == 'dpol':
                 merror = 'npol in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distributed polarizabilities ]!'
                 assert self.__npol == (N/9), merror
                 self.__dpol = data.reshape(self.__npol,3,3)
            elif key == 'dpol1':
                 merror = 'npol and nmodes in section [ molecule ] '
                 merror+= 'is not consistent with section [ Distributed polarizabilities ]!'
                 assert self.__npol == (N/(9*self.__nmodes)), merror
                 self.__dpol1 = data.reshape(self.__nmodes,self.__npol,3,3)
                 temp = sqrt(self.__redmass)[:,newaxis,newaxis,newaxis]
                 self.__dpol1 = temp * self.__dpol1
            # ------------------------------------ FREQ -------------------------------------------
            # Harmonic frequencies
            elif key == 'freq':
                 merror = 'nmodes in section [ molecule ] '
                 merror+= 'is not consistent with section [ Harmonic frequencies ]!'
                 assert self.__nmodes == N, merror
                 self.__freq = data
            # Reduced masses
            elif  key == 'redmass':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Reduced masses ]!'
                  assert self.__nmodes == N, merror
                  self.__redmass = data
            # Mass-weighted eigenvectors
            elif  key == 'lvec':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Mass-weighted eigenvectors ]!'
                  assert self.__nmodes == N/(self.__natoms*3), merror
                  self.__lvec = data.reshape(self.__nmodes,self.__natoms,3)
            # Cubic anharmonic constants
            elif  key == 'gijk':
                  merror = 'nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ Cubic anharmonic constants ]!'
                  assert 1==1, merror
                  gijk = zeros((self.__nmodes,self.__nmodes,self.__nmodes),dtype=float64)
                  K = 0
                  for i in xrange(self.__nmodes):
                      for j in xrange(i+1):
                          for k in xrange(j+1):
                              d = data[K]
                              gijk[i,j,k] = d
                              gijk[i,k,j] = d
                              gijk[j,i,k] = d
                              gijk[j,k,i] = d
                              gijk[k,i,j] = d
                              gijk[k,j,i] = d
                              K += 1                
                  self.__gijk = gijk
            # ------------------------------------ EFP --------------------------------------------
            # LMO centroids
            elif  key == 'lmoc':
                  merror = 'nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ LMO centroids ]!'
                  assert self.__nmos == N/3, merror
                  data = data.reshape(self.__nmos,3)
                  self.__lmoc = data
            # LMO centroids first derivatives
            elif  key == 'lmoc1':
                  merror = 'nmos and nmodes in section [ molecule ] '
                  merror+= 'is not consistent with section [ LMO centroids - first derivatives ]!'
                  assert self.__nmos == N/(3*self.__nmodes), merror
                  data = data.reshape(self.__nmodes,self.__nmos,3)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                  data = temp * data
                  self.__lmoc1 = data
            # Fock matrix
            elif  key == 'fock':
                  merror = 'nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Fock matrix ]!'
                  assert self.__nmos == (math.sqrt(1+8*N)-1)/2, merror
                  fock = zeros((self.__nmos,self.__nmos),dtype=float64)
                  K = 0
                  for i in xrange(self.__nmos):
                      for j in xrange(i+1):
                          d = data[K]
                          fock[i,j] = d
                          fock[j,i] = d
                          K += 1
                  self.__fock = fock
            # Canonical Fock matrix
            elif  key == 'fckc':
                  merror = 'ncmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Canonical Fock matrix ]!'
                  assert self.__ncmos == (math.sqrt(1+8*N)-1)/2, merror
                  fckc = zeros((self.__ncmos,self.__ncmos),dtype=float64)
                  K = 0
                  for i in xrange(self.__ncmos):
                      for j in xrange(i+1):
                          d = data[K]
                          fckc[i,j] = d
                          fckc[j,i] = d
                          K += 1
                  self.__fckc = fckc
            # Fock matrix first derivatives
            elif  key == 'fock1':
                  merror = 'nbasis and nmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Fock matrix - first derivatives ]!'
                  assert self.__nmos == (math.sqrt(1+8*(N/self.__nmodes))-1)/2, merror
                  fock1 = zeros((self.__nmodes,self.__nmos,self.__nmos),dtype=float64)
                  data = data.reshape(self.__nmodes,N/self.__nmodes)
                  for i in xrange(self.__nmodes):
                      K = 0
                      for j in xrange(self.__nmos):
                          for k in xrange(j+1):
                              d = data[i,K]
                              fock1[i,j,k] = d
                              fock1[i,k,j] = d
                              K += 1
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                  fock1= temp * fock1
                  self.__fock1 = fock1
            # Canonical Fock matrix first derivatives
            elif  key == 'fckc1':
                  merror = 'nbasis and ncmos in section [ molecule ] '
                  merror+= 'is not consistent with section [ Canonical Fock matrix - first derivatives ]!'
                  assert self.__ncmos == (math.sqrt(1+8*(N/self.__nmodes))-1)/2, merror
                  fckc1 = zeros((self.__nmodes,self.__ncmos,self.__ncmos),dtype=float64)
                  data = data.reshape(self.__nmodes,N/self.__nmodes)
                  for i in xrange(self.__nmodes):
                      K = 0
                      for j in xrange(self.__ncmos):
                          for k in xrange(j+1):
                              d = data[i,K]
                              fckc1[i,j,k] = d
                              fckc1[i,k,j] = d
                              K += 1
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                  fckc1= temp * fckc1
                  self.__fckc1 = fckc1
            # AO to LMO transformation matrix
            elif  key == 'vecl':
                  merror = 'nmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->LMO matrix ]!'
                  assert self.__nmos == N/self.__nbasis, merror
                  data = data.reshape(self.__nmos,self.__nbasis)
                  self.__vecl = data
            # AO to CMO transformation matrix
            elif  key == 'vecc':
                  merror = 'ncmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->CMO matrix ]!'
                  assert self.__ncmos == N/self.__nbasis, merror
                  data = data.reshape(self.__ncmos,self.__nbasis)
                  self.__vecc = data
            # AO to LMO transformation matrix first derivatives
            elif  key == 'vecl1':
                  merror = 'nmodes, nmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->LMO matrix - first derivatives ]!'
                  assert self.__nmos == N/(self.__nbasis*self.__nmodes)
                  data = data.reshape(self.__nmodes,self.__nmos,self.__nbasis)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                  data = temp * data
                  self.__vecl1 = data
            # AO to CMO transformation matrix first derivatives
            elif  key == 'vecc1':
                  merror = 'nmodes, ncmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->CMO matrix - first derivatives ]!'
                  assert self.__ncmos == N/(self.__nbasis*self.__nmodes)
                  data = data.reshape(self.__nmodes,self.__ncmos,self.__nbasis)
                  # multiply by sqrt(redmass)
                  assert self.__redmass is not None, 'No reduced masses supplied!'
                  temp = sqrt(self.__redmass)[:,newaxis,newaxis]
                  data = temp * data
                  self.__vecc1 = data
            # other not-programmed section
            else:
                raise Exception("Non-standard section name '%s' detected! Check the format of your file!"%key)
        return

    # --------------------------------------------------------- #
    #            W R I T I N G    P R O C E D U R E S           #
    # --------------------------------------------------------- #
    
    def _write_preambule(self,file):
        """write the preambule of the parameter file"""
        log    = ' %s\n' % self.__sec_names['mol'].ljust(40)
        log   += '   name       = %s\n'    % self.__name
        log   += '   basis      = %s/%s\n' % (self.__method, self.__basis)
        log   += '   natoms     = %s\n'    % self.__natoms
        log   += '   nbasis     = %s\n'    % self.__nbasis
        log   += '   nmodes     = %s\n'    % self.__nmodes
        if self.__mode is not None:
           log+= '   mode       = %s\n'    % self.__mode
        log   += '   nmos       = %s\n'    % self.__nmos
        log   += '   ncmos      = %s\n'    % self.__ncmos
        if self.__ndma is not None:
           log+= '   ndma       = %s\n'    % self.__ndma
        log   += ' \n'
        file.write(log)
        return
    
    def _write_pos(self,file):
        """write atomic coordinates"""
        natoms = self.__pos.shape[0]
        N = natoms * 3
        log = ' %s %s= %d\n' % (self.__sec_names['pos'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(natoms):
            for j in xrange(3):
                log+= "%20.10E" % self.__pos[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_atno(self,file):
        """write atomic numbers"""
        natoms = self.__atno.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['atno'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(natoms):
            log+= "%20d" % self.__atno[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_atms(self,file):
        """write atomic masses"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['atms'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(natoms):
            log+= "%20.10E" % self.__atms[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_origin(self,file):
        """write DMTP origin coordinates"""
        nsites = self.__origin.shape[0]
        N = nsites * 3
        log = ' %s %s= %d\n' % (self.__sec_names['origin'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nsites):
            for j in xrange(3):
                log+= "%20.10E" % self.__origin[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_esp(self,file):
        """write ESP charges"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['esp'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(natoms):
            log+= "%20.10E" % self.__esp[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_chlpg(self,file):
        """write ESP charges"""
        natoms = self.__atms.shape[0]
        N = natoms
        log = ' %s %s= %d\n' % (self.__sec_names['chlpg'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(natoms):
            log+= "%20.10E" % self.__chlpg[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_rdma(self,file):
        """write DMTP center coordinates"""
        ndma = self.__rdma.shape[0]
        N = ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['rdma'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(3):
                log+= "%20.10E" % self.__rdma[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_dmac(self,file):
        """write DMTP charges"""
        ndma = self.__dmac.shape[0]
        N = ndma
        log = ' %s %s= %d\n' % (self.__sec_names['dmac'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            log+= "%20.10E" % self.__dmac[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmac1(self,file):
        """write DMTP charges 1st der wrt nmodes"""
        nmodes, ndma = self.__dmac1.shape
        N = ndma * nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['dmac1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ndma):
                log+= "%20.10E" % self.__dmac1[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmac2(self,file):
        """write DMTP charges 2nd derivatives wrt <mode>"""
        ndma = self.__dmac2.shape[0]
        N = ndma
        log = ' %s %s= %d\n' % (self.__sec_names['dmac2'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            log+= "%20.10E" % self.__dmac2[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad(self,file):
        """write DMTP dipoles"""
        ndma = self.__dmad.shape[0]
        N = ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(3):
                log+= "%20.10E" % self.__dmad[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad1(self,file):
        """write DMTP dipoles 1st der wrt nmodes"""
        nmodes,ndma,c = self.__dmad1.shape
        N = nmodes * ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ndma):
                for k in xrange(3):
                    log+= "%20.10E" % self.__dmad1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmad2(self,file):
        """write DMTP dipoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmad2.shape[0]
        N = ndma * 3
        log = ' %s %s= %d\n' % (self.__sec_names['dmad2'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(3):
                log+= "%20.10E" % self.__dmad2[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq(self,file):
        """write DMTP quadrupoles"""
        ndma = self.__dmaq.shape[0]
        N = ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(6):
                log+= "%20.10E" % self.__dmaq[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq1(self,file):
        """write DMTP quadrupoles 1st der wrt nmodes"""
        nmodes,ndma,c = self.__dmaq1.shape
        N = nmodes * ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ndma):
                for k in xrange(6):
                    log+= "%20.10E" % self.__dmaq1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmaq2(self,file):
        """write DMTP quadrupoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmaq2.shape[0]
        N = ndma * 6
        log = ' %s %s= %d\n' % (self.__sec_names['dmaq2'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(6):
                log+= "%20.10E" % self.__dmaq2[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao(self,file):
        """write DMTP octupoles"""
        ndma = self.__dmao.shape[0]
        N = ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(10):
                log+= "%20.10E" % self.__dmao[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao1(self,file):
        """write DMTP octupoles 1st der wrt nmodes"""
        nmodes,ndma,c = self.__dmao1.shape
        N = nmodes * ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ndma):
                for k in xrange(10):
                    log+= "%20.10E" % self.__dmao1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_dmao2(self,file):
        """write DMTP octupoles 2nd derivatives wrt <mode>"""
        ndma = self.__dmao2.shape[0]
        N = ndma * 10
        log = ' %s %s= %d\n' % (self.__sec_names['dmao2'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ndma):
            for j in xrange(10):
                log+= "%20.10E" % self.__dmao2[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_rpol(self,file):
        """write polarizable center coordinates"""
        npol = self.__rpol.shape[0]
        N = npol * 3
        log = ' %s %s= %d\n' % (self.__sec_names['rpol'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(npol):
            for j in xrange(3):
                log+= "%20.10E" % self.__rpol[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_dpol(self,file):
        """write distributed polarizabilities"""
        nmos = self.__dpol.shape[0]
        N = nmos*9
        log = ' %s %s= %d\n' % (self.__sec_names['dpol'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmos):
            for j in [0,1,2]:
                for k in [0,1,2]:
                    log+= "%20.10E" % self.__dpol[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_dpol1(self,file):
        """write distributed polarizabilities first derivatives wrt modes"""
        nmodes = self.__dpol1.shape[0]
        nmos   = self.__dpol1.shape[1]
        N = nmodes*nmos*9
        log = ' %s %s= %d\n' % (self.__sec_names['dpol1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(nmos):
                for k in [0,1,2]:
                    for l in [0,1,2]:
                        log+= "%20.10E" % self.__dpol1[i,j,k,l]
                        if not n%5: log+= '\n'
                        n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_freq(self,file):
        """write harmonic frequencies"""
        nmodes = self.__freq.shape[0]
        N = nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['freq'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            log+= "%20.10E" % self.__freq[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_redmass(self,file):
        """write reduced masses"""
        nmodes = self.__redmass.shape[0]
        N = nmodes
        log = ' %s %s= %d\n' % (self.__sec_names['redmass'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            log+= "%20.10E" % self.__redmass[i]
            if not n%5: log+= '\n'
            n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_lvec(self,file):
        """write mass-weighted eigenvectors"""
        nmodes, natoms, x = self.__lvec.shape 
        N = nmodes * natoms * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lvec'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(natoms):
                for k in [0,1,2]:
                    log+= "%20.10E" % self.__lvec[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_gijk(self,file):
        """write cubic anharmonic constants"""
        nmodes = self.__gijk.shape[0]
        N = nmodes*(nmodes+1)*(nmodes+2)/6
        log = ' %s %s= %d\n' % (self.__sec_names['gijk'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(i+1):
                for k in xrange(j+1):
                    log+= "%20.10E" % self.__gijk[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_lmoc(self,file):
        """write LMO centroid coordinates"""
        nmos = self.__lmoc.shape[0]
        N = nmos * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lmoc'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmos):
            for j in xrange(3):
                log+= "%20.10E" % self.__lmoc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_lmoc1(self,file):
        """write LMO centroid first derivatives wrt nmodes"""
        nmodes, nmos, n = self.__lmoc1.shape
        N = nmodes * nmos * 3
        log = ' %s %s= %d\n' % (self.__sec_names['lmoc1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(nmos):
                for k in xrange(3):
                    log+= "%20.10E" % self.__lmoc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_fock(self,file):
        """write Fock matrix elements"""
        nmos = self.__fock.shape[0]
        N = (nmos**2 - nmos) / 2 + nmos
        log = ' %s %s= %d\n' % (self.__sec_names['fock'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmos):
            for j in xrange(i+1):
                log+= "%20.10E" % self.__fock[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_fckc(self,file):
        """write canonical Fock matrix elements"""
        ncmos = self.__fckc.shape[0]
        N = (ncmos**2 - ncmos) / 2 + ncmos
        log = ' %s %s= %d\n' % (self.__sec_names['fckc'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ncmos):
            for j in xrange(i+1):
                log+= "%20.10E" % self.__fckc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_fock1(self,file):
        """write Fock matrix element first derivatives wrt nmodes"""
        nmodes = self.__fock1.shape[0]
        nmos = self.__fock1.shape[1]
        N = nmodes * ( (nmos**2 - nmos) / 2 + nmos )
        log = ' %s %s= %d\n' % (self.__sec_names['fock1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(nmos):
                for k in xrange(j+1):
                    log+= "%20.10E" % self.__fock1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_fckc1(self,file):
        """write canonical Fock matrix element first derivatives wrt nmodes"""
        nmodes = self.__fckc1.shape[0]
        ncmos = self.__fckc1.shape[1]
        N = nmodes * ( (ncmos**2 - ncmos) / 2 + ncmos )
        log = ' %s %s= %d\n' % (self.__sec_names['fckc1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ncmos):
                for k in xrange(j+1):
                    log+= "%20.10E" % self.__fckc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
        
    def _write_vecl(self,file):
        """write AO-LMO transformation matrix elements"""
        nmos, nbasis = self.__vecl.shape
        N = nmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecl'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmos):
            for j in xrange(nbasis):
                log+= "%20.10E" % self.__vecl[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_vecc(self,file):
        """write AO-CMO transformation matrix elements"""
        ncmos, nbasis = self.__vecc.shape
        N = ncmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecc'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(ncmos):
            for j in xrange(nbasis):
                log+= "%20.10E" % self.__vecc[i,j]
                if not n%5: log+= '\n'
                n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
    def _write_vecl1(self,file):
        """write AO-LMO transformation matrix element first derivatives wrt nmodes"""
        nmodes, nmos, nbasis = self.__vecl1.shape
        N = nmodes * nmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecl1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(nmos):
                for k in xrange(nbasis):
                    log+= "%20.10E" % self.__vecl1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return

    def _write_vecc1(self,file):
        """write AO-CMO transformation matrix element first derivatives wrt nmodes"""
        nmodes, ncmos, nbasis = self.__vecc1.shape
        N = nmodes * ncmos * nbasis
        log = ' %s %s= %d\n' % (self.__sec_names['vecc1'].ljust(40),'N'.rjust(10),N)
        n = 1
        for i in xrange(nmodes):
            for j in xrange(ncmos):
                for k in xrange(nbasis):
                    log+= "%20.10E" % self.__vecc1[i,j,k]
                    if not n%5: log+= '\n'
                    n+=1
        log+= '\n'
        if N%5: log+= '\n'
        file.write(log)
        return
    
