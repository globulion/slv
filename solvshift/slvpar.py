# --------------------------------------------------------------------- #
#             SOLVATOCHROMIC PARAMETER FORMAT MODULE                    #
# --------------------------------------------------------------------- #

from numpy     import array, float64, zeros, newaxis, sqrt, dot,asfortranarray, transpose
from units     import *
from dma       import DMA
from utilities import order, SVDSuperimposer as svd_sup, MakeMol
from PyQuante.Ints import getbasis
from efprot    import vecrot, vc1rot
import sys, copy, os, re, math
sys.stdout.flush()

__all__ = ['SLVPAR',]
__version__ = '1.0.1'

class SLVPAR(object,UNITS):
    """
Represents Solvshift Solvatochromic Parameter Format System
-----------------------------------------------------------

It contains various sections necessary for evaluation of
frequency shifts using coarse-grained models. The major gro
ups of parameters are:

0) basic molecule specification
   a) the structure - atomic coordinates
   b) atomic numbers
   c) atomic masses
1) electrostatics - distributed multipole approximation
   a) first-order Coulomb forces
   b) second-order polarization forces
2) frequency analysis data
3) non-electrostatics:
   a) exchange-repulsion
   b) charge-transfer
Part 1a) is a part of Coulomb file format.

Beneath I list the various sections of both groups of parameters:

 - GROUP (0) -
 [ molecule ] - basic specifications of the molecule
 [ structure ]
 [ atomic numbers ]
 [ atomic masses ]
 
 - GROUP (1) -
 [ DMTP ] - distributed multipole moments
 [ SolDMTP ] - solvatochromic  moments
 [ origins ] - origins of distributed moments
 
 - GROUP (3) -
 [ Harmonic frequencies ]
 [ Reduced masses ]
 [ Mass-weighted eigenvectors ]
 [ Cubic anharmonic constants ]
 
 - GROUP (2) -
 [ AO to LMO transformation matrix ] = [ AO->LMO matrix ]
 [ LMO centroids ]
 [ Fock matrix ] - written in LMO basis
 [ AO to LMO transformation matrix - first derivatives ] = [ AO->LMO matrix - first derivatives ]
 [ LMO centroids - first derivatives ]
 [ Fock matrix - first derivatives ]
 
Usage:

Notes:
1) all the data in the instance of SLVPAR are stored
   in atomic units!!!
"""
    def __init__(self,file=None):
        self.__file = file
        self._create()
        self._nlines = lambda n: n/5+bool(n%5)
        if file is not None: self.__call__(file)
        
    def __call__(self,file):
        """parse the parameters from a file"""
        self.__file = file
        filef = open(file,'r')
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
            dma=None,chelpg=None,esp=None):
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
           assert dma.if_traceless(), "DMA object is NOT traceless!!!"
           self.__pos = dma.get_pos()
           self.__origin = dma.get_origin()
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
        # other
        self.__chlpg  = chelpg
        self.__esp    = esp
        return
        
    def get(self):
        """returns dictionary with parameters"""
        return self._make_dict()
    
    def get_bfs(self):
        """return basis set object"""
        mol = MakeMol(self.__atno,self.__pos,name=self.__name)
        bfs = getbasis(mol,self.__basis)
        return bfs
    
    def get_pos(self):
        """return position of a fragment"""
        return self.__pos
    
    def write(self,file='slv.par',par=None):
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
        if 0: pass
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
    
    def sup(self,str):
        """superimpose structures <str> and <self.__pos>. Return rms in A.U."""
        s = svd_sup()
        s.set(str,self.__pos)
        s.run()
        rms = s.get_rms()
        rot, transl = s.get_rotran()
        # perform transformations
        self.__pos    = s.get_transformed()
        if self.__lmoc  is not None: self.__lmoc   = dot(self.__lmoc , rot) + transl
        if self.__rpol  is not None: self.__rdpol  = dot(self.__rpol , rot) + transl
        if self.__rdma  is not None: self.__rdma   = dot(self.__rdma , rot) + transl
        if self.__lmoc1 is not None: self.__lmoc1  = dot(self.__lmoc1, rot)
        if self.__lvec  is not None: self.__lvec   = dot(self.__lvec , rot)
        # transform dipoles, quadrupoles and octupoles!
        if 0: pass
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
    
    # protected
    
    def _create(self):
        """creates the namespace of parameter variables"""
        self.__name ,self.__natoms, self.__nbasis = None, None, None
        self.__nmos ,self.__nmodes, self.__basis  = None, None, None
        self.__atoms,self.__shortname             = None, None
        self.__pos  ,self.__origin, self.__nsites = None, None, None
        self.__atno ,self.__atms                  = None, None
        #
        self.__redmass, self.__freq, self.__lvec = None, None, None
        self.__gijk                              = None
        #
        self.__fock ,self.__lmoc ,self.__vecl = None, None, None
        self.__fock1,self.__lmoc1,self.__vecl1= None, None, None
        self.__ncmos,self.__vecc ,self.__vecc1= None, None, None
        self.__fckc ,self.__fckc1             = None, None
        #
        self.__esp  ,self.__chlpg             = None, None
        self.__dpol ,self.__dpol1,self.__rpol = None, None, None
        #
        self.__rdma, self.__dmac, self.__dmad = None, None, None
        self.__dmaq, self.__dmao              = None, None
        #
        mol_names = ('name','basis','method','natoms','nbasis',
                     'nmos','nmodes','atoms','shortname','nsites',
                     'ncmos','npol','ndma',)
        sec_names = {'lmoc' : '[ LMO centroids ]',
                     'lmoc1': '[ LMO centroids - first derivatives ]',
                     'fock' : '[ Fock matrix ]',
                     'fock1': '[ Fock matrix - first derivatives ]',
                     'fckc' : '[ Canonical Fock matrix ]',
                     'fckc1': '[ Canonical Fock matrix - first derivatives ]',
                     'vecl' : '[ AO->LMO matrix ]',
                     'vecl1': '[ AO->LMO matrix - first derivatives ]',
                     'vecc' : '[ AO->CMO matrix ]',
                     'vecc1': '[ AO->CMO matrix - first derivatives ]',
                     'mol'  : '[ molecule ]',
                     'freq' : '[ Harmonic frequencies ]',
                   'redmass': '[ Reduced masses ]',
                     'lvec' : '[ Mass-weighted eigenvectors ]',
                     'gijk' : '[ Cubic anharmonic constants ]',
                      'pos' : '[ Atomic coordinates ]',
                   'origin' : '[ DMTP origins ]',
                     'atno' : '[ Atomic numbers ]',
                     'atms' : '[ Atomic masses ]',
                     'esp'  : '[ ESP charges ]',
                     'chlpg': '[ ChelpG charges ]',
                      'dpol': '[ Distributed polarizabilities ]',
                     'dpol1': '[ Distributed polarizabilities - first derivatives ]',
                      'rpol': '[ Polarizable centers ]',
                      'rdma': '[ DMTP centers ]',
                      'dmac': '[ DMTP charges ]',
                      'dmad': '[ DMTP dipoles ]',
                      'dmaq': '[ DMTP quadrupoles ]',
                      'dmao': '[ DMTP octupoles ]',}
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
        if self.__rpol  is not None: par['rpol' ] = self.__rpol
        if self.__dpol  is not None: par['dpol' ] = self.__dpol
        if self.__dpol1 is not None: par['dpol1'] = self.__dpol1
        #
        return par
        
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
        log = ' %s\n' % self.__sec_names['mol'].ljust(40)
        log+= '   name       = %s\n'    % self.__name
        log+= '   basis      = %s/%s\n' % (self.__method, self.__basis)
        log+= '   natoms     = %s\n'    % self.__natoms
        log+= '   nbasis     = %s\n'    % self.__nbasis
        log+= '   nmodes     = %s\n'    % self.__nmodes
        log+= '   nmos       = %s\n'    % self.__nmos
        log+= '   ncmos      = %s\n'    % self.__ncmos
        log+= ' \n'
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
    