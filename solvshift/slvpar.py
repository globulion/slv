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
    
    def set(self,mol=None,anh=None,dma=None,frag=None):
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
           self.__atno   = mol.get_atno()
           self.__atms   = mol.get_atms() * self.AmuToElectronMass
        # electrostatic data
        if dma is not None:
           self.__pos = dma.get_pos()
           self.__origin = dma.get_origin()
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
        f.close()
        return
    
    def rotate(self,rot):
        """rotate the tensors by <rot> unitary matrix"""
        # transform the atomic position static and dynamic information
        if self.__pos   is not None:
           self.__pos    = dot(self.__pos, rot)
        if self.__lmoc  is not None:
           self.__lmoc   = dot(self.__lmoc, rot)
        if self.__lmoc1 is not None:
           self.__lmoc1  = dot(self.__lmoc1, rot)
        if self.__lvec  is not None:
           self.__lvec   = dot(self.__lvec, rot)
        # transform the wave function!
        if self.__vecl  is not None:
           bfs = self.get_bfs()
           typs= bfs.get_bfst().sum(axis=1)
           self.__vecl = vecrot(self.__vecl,rot,typs)
           if self.__vecl1 is not None:
              self.__vecl1 = vc1rot(self.__vecl1,rot,typs)
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
        if self.__lmoc1 is not None: self.__lmoc1  = dot(self.__lmoc1, rot)
        if self.__lvec  is not None: self.__lvec   = dot(self.__lvec , rot)
        #
        if self.__vecl  is not None:
           bfs = self.get_bfs()
           typs= bfs.get_bfst().sum(axis=1)
           self.__vecl = vecrot(self.__vecl,rot,typs)
           if self.__vecl1 is not None:
              self.__vecl1 = vc1rot(self.__vecl1,rot,typs)     
        return rms
    
    def copy(self):
        """return a deepcopy of me!"""
        return copy.deepcopy(self)
    
    # protected
    
    def _create(self):
        """creates the namespace of parameter variables"""
        self.__name, self.__natoms, self.__nbasis = None, None, None
        self.__nmos, self.__nmodes, self.__basis  = None, None, None
        self.__atoms,self.__shortname = None, None
        self.__pos, self.__origin, self.__nsites = None, None, None
        self.__atno, self.__atms = None, None
        #
        self.__redmass, self.__freq, self.__lvec = None, None, None
        self.__gijk = None
        #
        self.__fock, self.__lmoc, self.__vecl = None, None, None
        self.__fock1,self.__lmoc1,self.__vecl1= None, None, None
        #
        mol_names = ('name','basis','method','natoms','nbasis',
                     'nmos','nmodes','atoms','shortname','nsites',)
        sec_names = {'lmoc' : '[ LMO centroids ]',
                     'lmoc1': '[ LMO centroids - first derivatives ]',
                     'fock' : '[ Fock matrix ]',
                     'fock1': '[ Fock matrix - first derivatives ]',
                     'vecl' : '[ AO->LMO matrix ]',
                     'vecl1': '[ AO->LMO matrix - first derivatives ]',
                     'mol'  : '[ molecule ]',
                     'freq' : '[ Harmonic frequencies ]',
                   'redmass': '[ Reduced masses ]',
                     'lvec' : '[ Mass-weighted eigenvectors ]',
                     'gijk' : '[ Cubic anharmonic constants ]',
                      'pos' : '[ Atomic coordinates ]',
                   'origin' : '[ DMTP origins ]',
                     'atno' : '[ Atomic numbers ]',
                     'atms' : '[ Atomic masses ]',}
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
        return par
        
    def _tr_lvec(self,lvec,nmodes,natoms):
        """transpose axis and reshape"""
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
                    if name == 'nmodes':
                        self.__nmodes = int(arg)
                
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
                 data = data
                 self.__atno = array(data,int)
            # Atomic masses
            elif key == 'atms':
                 merror = 'natoms in section [ molecule ] '
                 merror+= 'is not consistent with section [ Atomic masses ]!'
                 assert self.__natoms == N, merror
                 data = data
                 self.__atms = data
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
            # AO to LMO transformation matrix
            elif  key == 'vecl':
                  merror = 'nmos and nbasis in section [ molecule ] '
                  merror+= 'is not consistent with section [ AO->LMO matrix ]!'
                  assert self.__nmos == N/self.__nbasis, merror
                  data = data.reshape(self.__nmos,self.__nbasis)
                  self.__vecl = data
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
    