# --------------------------------------------------------------------- #
#             SOLVATOCHROMIC PARAMETER FORMAT MODULE                    #
# --------------------------------------------------------------------- #

from numpy     import array, float64, zeros
from units     import *
from dma       import DMA
from utilities import order
import sys, copy, os, re, math
sys.stdout.flush()

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
2) non-electrostatics:
   a) exchange-repulsion
   b) charge-transfer
Part 1a) is a part of Coulomb file format.

Beneath I list the various sections of both groups of param
eters:

 - GROUP (0) -
 [ molecule ] - basic specifications of the molecule
 [ structure ]
 [ atomic numbers ]
 [ atomic masses ]
 
 - GROUP (1) -
 [ DMTP ] - distributed multipole moments
 [ SolDMTP ] - solvatochromic  moments
 [ origins ] - origins of distributed moments
 
 - GROUP (2) -
 [ AO to LMO transformation matrix ] = [ AO->LMO matrix ]
 [ LMO centroids ]
 [ Fock matrix ] - written in LMO basis
 [ AO to LMO transformation matrix - first derivatives ] = [ AO->LMO matrix - first derivatives ]
 [ LMO centroids - first derivatives ]
 [ Fock matrix - first derivatives ]
"""
    def __init__(self,file=None):
        self.__file = file
        self._create()
        self._nlines = lambda n: n/5+bool(n%5)
        if file is not None: self.__call__(file)
        
    def __call__(self,file):
        """parse the parameters from a file"""
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
        log = '\n'
        log+= ' ================================= \n'
        log+= ' SLV SOLVATOCHROMIC EFP PARAMETERS \n'
        log+= ' ================================= \n\n'
        log+= ' MOLECULE SPECIFICATION \n'
        log+= ' name   %s \n' % self.__name
        log+= ' method %s \n' % self.__basis
        log+= ' natoms %s \n' % self.__natoms
        log+= ' nbasis %s \n' % self.__nbasis
        log+= ' nmos   %s \n' % self.__nmos
        log+= ' nmodes %s \n' % self.__nmodes
        log+= '\n'
        return str(log)
    
    # public
    
    def set(self,fock=None,lmoc=None,vecl=None,
                 fock1=None,lmoc1=None,vecl1=None):
        """set the parameters"""
        self.__fock = fock; self.__lmoc = lmoc; self.__vecl = vecl
        self.__fock1 = fock1; self.__lmoc1 = lmoc1; self.__vecl1 = vecl1
        return

    def get(self):
        """returns dictionary with parameters"""
        par = {}
        if self.__lmoc  is not None: par['lmoc' ] = self.__lmoc
        if self.__lmoc1 is not None: par['lmoc1'] = self.__lmoc1
        if self.__fock  is not None: par['fock' ] = self.__fock
        if self.__fock1 is not None: par['fock1'] = self.__fock1
        if self.__vecl  is not None: par['vecl' ] = self.__vecl
        if self.__vecl1 is not None: par['vecl1'] = self.__vecl1
        return par
        
    def write(self,file='slv.par'):
        """writes the parameters into a file"""
        f = open(file,'w')
        if self.__lmoc  is not None: self._write_lmoc(f)
        if self.__lmoc1 is not None: self._write_lmoc1(f)
        if self.__fock  is not None: self._write_fock(f)
        if self.__fock1 is not None: self._write_fock1(f)
        if self.__vecl  is not None: self._write_vecl(f)
        if self.__vecl1 is not None: self._write_vecl1(f)
        f.close()
        return
    
    # protected
    
    def _create(self):
        """creates the namespace of parameter variables"""
        self.__name, self.__natoms, self.__nbasis = None, None, None
        self.__nmos, self.__nmodes, self.__basis  = None, None, None
        #
        self.__fock, self.__lmoc, self.__vecl = None, None, None
        self.__fock1,self.__lmoc1,self.__vecl1= None, None, None
        #
        mol_names = ('name','basis','method','natoms','nbasis',
                     'nmos','nmodes',)
        sec_names = {'lmoc' : '[ LMO centroids ]',
                     'lmoc1': '[ LMO centroids - first derivatives ]',
                     'fock' : '[ Fock matrix ]',
                     'fock1': '[ Fock matrix - first derivatives ]',
                     'vecl' : '[ AO->LMO matrix ]',
                     'vecl1': '[ LMO centroids - first derivatives ]',
                     'mol'  : '[ molecule ]',}
        self.__mol_names = mol_names
        self.__sec_names = sec_names
        return
        
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
                    if name == 'basis':
                        self.__basis = arg
                    if name == 'natoms':
                        self.__natoms = int(arg)
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
            # LMO centroids
            if   key == 'lmoc':
                merror = 'nmos in section [ molecule ] '
                merror+= 'is not consistent with section [ LMO centroids ]!'
                assert self.__nmos == N/3, merror
                data = data.reshape(self.__nmos,3)
                self.__lmoc = data
            # LMO centroids first derivatives
            elif key == 'lmoc1':
                merror = 'nmos and nmodes in section [ molecule ] '
                merror+= 'is not consistent with section [ LMO centroids - first derivatives ]!'
                assert self.__nmos == N/(3*self.__nmodes), merror
                data = data.reshape(self.__nmodes,self.__nmos,3)
                self.__lmoc1 = data
            # Fock matrix
            elif key == 'fock':
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
            elif key == 'fock1':
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
                self.__fock1 = fock1
            # AO to LMO transformation matrix
            elif key == 'vecl':
                merror = 'nmos and nbasis in section [ molecule ] '
                merror+= 'is not consistent with section [ AO->LMO matrix ]!'
                assert self.__nmos == N/self.__nbasis, merror
                data = data.reshape(self.__nmos,self.__nbasis)
                self.__vecl = data
            # AO to LMO transformation matrix first derivatives
            elif key == 'vecl1':
                merror = 'nmodes, nmos and nbasis in section [ molecule ] '
                merror+= 'is not consistent with section [ AO->LMO matrix - first derivatives ]!'
                assert self.__nmos == N/(self.__nbasis*self.__nmodes)
                data = data.reshape(self.__nmodes,self.__nmos,self.__nbasis)
                self.__vecl1 = data
            # other not-programmed section
            else:
                raise Exception('Non-standard section name detected! Check the format of your file!')
        return
    
    def _write_lmoc(self,file):
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
        file.write(log)
        return
    
    def _write_lmoc1(self,file):
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
        file.write(log)
        return
        
    def _write_fock(self,file):
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
        file.write(log)
        return
    
    def _write_fock1(self,file):
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
        file.write(log)
        return
    
    def _write_vecl(self,file):
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
        file.write(log)
        return
    
    def _write_vecl1(self,file):
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
        file.write(log)
        return
    