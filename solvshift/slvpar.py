# --------------------------------------------------------------------- #
#             SOLVATOCHROMIC PARAMETER FORMAT MODULE                    #
# --------------------------------------------------------------------- #

from numpy     import array
from units     import *
from dma       import DMA
from utilities import order
import sys, copy, os
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
        self.__fock, self.__lmoc, self.__vecl = None, None, None
        self.__fock1,self.__lmoc1,self.__vecl1= None, None, None
    
    # public
    
    def set(self,fock=None,lmoc=None,vecl=None,
            fock1=None,lmoc1=None,vecl1=None):
        """set the parameters"""
        self.__fock = fock; self.__lmoc = lmoc; self.__vecl = vecl
        self.__fock1 = fock1; self.__lmoc1 = lmoc1; self.__vecl1 = vecl1
        return
    
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
    
    def _write_fock(self,file):
        nmos = self.__fock.shape[0]
        N = (nmos**2 - nmos) / 2 + nmos
        log = ' %s %s= %d\n' % ('[ Fock matrix ]'.ljust(40),'N'.rjust(10),N)
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
        log = ' %s %s= %d\n' % ('[ Fock matrix - first derivatives ]'.ljust(40),'N'.rjust(10),N)
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
    def _write_lmoc(self,file):
        nmos = self.__lmoc.shape[0]
        N = nmos * 3
        log = ' %s %s= %d\n' % ('[ LMO centroids ]'.ljust(40),'N'.rjust(10),N)
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
        log = ' %s %s= %d\n' % ('[ LMO centroids - first derivatives ]'.ljust(40),'N'.rjust(10),N)
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
    def _write_vecl(self,file):
        nmos, nbasis = self.__vecl.shape
        N = nmos * nbasis
        log = ' %s %s= %d\n' % ('[ AO->LMO matrix ]'.ljust(40),'N'.rjust(10),N)
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
        log = ' %s %s= %d\n' % ('[ AO->LMO matrix - first derivatives ]'.ljust(40),'N'.rjust(10),N)
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
    