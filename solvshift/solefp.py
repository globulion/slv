# -------------------------------------------------------------------------- #
#            SOLVATOCHROMIC EFFECTIVE FRAGMENT POTENTIAL MODULE              #
# -------------------------------------------------------------------------- #

from numpy     import tensordot, dot, transpose, array, float64
from units     import *
from utilities import Read_xyz_file, get_pmloca, ParseFockFromGamessLog, \
                      ParseDmatFromFchk, ParseVecFromFchk
from diff      import DIFF
from PyQuante.Ints import getSAB, getTAB, getSA1B, getTA1B
from shftex import shftex
import sys, copy, os, re, math, glob, PyQuante.Ints, coulomb.multip
sys.stdout.flush()

__all__ = ['Frag','EFP',]
__version__ = '1.0.1'

class EFP(object,UNITS):
    """
=============================================================================
                     EFFECTIVE FRAGMENT POTENTIAL METHOD                     
=============================================================================

Usage:
A = EFP(a,b)
result = A(shift=False,nmode=None,cunit=False)
rms_a, rms_b = A.sup(str_a=None,str_b=None)

Notes:
1) a and b are SLVPAR instances. If shift=True, a denotes IR-active molecule
   and should contain all necessary parameters for it
2) the a and b objects are assumed to be appropriately transformed in space 
   by rotations and translations
3) a and b in sup argument list are ndarray structures of dimension (natoms,3)
   The coordinates are assumed to be in A.U.
4) cunit - change units. If <True> then energies are returned in [kcal/mole]
   and frequency shifts in [cm-1]. Otherwise all is returned in A.U. units
"""
    def __init__(self,a,b):
        self.__molA = a
        self.__molB = b
        return
    
    def __call__(self,shift=False,nmode=None,cunit=False):
        """perform all the calculation"""
        ma, ea = self._ex(shift,nmode,cunit)
        return ma, ea
    
    def __repr__(self):
        """print the two fragments"""
        log = '\n'
        log+= str(self.__molA.get_pos()*self.BohrToAngstrom)
        log+= '\n'
        log+= str(self.__molB.get_pos()*self.BohrToAngstrom)
        log+= '\n'
        return str(log)
    
    def sup(self,str_a=None,str_b=None):
        """superimpose the a and b objects to given structures"""
        rms_a, rms_b = None, None
        if str_a is not None:
           rms_a = self.__molA.sup(str_a)
        if str_b is not None:
           rms_b = self.__molB.sup(str_b)
        return rms_a, rms_b
    
    # protected
    
    def _ex(self,shift,nmode,cunit):
        """calculate exchange-repulsion property"""
        # variables
        varA = self.__molA.get()
        varB = self.__molB.get()
        # basis sets
        bfsA = self.__molA.get_bfs()
        bfsB = self.__molB.get_bfs()
        # instantaneous integrals
        skm  = getSAB(bfsA,bfsB)
        tkm  = getTAB(bfsA,bfsB)
        sk1m = getSA1B(bfsA,bfsB)
        tk1m = getTA1B(bfsA,bfsB)
        # parameters
        ### molecule A
        faij   = varA['fock']
        faij1  = varA['fock1']
        cika   = varA['vecl']
        cika1  = varA['vecl1']
        za     = varA['atno']
        rna    = varA['pos']
        ria    = varA['lmoc']
        ria1   = varA['lmoc1']
        redmss = varA['redmass']
        freq   = varA['freq']
        gijj   = varA['gijk'][nmode-1,nmode-1,:]
        lvec   = varA['lvec']
        mlist  = bfsA.get_bfsl() + 1
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # transform the integrals
        sij = dot(dot(cika,skm),transpose(cikb))
        tij = dot(dot(cika,tkm),transpose(cikb))
        # calculate the properties!
        shftma,shftea = shftex(redmss,freq,gijj,lvec,
                               ria,rib,rna,rnb,ria1,
                               cika,cikb,cika1,
                               skm,tkm,sk1m,tk1m,
                               za,zb,mlist,
                               faij,fbij,faij1,nmode)
        if cunit:
           shftma *= self.HartreePerHbarToCmRec
           shftea *= self.HartreePerHbarToCmRec

        return shftma, shftea
    
class Frag(object,DIFF):
    """
Solvatochromic Effective Fragment Potential Fragment
----------------------------------------------------

Usage:
"""
    def __init__(self,anh=None,basis=None,nae=None,
                      fchk=None,gmslog=None,):
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
        self._init()
        self._create()

    # public
    
    def set(self,anh=None,basis=None,nae=None,
                 fchk=None,gmslog=None,):
        """set the properties to the object"""
        if self.__anh    is not None: self.__anh    = anh
        if self.__fchk   is not None: self.__fchk   = fchk
        if self.__gmslog is not None: self.__gmslog = gmslog
        if self.__basis  is not None: self.__basis  = basis
        if self.__nae    is not None: self.__nae    = nae
        return

    def reset(self,anh=None,basis=None,nae=None,
                   fchk=None,gmslog=None,):
        """reset the properties"""
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
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
    
    def eval(self):
        """Parses AO-LMO transformation matrix and Fock matrix.
Transforms the latter from AO to LMO space. Computes also 
overlap integrals and parses density matrix."""
        assert self.__mol is not None, 'molecule not specified! (no fchk file)'
        # evaluate transformation matrices and LMO centroids
        SAO   = PyQuante.Ints.getS(self.__bfs)
        dmat = ParseDmatFromFchk(self.__fchk,self.__basis_size)
        veccmo= ParseVecFromFchk(self.__fchk)[:self.__nae,:]
        tran, veclmo = get_pmloca(self.__natoms,mapi=self.__bfs.LIST1,sao=SAO,
                                  vecin=veccmo,nae=self.__nae,
                                  maxit=100000,conv=1.0E-14,
                                  lprint=True,
                                  freeze=None)
        # calculate LMTPs
        camm = coulomb.multip.MULTIP(molecule=self.__mol,
                                     basis=self.__basis,
                                     method='b3lyp',
                                     matrix=dmat,
                                     transition=False,
                                     bonds=None,vec=veclmo)
        camm.camms()
        dma = camm.get()[0]
        # parse Fock matrix
        fock = ParseFockFromGamessLog(self.__gmslog,interpol=False)
        fock = tensordot(veclmo,tensordot(veclmo,fock,(1,0)),(1,1))
        # save
        self.__lmoc = dma.get_origin()[self.__natoms:]
        self.__tran = tran
        self.__vecl = veclmo
        self.__fock = fock
        self.__sao  = SAO
        self.__dmat = dmat
        self.__dma  = dma
        return
    
    # protected
    
    def _init(self):
        """initialize the other memorials"""
        self.__mol    = None; self.__basis_size = None
        self.__lmoc   = None; self.__lmoc1  = None; self.__natoms=None
        self.__fock   = None; self.__fock1  = None; self.__bfs  = None
        self.__vecl   = None; self.__vecl1  = None
        return
    
    def _create(self):
        """creates the molecule"""
        if self.__fchk is not None:
           mol  = Read_xyz_file(self.__fchk,mol=True,
                                mult=1,charge=0,
                                name='happy dummy molecule')
           bfs        = PyQuante.Ints.getbasis(mol,self.__basis)
           basis_size = len(bfs)
           natoms= len(mol.atoms)
           # save
           self.__mol = mol
           self.__bfs = bfs
           self.__basis_size = basis_size
           self.__natoms = natoms
        return
    