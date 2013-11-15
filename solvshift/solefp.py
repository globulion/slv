# -------------------------------------------------------------------------- #
#            SOLVATOCHROMIC EFFECTIVE FRAGMENT POTENTIAL MODULE              #
# -------------------------------------------------------------------------- #

from numpy     import tensordot, dot, transpose, array, float64, zeros, \
                      concatenate as con
from units     import *
from utilities import Read_xyz_file, get_pmloca, ParseFockFromGamessLog, \
                      ParseDmatFromFchk, ParseVecFromFchk, MakeMol
from diff      import DIFF
from PyQuante.Ints import getT, getSAB, getTAB, getSA1B, getTA1B, getVEFP, \
                          getV, getbasis
from shftex import shftex
from shftce import shftce
from cosik  import solpol, mollst
from efprot import tracls
from exrep  import exrep
import sys, copy, os, re, math, glob, PyQuante.Ints, coulomb.multip, clemtp
sys.stdout.flush()

__all__ = ['EFP','FragFactory',]
__version__ = '1.0.1'

class EFP(object,UNITS):
    """"""
    def __init__(self,ccut=None,pcut=None,pairwise_all=False,
                      elect=True,pol=False,rep=False,ct=False,disp=False,all=False,
                      nlo=False,freq=False,mode=None):
        """The global settings for a computation:
ccut         - Coulomb cutoff
pcut         - Polarization cutoff
pairwise_all - count all pairwise interactions (otherwise count
               only these that include central molecule, False 
               when frequency shift mode is switched on in freq and mode)
elect        - evaluate electrostatics
pol          - evaluate polarization
rep          - evaluate exchange-repulsion
ct           - evaluate charge-transfer
disp         - evaluate dispersion
all          - evaluate all these interactions"""
        self._create()
        #
        self.__ccut = ccut
        self.__pcut = pcut
        #
        if all:
           self.__eval_elect = True
           self.__eval_pol   = True
           self.__eval_rep   = True
           self.__eval_ct    = True
           self.__eval_disp  = True
        else:
           self.__eval_elect = elect
           self.__eval_pol   = pol
           self.__eval_rep   = rep
           self.__eval_ct    = ct
           self.__eval_disp  = disp
        #
        if freq: 
           self.__pairwise_all = False
           self.__eval_freq = True
           self.__mode = mode
        else: self.__pairwise_all = pairwise_all
        #
        if nlo: 
           self.__eval_nlo = True
        #
        return
    
    def __call__(self,lwrite=False):
        """evaluate the property"""
        self.eval(lwrite)
        return
    
    # P U B L I C

    def set(self,pos,ind,nmol,bsm=None):
        """Set the initial molecular coordinates and the moltype index list. 
Also set the BSM parameters if not done in set_bsm. 
<pos>     -  an 2D array of atomic coordinates of dimension natoms, 3
<ind>     -  a list of moltype ids in bsm list. Has length equal to the 
          -  total number of molecules
<nmol>    - a list of number of atoms (integers) in each molecule
<bsm>     - a list of a form: bsm = list(<Frag object>)"""
        self.__ind = array(ind,int)
        self.__nmol= array(nmol,int)
        if bsm is not None: self.__bsm = bsm
        self._update(pos)
        return    
    
    def set_cut(self,ccut,pcut):
        """Set the cutoff radii for Coulomb and polarization forces"""
        self.__ccut = ccut
        self.__pcut = pcut
        return
    
    def set_bsm(self,bsm):
        """Set the BSM parameters. <bsm> is a list of a form:
  bsm = list(<Frag object>)"""
        self.__bsm = bsm
        return
        
    def set_pos(self,pos):
        """update position array"""
        self._update(pos)
        return
    
    def eval(self,lwrite=False):
        """evaluate the properties"""
        #                                        #
        #           PAIRWISE ALL MODE            #
        #   (interatction energy, NLO property   #
        #                                        #
        if self.__pairwise_all:
           N = len(self.__nmol)
           PAR = []
           QO  = []
           for i in range(N):
               nm = self.__nmol[i]
               im = self.__ind[i]
               #
               STR = self.__rcoordc[nm*i:nm*(i+1)]
               frg = self.__bsm[im].copy()
               rms = frg.sup( STR )
               par = frg.get()
               PAR.append( par )
               #
               qad, oct = tracls( par['dmaq'], par['dmao'] )
               QO.append( (qad,oct) )
           
           # ----------------------------------- ELECT --------------------------------- #
           if self.__eval_elect:
              ndma = [ x['ndma'] for x in PAR ]
              ndmas= sum(ndma)
              
              rdma = con  ([ x['rdma'] for x in PAR ]).reshape(ndmas*3)
              chg  = con  ([ x['dmac'] for x in PAR ]).reshape(ndmas)
              dip  = con  ([ x['dmad'] for x in PAR ]).reshape(ndmas*3)
              qad  = con  ([ QO[x][0]  for x in range(N)   ]).reshape(ndmas*6)
              oct  = con  ([ QO[x][1]  for x in range(N)   ]).reshape(ndmas*10)
              
              e_el = clemtp.edmtpa(rdma,chg,dip,qad,oct,ndma,lwrite)
              if lwrite: print " Electrostatic      energy: %10.6f"%e_el
          # ------------------------------------- POL --------------------------------- #
              if self.__eval_pol:
                 ### ENERGY
              
                 npol = [ x['npol'] for x in PAR ]
                 npols= sum(npol)
                 #
                 rpol = con  ([ x['rpol'] for x in PAR ]).reshape(npols*3)
                 pol  = con  ([ x['dpol'] for x in PAR ]).reshape(npols*9)
                 #
                 DIM  = npols*3
                 dmat = zeros((DIM,DIM),float64)
                 flds = zeros( DIM,float64)
                 dipind=zeros( DIM,float64)
                 #
                 e_pol = solpol(rdma,chg,dip,qad,oct,
                                rpol,pol,dmat,
                                flds,dipind,
                                ndma,npol,lwrite=False)
                 if lwrite: print " Polarization       energy: %10.6f"%e_pol
              
              ### FREQUENCY SHIFTS
              if self.__eval_freq:
                  pass
              
              ### NLO PROPERTY
              if self.__eval_nlo:
                  pass
           # ---------------------------------- EX-REP --------------------------------- #
           if  self.__eval_rep:
               e_rep = 0.0
               for i in xrange(N):
                   for j in xrange(i):
                       varA = PAR[i]
                       varB = PAR[j]
                       e_rep += self._pair_rep(varA,varB)
                       #print e_rep
               if lwrite: print " Exchange-repulsion energy: %10.6f"%e_rep
                              
        else:
        #                                            #
        #           CENTRAL MOLECULE MODE            #
        #    (interaction energy,frequency shift)    #
        #                                            #
           print " NOT IMPLEMENTED YET. QUITTING..."
        #
        return
    
    def __repr__(self):
        """print the status"""
        log = '\n'
        log+= " STORED MOLECULES\n"
        if self.__bsm is None:
           log+= " --- no molecules added ---"
        else:
           for i,bsm in enumerate(self.__bsm):
               log+= " %5i %15s\n" % (i,bsm.get_name())
        return str(log)
    
    # P R O T E C T E D
    
    def _create(self):
        """namespace of objects"""
        self.__bsm = None
        self.__eval_freq = False
        self.__eval_nlo  = None
        return
    
    def _update(self,pos):
        """update the neighbour/in-sphere lists based on actual <pos> coordinate array"""
        if not self.__pairwise_all:
           nm = len(self.__nmol)-1
           nac= self.__nmol[0]
           natoms = len(pos)
           rc = pos[:nac ]
           rm = pos[ nac:].reshape((natoms-nac)*3)
           #
           ic = zeros(natoms-nac,dtype=bool)
           ip = zeros(natoms-nac,dtype=bool)
           icm = zeros(nm,dtype=bool)
           icp = zeros(nm,dtype=bool)
           #
           nccut, npcut, mccut, mpcut = mollst(rc,rm,ic,ip,icm,icp,
                                               self.__nmol[1:],
                                               self.__ccut,self.__pcut)
                                               #
           # molecule coordinate lists
           rcoordc = pos[nac:][array(nccut,dtype=bool)]
           rcoordp = pos[nac:][array(npcut,dtype=bool)]
           # molecule type lists
           mtc= self.__ind[1:][array(mccut,dtype=bool)]
           mtp= self.__ind[1:][array(mpcut,dtype=bool)]
           # moleculear atom number lists
           ntc= self.__nmol[1:][array(mccut,dtype=bool)]
           ntp= self.__nmol[1:][array(mpcut,dtype=bool)]
           #
           self.__rc = rc
           self.__rcoordc = rcoordc
           self.__rcoordp = rcoordp
           self.__mtc = mtc
           self.__mtp = mtp
           self.__ntc = ntc
           self.__ntp = ntp
        else:
           self.__rcoordc = pos
        return

    # PAIR ENERGIES
    
    def _pair_elect(self,varA,varB):
        """MTP electrostatic pair energy"""
        return
    
    def _pair_rep(self,varA,varB):
        """exchange-repulsion pair energy"""
        # basis sets
        molA = MakeMol(varA['atno'],varA['pos'])
        molB = MakeMol(varB['atno'],varB['pos'])
        bfsA = getbasis(molA,varA['basis'])
        bfsB = getbasis(molB,varB['basis'])
        # instantaneous integrals
        skm  = getSAB(bfsA,bfsB)
        tkm  = getTAB(bfsA,bfsB)
        # parameters
        ### molecule A
        faij = varA['fock']
        cika = varA['vecl']
        za   = varA['atno']
        rna  = varA['pos']
        ria  = varA['lmoc']
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # calculate the properties!
        eint = exrep(ria,rib,rna,rnb,faij,fbij,cika,cikb,skm,tkm,za,zb)
        return eint
        


class EFP_pair(object,UNITS):
    """
=============================================================================
              EFFECTIVE FRAGMENT POTENTIAL METHOD FOR A DIMER                
=============================================================================

Usage:
A = EFP(a,b)
result = A(ct=True,nlo=False,shift=False,nmode=None,cunit=False)
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
    
    def __call__(self,ct=True,nlo=False,shift=False,nmode=None,cunit=False):
        """perform all the calculation"""
        if ct: 
             ma, ea = self._ex_rep_ct(nlo,shift,nmode,cunit)
        else:
             ma, ea = self._ex_rep(nlo,shift,nmode,cunit)
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
        if str_a is not None: rms_a = self.__molA.sup(str_a)
        if str_b is not None: rms_b = self.__molB.sup(str_b)
        return rms_a, rms_b
    
    def get(self):
        """return a tuple of parameters for two fragments"""
        return self.__molA.get(), self.__molB.get()
        
    def get_mol(self):
        """returns molecular fragments"""
        return self.__molA, self.__molB
    
    # protected
    def _ex_rep_ct(self,nlo,shift,nmode,cunit):
        """calculate exchange-repulsion and charge-transfer properties"""
        # variables
        varA = self.__molA.get()
        varB = self.__molB.get()
        # basis sets
        bfsA = self.__molA.get_bfs()
        bfsB = self.__molB.get_bfs()
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
        faijc  = varA['fckc'].diagonal()
        cikca  = varA['vecc']
        qa     = varA['chlpg']
        redmss = varA['redmass']
        freq   = varA['freq']
        gijj   = varA['gijk'][nmode-1,nmode-1,:]
        eiglvc = varA['lvec']
        mlist  = bfsA.get_bfsl() + 1
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        fbijc= varB['fckc'].diagonal()
        cikcb= varB['vecc']
        qb   = varB['chlpg']
        # instantaneous integrals
        #import time
        #t0 = time.time()
        skm  = getSAB(bfsA,bfsB)
        tkm  = getTAB(bfsA,bfsB)
        #t1 = time.time()
        sk1m = getSA1B(bfsA,bfsB)
        tk1m = getTA1B(bfsA,bfsB)
        tkk =  getT(bfsA)
        tll =  getT(bfsB)
        #t2 = time.time()
        vkl =  getVEFP(bfsA,bfsB,qb,rnb)
        vlk =  getVEFP(bfsB,bfsA,qa,rna)
        vkm =  getVEFP(bfsA,bfsA,qb,rnb)
        vln =  getVEFP(bfsB,bfsB,qa,rna)
        #t3 = time.time()
        #vvv = getV(bfsA,self.__molA.atoms)
        #t4 = time.time()
        #print "STAB %.2f" % (t1-t0)
        #print "XA1B %.2f" % (t2-t1)
        #print "VXLN %.2f" % (t3-t2)
        #print "TOTL %.2f" % (t3-t0)
        #print "DUPA %.2f" % (t4-t3)
        # calculate the properties!
        shftma,shftea = shftce(redmss,freq,gijj,eiglvc,
                               rna,rnb,ria,rib,ria1,
                               cika,cikb,cikca,cikcb,cika1,
                               skm,tkm,tkk,tll,vkm,vkl,vlk,vln,sk1m,tk1m,
                               faij,fbij,faijc,fbijc,faij1,
                               za,zb,mlist,nmode)
        if cunit:
           shftma *= self.HartreePerHbarToCmRec
           shftea *= self.HartreePerHbarToCmRec
           
        return shftma, shftea

    def _ex_rep(self,nlo,shift,nmode,cunit):
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
    
class FragFactory(object,DIFF):
    """
Solvatochromic Effective Fragment Potential Fragment
----------------------------------------------------

Usage:
"""
    def __init__(self,anh=None,basis=None,nae=None,
                      fchk=None,gmslog=None,
                      chelpg=None,esp=None,):
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
        self.__chlpg  = chelpg
        self.__esp    = esp
        self._init()
        self._create()

    # public
    
    def set(self,anh=None,basis=None,nae=None,
                 fchk=None,gmslog=None,
                 chelpg=None,esp=None,
                 dpol=None,dpol1=None,):
        """set the properties to the object"""
        if self.__anh    is not None: self.__anh    = anh
        if self.__fchk   is not None: self.__fchk   = fchk
        if self.__gmslog is not None: self.__gmslog = gmslog
        if self.__basis  is not None: self.__basis  = basis
        if self.__nae    is not None: self.__nae    = nae
        if self.__chlpg  is not None: self.__chlpg  = chelpg
        if self.__esp    is not None: self.__esp    = esp
        if self.__dpol   is not None: self.__dpol   = dpol
        if self.__dpol1  is not None: self.__dpol1  = dpol1
        return

    def reset(self,anh=None,basis=None,nae=None,
                   fchk=None,gmslog=None,
                   chelpg=None,esp=None,
                   dpol=None,dpol1=None,):
        """reset the properties"""
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
        self.__chlpg  = chelpg
        self.__esp    = esp
        self.__dpol   = dpol
        self.__dpol1  = dpol1
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
        if self.__vecc  is not None: par['vecc' ] = self.__vecc
        if self.__vecc1 is not None: par['vecc1'] = self.__vecc1
        if self.__fckc  is not None: par['fckc' ] = self.__fckc
        if self.__fckc1 is not None: par['fckc1'] = self.__fckc1
        if self.__esp   is not None: par['esp'  ] = self.__esp
        if self.__chlpg is not None: par['chlpg'] = self.__chlpg
        if self.__dpol  is not None: par['dpol' ] = self.__dpol
        if self.__dpol1 is not None: par['dpol1'] = self.__dpol1
        return par
    
    def eval(self,ct=False):
        """Parses AO-LMO transformation matrix and Fock matrix.
Transforms the latter from AO to LMO space. Computes also 
overlap integrals and parses density matrix. If ct is True
the canonical Fock matrix and vectors will be saved."""
        assert self.__mol is not None, 'molecule not specified! (no fchk file)'
        # evaluate transformation matrices and LMO centroids
        SAO   = PyQuante.Ints.getS(self.__bfs)
        dmat = ParseDmatFromFchk(self.__fchk,self.__basis_size)
        vecc = ParseVecFromFchk(self.__fchk)
        veccocc= vecc[:self.__nae,:]
        tran, veclmo = get_pmloca(self.__natoms,mapi=self.__bfs.LIST1,sao=SAO,
                                  vecin=veccocc,nae=self.__nae,
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
        Fock = ParseFockFromGamessLog(self.__gmslog,interpol=False)
        fock = tensordot(veclmo,tensordot(veclmo,Fock,(1,0)),(1,1))
        if ct: fckc = tensordot(vecc  ,tensordot(vecc  ,Fock,(1,0)),(1,1))
        # save
        self.__lmoc = dma.get_origin()[self.__natoms:]
        self.__tran = tran
        self.__vecl = veclmo
        if ct: self.__vecc = vecc
        self.__fock = fock
        if ct: self.__fckc = fckc
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
        self.__vecl   = None; self.__vecl1  = None; self.__vecc = None
        self.__vecc1  = None; self.__fckc   = None; self.__fckc1= None
        self.__chlpg  = None; self.__esp    = None; self.__dpol = None
        self.__dpol1  = None
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
    