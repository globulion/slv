# -------------------------------------------------------------------------- #
#            SOLVATOCHROMIC EFFECTIVE FRAGMENT POTENTIAL MODULE              #
# -------------------------------------------------------------------------- #

from numpy     import tensordot, dot, transpose, array, float64, zeros, \
                      concatenate as con, where
from units     import *
from utilities import Read_xyz_file, get_pmloca, ParseFockFromGamessLog, \
                      ParseDmatFromFchk, ParseVecFromFchk, MakeMol
from diff      import DIFF
from PyQuante.Ints import getT, getSAB, getTAB, getSA1B, getTA1B, getVEFP, \
                          getV, getbasis
from shftex import shftex
from shftce import shftce
from solpol import solpol, mollst, sftpol
from efprot import tracls
from exrep  import exrep
import sys, copy, os, re, math, glob, PyQuante.Ints, coulomb.multip, clemtp
sys.stdout.flush()

__all__ = ['EFP','FragFactory',]
__version__ = '1.0.2'

class EFP(object,UNITS):
    """"""
    def __init__(self,ccut=None,pcut=None,ecut=None,#pairwise_all=False,
                      elect=True,pol=False,rep=False,ct=False,disp=False,all=False,
                      corr=False,
                      nlo=False,freq=False,cunit=False):
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
        self.__cunit = cunit
        #
        if self.__ccut is None: self.__ccut = 1.E+10
        else:                   self.__ccut = ccut
        if self.__pcut is None: self.__pcut = 1.E+10
        else:                   self.__pcut = pcut
        if self.__ecut is None: self.__ecut = 1.E+10
        else:                   self.__ecut = ecut
        #
        if all:
           self.__eval_elect = True
           self.__eval_pol   = True
           self.__eval_rep   = True
           self.__eval_ct    = True
           self.__eval_disp  = True
           self.__eval_corr  = True
        else:
           self.__eval_elect = elect
           self.__eval_pol   = pol
           self.__eval_rep   = rep
           self.__eval_ct    = ct
           self.__eval_disp  = disp
           self.__eval_corr  = corr
        #
        if freq: 
           self.__pairwise_all = False
           self.__eval_freq = True

        else: self.__pairwise_all = True
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

    def set(self,pos,ind,nmol,bsm=None, supl=None):
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
        if supl is not None: 
           self.__suplist = supl
           self.__suplist_c = supl[0]
        else: self.__suplist_c = None
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
    
    def get_shift(self):
        """return frequency shift data"""
        return self.__shift
    
    def get_rms(self):
        """return RMS of superimposition of central molecule with its parameters (relevant for central molecule mode)"""
        return self.__rms_central
    
    def eval(self,lwrite=False):
        """evaluate the properties"""
        if self.__pairwise_all:
        # ============================================================================== #
        #                                                                                #
        #                               PAIRWISE ALL MODE                                #
        #                       (interatction energy, NLO property                       #
        #                                                                                #
        # ============================================================================== #
           N = len(self.__nmol)
           PAR = []
           QO  = []
           nm_sum = 0
           for i in range(N):
               nm = self.__nmol[i]
               im = self.__ind[i]
               nm_sum += nm
               #
               STR = self.__rcoordc[nm_sum-nm:nm_sum]
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
              if self.__cunit: e_el  *= self.HartreeToKcalPerMole
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
                 if self.__cunit: e_pol  *= self.HartreeToKcalPerMole
                 if lwrite: print " Polarization       energy: %10.6f"%e_pol
              
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
               if self.__cunit: e_rep  *= self.HartreeToKcalPerMole
               if lwrite: print " Exchange-repulsion energy: %10.6f"%e_rep
                              
        else:
        # ============================================================================== #
        #                                                                                #
        #                             CENTRAL MOLECULE MODE                              #
        #                      (interaction energy,frequency shift)                      #
        #                                                                                #
        # ============================================================================== #
           shift_total = 0.0; corr = 0.0
           rf2 = 0.0; rf3 = 0.0; rf4 = 0.0
           rk2 = 0.0; rk3 = 0.0; rk4 = 0.0
           corr_b, corr_c, corr_d = 0.0, 0.0, 0.0
           #
           N = len(self.__ntc)
           nmols = N+1
           PAR = []
           QO  = []
           ### central molecule
           frg = self.__bsm[0].copy()
           self.__rms_central = frg.sup( self.__rc, self.__suplist_c )
           if lwrite: print "Central rms: ",self.__rms_central
           parc= frg.get()
           #
           self.__mode = parc['mode']
           #
           PAR.append( parc )
           #
           qadc, octc = tracls( parc['dmaq'], parc['dmao'] )
           QO.append( (qadc,octc) )
           #
           chgc1 = parc['dmac1'].ravel()
           dipc1 = parc['dmad1'].ravel()
           qadc1, octc1 = frg.get_traceless_1(ravel=True)
           #
           chgc2 = parc['dmac2'].ravel()
           dipc2 = parc['dmad2'].ravel()
           qadc2, octc2 = frg.get_traceless_2()
           qadc2 = qadc2.ravel()
           octc2 = octc2.ravel()
           #
           gijj   = parc['gijk'][:,self.__mode-1,self.__mode-1]
           freq   = parc['freq']
           redmss = parc['redmass']
           lvec   = parc['lvec'].ravel()
           nmodes = parc['nmodes']
           #
           freqc = freq[self.__mode-1]
           redmssc=redmss[self.__mode-1]
           #
           ### other molecules
           for i in range(N):
               nm = self.__ntc[i]
               im = self.__mtc[i]
               #
               STR = self.__rcoordc[nm*i:nm*(i+1)]
               frg = self.__bsm[im].copy()
               rms = frg.sup( STR , suplist= self.__suplist[self.__ind[im]] )
               if lwrite: print "rms C: ",rms
               par = frg.get()
               PAR.append( par )
               #
               qad, oct = tracls( par['dmaq'], par['dmao'] )
               QO.append( (qad,oct) )
               #
           # ----------------------------------- ELECT --------------------------------- #
           if self.__eval_elect:
              ndma = [ x['ndma'] for x in PAR ]
              ndmac= parc['ndma']
              ndmas= sum(ndma)
              
              rdma = con  ([ x['rdma'] for x in PAR ]).reshape(ndmas*3)
              chg  = con  ([ x['dmac'] for x in PAR ]).reshape(ndmas)
              dip  = con  ([ x['dmad'] for x in PAR ]).reshape(ndmas*3)
              qad  = con  ([ QO[x][0]  for x in range(N+1)   ]).reshape(ndmas*6)
              oct  = con  ([ QO[x][1]  for x in range(N+1)   ]).reshape(ndmas*10)
              
              # mechanical anharmonicity
              mea,a,b,c,d,e = clemtp.sdmtpm(rdma,ndma,chg,dip,qad,oct,
                                              chgc1,dipc1,qadc1,octc1,
                                              redmss,freq,gijj,
                                              ndmac,self.__mode,lwrite)
              # electronic anharmonicity
              ea ,a,b,c,d,e = clemtp.sdmtpe(rdma,ndma,chg,dip,qad,oct,
                                              chgc2,dipc2,qadc2,octc2,
                                              redmss,freq,
                                              self.__mode,lwrite)
              # correction terms
              if self.__eval_corr:
                 corr,rf2,rf3,rf4,rk2,rk3,rk4,corr_b,corr_c,corr_d = \
                              clemtp.dmtcor(rdma,ndma,chg,dip,qad,oct,
                                            chgc1,dipc1,qadc1,octc1,
                                            redmss,freq,gijj,lvec,
                                            ndmac,self.__mode,lwrite)
              # change units from A.U. to specific units
              if self.__cunit:
                    #eel  *= self.HartreeToKcalPerMole
                    mea *= self.HartreePerHbarToCmRec
                    ea  *= self.HartreePerHbarToCmRec
                    corr_b   *= self.HartreePerHbarToCmRec
                    corr_c   *= self.HartreePerHbarToCmRec
                    corr_d   *= self.HartreePerHbarToCmRec
                    rf2 *= self.HartreePerHbarToCmRec
                    rk2 *= self.HartreePerHbarToCmRec
                    rf3 *= self.HartreePerHbarToCmRec
                    rk3 *= self.HartreePerHbarToCmRec
                    rf4 *= self.HartreePerHbarToCmRec
                    rk4 *= self.HartreePerHbarToCmRec
                    corr*= self.HartreePerHbarToCmRec
                    
              shift_total    += mea+ea+corr_d
              self.__shift[0] = mea
              self.__shift[1] =  ea
              #self.__shift[7] = rf2+rf3+rf4
              #self.__shift[8] = rk2+rk3+rk4
              self.__shift[4] = rf2+rf3+rf4
              self.__shift[5] = rk2+rk3+rk4

              if lwrite: 
                 print " Electrostatic  MEA frequency shift: %10.2f"%mea
                 print " Electrostatic   EA frequency shift: %10.2f"% ea
                 print " Electrostatic  TOT frequency shift: %10.2f"% shift_total
                 print " Electrostatic CORR frequency shift: %10.2f"% corr
              del PAR, QO
           # ------------------------------------- POL --------------------------------- #
              if self.__eval_pol:
                 npolc = parc['npol']
                 ### central molecule
                 N = len(self.__ntp)
                 PAR = []
                 QO  = []
                 PAR.append( parc )
                 QO.append( (qadc,octc) )
                 ### other molecules
                 for i in range(N):
                     nm = self.__ntp[i]
                     im = self.__mtp[i]
                     #
                     STR = self.__rcoordp[nm*i:nm*(i+1)]
                     frg = self.__bsm[im].copy()
                     rms = frg.sup( STR, suplist= self.__suplist[self.__ind[im]])
                     #if lwrite: print "rms P: ",rms
                     par = frg.get()
                     PAR.append( par )
                     #
                     qad, oct = tracls( par['dmaq'], par['dmao'] )
                     QO.append( (qad,oct) )
                     #
                 ndma = [ x['ndma'] for x in PAR ]
                 ndmas= sum(ndma)
                 #
                 rdma = con  ([ x['rdma'] for x in PAR ]).reshape(ndmas*3)
                 chg  = con  ([ x['dmac'] for x in PAR ]).reshape(ndmas)
                 dip  = con  ([ x['dmad'] for x in PAR ]).reshape(ndmas*3)
                 qad  = con  ([ QO[x][0]  for x in range(N+1)   ]).reshape(ndmas*6)
                 oct  = con  ([ QO[x][1]  for x in range(N+1)   ]).reshape(ndmas*10)
                 #
                 npol = [ x['npol'] for x in PAR ]
                 npols= sum(npol)
                 #
                 rpol = con  ([ x['rpol'] for x in PAR ]).reshape(npols*3)
                 pol  = con  ([ x['dpol'] for x in PAR ]).reshape(npols*9)
                 #
                 rpol1= parc['lmoc1'].reshape(nmodes*npolc*3)
                 pol1 = parc['dpol1'].reshape(nmodes*npolc*9)
                 #
                 DIM  = npols*3
                 #rpol1.fill(0.0)
                 #pol1.fill(0.0)
                 #
                 dmat = zeros((DIM,DIM),float64)
                 dimat= zeros((DIM,DIM),float64)
                 mat1 = zeros((DIM,DIM),float64)
                 #
                 flds = zeros( DIM,float64)
                 dipind=zeros( DIM,float64)
                 sdipnd=zeros( DIM,float64)
                 vec1  =zeros( DIM,float64)
                 fivec =zeros( DIM,float64)
                 avec  =zeros( DIM,float64)
                 #
                 #pol1.fill(0.0)
                 epol, spol = sftpol(rdma,chg,dip,qad,oct,
                                      chgc1,dipc1,qadc1,octc1,
                                      rpol,pol,dmat,flds,dipind,dimat,fivec,
                                      sdipnd,avec,vec1,mat1,
                                      redmss,freq,gijj,rpol1,pol1,lvec,
                                      ndma,npol,self.__mode,ndmac,npolc,lwrite=False)
                 if self.__cunit:
                    epol  *= self.HartreeToKcalPerMole
                    spol  *= self.HartreePerHbarToCmRec
                 shift_total += spol
                 self.__shift[2] = spol
                 if lwrite: 
                    print " Polarization       frequency shift: %10.2f"%spol
                    #print " Polarization energy         : %10.6f"%epol
                 del PAR, QO
           # ---------------------------------- EX-REP --------------------------------- #
           if  self.__eval_rep:
               N = len(self.__nte)
               PAR = []
               ### other molecules
               for i in range(N):
                   nm = self.__nte[i]
                   im = self.__mte[i]
                   #
                   STR = self.__rcoorde[nm*i:nm*(i+1)]
                   frg = self.__bsm[im].copy()
                   rms = frg.sup( STR , suplist= self.__suplist[self.__ind[im]] )
                   if lwrite: print "rms E: ",rms
                   par = frg.get()
                   PAR.append( par )
                   #
               serp = 0
               #for par in PAR: shift += self._pair_rep_freq(parc,par)
               # basis sets
               molA = MakeMol(parc['atno'],parc['pos'])
               bfsA = getbasis(molA,parc['basis'])
               nbsa = parc['nbasis']
               nmosa= parc['nmos']
               nmodes = parc['nmodes']
               # parameters for central molecule
               faij = parc['fock' ]               .ravel()
               faij1= parc['fock1']               .ravel()
               cika = parc['vecl' ]               .ravel()
               cika1= parc['vecl1']               .ravel()
               za   = parc['atno' ]
               rna  = parc['pos'  ]               .ravel()
               ria  = parc['lmoc' ]               .ravel()
               ria1 = parc['lmoc1']               .ravel()
               mlist= bfsA.get_bfsl() + 1
               redmss= parc['redmass']
               gijj = parc['gijk'][:,self.__mode-1,self.__mode-1]
               freq = parc['freq']
               lvec = parc['lvec']                .ravel()
               #
               for par in PAR:
                   molB = MakeMol(par['atno'],par['pos'])
                   bfsB = getbasis(molB,par['basis'])
                   nbsb = par['nbasis']
                   nmosb = par['nmos']
                   # instantaneous integrals
                   skm  = getSAB(bfsA,bfsB)       .ravel()
                   tkm  = getTAB(bfsA,bfsB)       .ravel()
                   sk1m = getSA1B(bfsA,bfsB)      .ravel()
                   tk1m = getTA1B(bfsA,bfsB)      .ravel()
                   ### molecule B
                   fbij = par['fock']             .ravel()
                   cikb = par['vecl']             .ravel()
                   zb   = par['atno']
                   rnb  = par['pos' ]             .ravel()
                   rib  = par['lmoc']             .ravel()
                   # calculate the properties!
                   #sma ,shftea = shftex(redmss,freq,gijj,lvec,
                   #                     ria,rib,rna,rnb,ria1,
                   #                     cika,cikb,cika1,
                   #                     skm,tkm,sk1m,tk1m,
                   #                     za,zb,mlist,
                   #                     faij,fbij,faij1,self.__mode)
                   sij = zeros(nmosa*nmosb,float64)
                   tij = zeros(nmosa*nmosb,float64)
                   smij= zeros(nmodes*nmosa*nmosb,float64)
                   tmij= zeros(nmodes*nmosa*nmosb,float64)
                   fi  = zeros(nmodes,float64)
                   #
                   sma, shftea = shftex(redmss,freq,gijj,lvec,
                                        ria,rib,rna,rnb,ria1,
                                        cika,cikb,cika1,
                                        skm,tkm,sk1m,tk1m,
                                        za,zb,nbsb,mlist,
                                        faij,fbij,faij1,self.__mode,
                                        sij,tij,smij,tmij,fi)
                   serp+=sma
                   #
               if self.__cunit:
                  serp *= self.HartreePerHbarToCmRec
               shift_total += serp
               self.__shift[3] = serp
               self.__shift[6] = shift_total
               if lwrite: 
                  print " Exchange-repulsion frequency shift: %10.2f"%serp
                  print " ----------------------------------------------- "
                  print " TOTAL FREQUENCY SHIFT             : %10.2f"%shift_total
                  print " TOTAL FREQUENCY SHIFT (MEA)       : %10.2f"%(shift_total-ea)
        #
        #
        return # ========================================================================== # 
    
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
        self.__shift = zeros(9,dtype=float64)
        self.__rms_central = None
        self.__suplist = None
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
           ie = zeros(natoms-nac,dtype=bool)
           icm = zeros(nm,dtype=bool)
           ipm = zeros(nm,dtype=bool)
           iem = zeros(nm,dtype=bool)
           #
           nccut, npcut, necut, mccut, mpcut, mecut = \
                                        mollst(rc,rm,ic,ip,ie,icm,ipm,iem,
                                               self.__nmol[1:],
                                               self.__ccut,
                                               self.__pcut,
                                               self.__ecut)
                                               #
           # molecule coordinate lists
           rcoordc = pos[nac:][array(nccut,dtype=bool)]
           rcoordp = pos[nac:][array(npcut,dtype=bool)]
           rcoorde = pos[nac:][array(necut,dtype=bool)]
           # molecule type lists
           mtc= self.__ind[1:][array(mccut,dtype=bool)]
           mtp= self.__ind[1:][array(mpcut,dtype=bool)]
           mte= self.__ind[1:][array(mecut,dtype=bool)]
           # moleculear atom number lists
           ntc= self.__nmol[1:][array(mccut,dtype=bool)]
           ntp= self.__nmol[1:][array(mpcut,dtype=bool)]
           nte= self.__nmol[1:][array(mecut,dtype=bool)]
           #
           self.__rc = rc
           self.__rcoordc = rcoordc
           self.__rcoordp = rcoordp
           self.__rcoorde = rcoorde
           self.__mtc = mtc
           self.__mtp = mtp
           self.__mte = mte
           self.__ntc = ntc
           self.__ntp = ntp
           self.__nte = nte
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

    def _pair_rep_freq(self,varA,varB):
        """exchange-repulsion pair frequency shift"""
        # basis sets
        molA = MakeMol(varA['atno'],varA['pos'])
        molB = MakeMol(varB['atno'],varB['pos'])
        bfsA = getbasis(molA,varA['basis'])
        bfsB = getbasis(molB,varB['basis'])
        # instantaneous integrals
        skm  = getSAB(bfsA,bfsB)
        tkm  = getTAB(bfsA,bfsB)
        sk1m = getSA1B(bfsA,bfsB)
        tk1m = getTA1B(bfsA,bfsB)
        # parameters
        ### molecule A
        faij = varA['fock']
        faij1= varA['fock1']
        cika = varA['vecl']
        cika1= varA['vecl1']
        za   = varA['atno']
        rna  = varA['pos']
        ria  = varA['lmoc']
        ria1 = varA['lmoc1']
        mlist= bfsA.get_bfsl() + 1
        redmss= varA['redmass']
        gijj = varA['gijk'][:,self.__mode-1,self.__mode-1]
        freq = varA['freq']
        lvec = varA['lvec']
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # calculate the properties!
        shift ,shftea = shftex(redmss,freq,gijj,lvec,
                               ria,rib,rna,rnb,ria1,
                               cika,cikb,cika1,
                               skm,tkm,sk1m,tk1m,
                               za,zb,mlist,
                               faij,fbij,faij1,self.__mode)
        return shift
        


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
    
