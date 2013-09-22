# --------------------------------------------------------------------- #
#                           HESSIAN MODULE                              #
# --------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
sys.stdout.flush()


class HESSIAN(UNITS):
    """\
It contains the implementation of Hessian. Here it is assumed,
that solvent DMA object, L-eigenvectors as well as all of the 
derivative tensors are properly rotated and translated!
Reduced masses should be given in AMU, harmonic frequencies
in cm-1 and other quantities in AU."""
    
    def __init__(self,fderiv=0,redmass=0,freq=0,L=[],
                 solute=0,solvent=0,mode_id=0,gijj=[],
                 ua_list=None,mol=0):
        # first and second DMA derivatives wrt normal mode
        self.__fderiv = copy.deepcopy(array(fderiv))
        # reduced masses in helico
        self.__redmass = array(redmass)
        self.__nModes = len(redmass)
        self.__nAtoms = (self.__nModes + 6)/3
        # harmonic frequencies in helico
        self.__freq = array(freq)
        self.__freq_mcho = array(freq) * self.CmRecToHz * self.HzToAuAngFreq
        # cubic anharmonic constants
        self.__gijj = array(gijj)
        # L matrix
        self.__L = array(L)
        # mode of interest for frequency shift
        self.__mode_id = mode_id
        # solute DMA (electrostatic nnormal moments!)
        self.__solute = solute.copy()
        # solvent DMA
        self.__solvent = solvent.copy()
        # ua_list
        self.__ua_list = ua_list
        # 
        self.__hess = None
        self.__mol = mol
        ### prepare solute and solvent DMA objects
        self.__prepare_dma()
                
        ### weight L-matrix and derivative tensors elements
        self.__weight()
        
        ### inline function withdrawing eigenvector for given mode i and atom x
        self.__Lvec  = lambda i,x: self.__L[:,i].copy().reshape(self.__nAtoms,3)[x]

    def __prepare_dma(self):
        fder = []
        for i in range(self.__nModes):
            dma1 = self.__fderiv[i].copy()
            if self.__ua_list is not None:
                dma1.MakeUa(self.__ua_list,change_origin=True)
            dma1.set_structure(pos=self.__solute.pos, origin=self.__solute.pos)
            fder.append(dma1)
        self.__fderiv = fder
        return
    
    def diag(self):
        """diagonalize Hessian. Returns new frequencies in [cm-1]"""
        if self.__hess is None: return 0
        else:
            print self.__freq
            n = 0
            L = self.__L#[:,:-n]
            H = self.__hess#[:-n,:-n]
            self.__hess_mwc = dot(dot(L,H),transpose(L))
            VecAnal = VIB(self.__mol,self.__hess_mwc)
            VecAnal.eval()
            self.__freq_new, self.__u = VecAnal.get()
            #self.__eigh,self.__u = linalg.eigh(self.__hess_mwc)
            #self.__freq_new = where(self.__eigh>0.,
            #                  sqrt(self.__eigh)\
            #                 *self.HartreePerHbarToCmRec,
            #                 -sqrt(-self.__eigh)\
            #                 *self.HartreePerHbarToCmRec)
        #print PUPA(self.__u)
        #a = dot(transpose(self.__u),dot(self.__hess_mwc,self.__u))
        #print PUPA(sqrt(a)*self.HartreePerHbarToCmRec)
        return self.__freq_new
    
    def approx(self):
        """evaluate approximate frequencies without Hessian diagonalization"""
        freq = zeros(self.__nModes,dtype=float64)
        for i in xrange(self.__nModes):
            f = sqrt(self.__hess[i,i]/self.__redmass[i])*self.HartreePerHbarToCmRec
            freq[i] =  f
        self.__freq_approx = freq
        return freq
    def eval(self):
        """evaluate the Hessian matrix"""
        self.__hess = zeros((self.__nModes,self.__nModes),dtype=float64)
        self.__fi = self.fi()
        for j in xrange(self.__nModes):
            for k in xrange(j+1):
                SUM = 0.0
                for i in xrange(self.__nModes):
                    SUM += self.__gijj[i,j,k] * self.__fi[i]/\
                    (self.__redmass[i]*self.__freq_mcho[i]**2)
                self.__hess[j,k]-= SUM
                if j==k:
                   self.__hess[j,k]+= self.__redmass[j]*self.__freq_mcho[j]**2
                   #self.__hess[j,k] /= sqrt(self.__redmass[j]*self.__redmass[k])
                   #print sqrt(self.__freq_mcho[j]**2 ) *self.HartreePerHbarToCmRec
                #else:
                   #self.__hess[j,k] /= sqrt(self.__redmass[j]*self.__redmass[k]) 
                   #self.__hess[k,j] = self.__hess[j,k]
        
        self.__freq_new  = self.diag()           
        self.__dq = self.dq()
        return
    
    def fi(self):
        """calculate fi tensors. result in AU"""
        f_i = []
        for i in xrange(self.__nModes):
            f = FrequencyShift(solute=self.__fderiv[i],
                               solvent=self.__solvent,
                               solute_structure=self.__solute.get_origin())[3]
            f_i.append(f)
        f_i = array(f_i,dtype=float64)
        f_i/= self.HartreePerHbarToCmRec
        return f_i
    
    def dq(self):
        """evaluate structural distortion"""
        dQ = []
        f_i = self.fi()
        for i in range(self.__nModes): 
            dq = f_i[i] / (self.__redmass[i]*self.__freq_mcho[i]**2)
            dQ.append(dq)
        dQ = array(dQ,float64)
        return dQ
    
    def __repr__(self):
        """print the correction terms to the frequency shift"""
        log = '\n'
        log+= ' NORMAL MODE DISPLACEMENTS DUE TO SOLVATION\n'
        log+= ' %s %s\n'%('Mode'.rjust(4),'dQ [A.U.]'.rjust(10))
        for i in range(self.__nModes):
            log+= ' %4i %10.6f\n' % ((i+1),self.__dq[i])
        log+= ' NEW HARMONIC FREQUENCIES\n'
        log+= ' %s %s\n'%('Mode'.rjust(4),'Freq [cm-1]'.rjust(10))
        for i in xrange(len(self.__freq_new)):
            log+= ' %4i %10.2f\n' % ((i+1),self.__freq_new[i])
        return str(log)
        
    def get(self):
        """return the Hessian in MWC"""
        return self.__hess_mwc
        
    def __weight(self):
        """weights L-matrix and derivative tensor elements
using reduced masses"""
        # L-matrix
        temp = sqrt(self.__redmass*self.AmuToElectronMass)[newaxis,:]
        #self.__L = self.__L / temp
        # fderiv
        for i in xrange(self.__nModes):
            self.__fderiv[i] *= sqrt(self.__redmass[i]*self.AmuToElectronMass)
        # cubic anharmonic constants
        temp = sqrt(self.__redmass)[:,newaxis,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,:,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,newaxis,:,]
        self.__gijj = temp * self.__gijj
        # reduced masses
        self.__redmass *= self.AmuToElectronMass