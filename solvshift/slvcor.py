# --------------------------------------------------------------------- #
#           SOLVATOCHROMIC FREQUENCY SHIFT CORRECTION MODULE            #
# --------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
sys.stdout.flush()

__all__ = ["SLVCOR"]

class SLVCOR(UNITS):
    """It contains the implementation of additional terms arising from 
    taking into account potential derivatives and applying local 
    approximation to normal coordinate gradient operator:
    grad(normal_coord) \approx L \dot \nabla. Here it is assumed,
    that solvent DMA object, L-eigenvectors as well as all of the 
    derivative tensors are properly rotated and translated!
    Reduced masses should be given in AMU, harmonic frequencies
    in cm-1 and other quantities in AU."""
    
    def __init__(self,fderiv=0,sderiv=0,redmass=0,freq=0,L=[],
                 solute=0,solvent=0,mode_id=0,gijj=[]):
        # first and second DMA derivatives wrt normal mode
        self.__fderiv = array(fderiv)
        self.__sderiv = array(sderiv)
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
        
        ### weight L-matrix and derivative tensors elements
        self.__weight()
        
        ### prepare solute and solvent DMA objects
        self.__prepare_dma()
        
        ### inline function withdrawing eigenvector for given mode i and atom x
        self.__Lvec  = lambda i,x: self.__L[:,i].copy().reshape(self.__nAtoms,3)[x]
                    
    def fi_kjj(self):
            """evaluate first- and second-derivative terms"""
            # R-1, R-2, R-3 and R-4 array of terms
            Rn_fi = zeros((self.__nModes,4),dtype=float64)
            Rn_kjj= zeros(4,dtype=float64)
            
            Ra,qa,Da,Qa,Oa = self.__sx
            Rb,qb,Db,Qb,Ob = self.__sy
            J = self.__mode_id
            
            for i in xrange(len(Ra)):
                for j in xrange(len(Rb)):
                    R    = Rb[j]-Ra[i]
                    Rab=sqrt(sum(R**2,axis=0))
                    
                    rf2 = 0; rf3 = 0; rf4 = 0
                    rk2 = 0; rk3 = 0; rk4 = 0
                    ### fi terms
                    for I in xrange(self.__nModes):
                        rf2 -=   qa[i] * qb[j] * dot(self.__Lvec(I,i),R)        / Rab**3
                        rf3 -= 3*qa[i] * dot(self.__Lvec(I,i),R) * dot(Db[j],R) / Rab**5
                        rf3 +=   qa[i] * dot(self.__Lvec(I,i),Db[j])            / Rab**3
                        
                        Rn_fi[I,1]+= rf2
                        Rn_fi[I,2]+= rf3
                        
                    ### kjj terms
                    fqa,fDa,fQa,fOa = self.__fdx[J]
                    
                    rk2 -= 2*fqa[i] * qb[j] * dot(self.__Lvec(J,i),R)        / Rab**3
                    rk3 -= 6*fqa[i] * dot(self.__Lvec(J,i),R) * dot(Db[j],R) / Rab**5
                    rk3 += 2*fqa[i] * dot(self.__Lvec(J,i),Db[j])            / Rab**3
                    rk3 += 6* qa[i] * qb[j] * dot(self.__Lvec(J,i),R)        / Rab**5
                    rk3 -=    qa[i] * qb[j] * dot(self.__Lvec(J,i),self.__Lvec(J,i)) / Rab**3
                    
                    Rn_kjj[1] += rk2
                    Rn_kjj[2] += rk3
                    
            return Rn_fi, Rn_kjj
        
    def eval(self):
            """evaluate the corrections to shifts"""
            self.__fi,self.__kjj  = self.fi_kjj()
            self.__ma = self.ma()
            self.__ea = self.ea()
            return
            
    def __repr__(self):
            """print the correction terms to the frequency shift"""
            log  = '\n'
            a,b = self.get()
            log += 4*'%13.2f'% tuple( a * self.HartreePerHbarToCmRec )
            log += '\n'
            log += 4*'%13.2f'% tuple( b * self.HartreePerHbarToCmRec )
            log += '\n'
            return str(log)
        
    def get(self):
            """return the corrections to shifts"""
            return self.__ma, self.__ea
        
    def ma(self):
            """evaluates the contributions from mechanical anharmonicity
            to frequency shift"""
            j = self.__mode_id
            sum = zeros(4,dtype=float64)
            for i in range(self.__nModes):#sqrt(self.__redmass[i]*self.AmuToElectronMass)*
                sum -= self.__fi[i]*(self.__gijj[i,j,j]/
                (self.__redmass[i]*self.AmuToElectronMass*self.__freq_mcho[i]**2))
            sum /= 2.0 * self.__redmass[j]*self.AmuToElectronMass * self.__freq_mcho[j]
            
            return sum
        
    def ea(self):
            """evaluates the contributions from electronic anharmonicity
            to frequency shift"""
            j = self.__mode_id
            sum = self.__kjj.copy() #* self.__redmass[j]*self.AmuToElectronMass
            sum/= 2.0 * self.__redmass[j]*self.AmuToElectronMass * self.__freq_mcho[j]
            
            return sum
        
    def __prepare_dma(self):
            """preparation of DMA: making them in traceless forms and 
            projection to explicit array format from DMA format"""
            ### normal electrostatic DMA objects
            dma1=self.__solute.copy()
            dma2=self.__solvent.copy()
            # make FULL format of DMA distribution
            dma1.MAKE_FULL()
            dma2.MAKE_FULL()
            # transform FULL format to fraceless forms for quadrupoles and octupoles
            dma1.MakeTraceless()
            #
            Ra,qa,Da,Qa,Oa = dma1.DMA_FULL
            Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
            self.__sx = (Ra,qa,Da,Qa,Oa)
            self.__sy = (Rb,qb,Db,Qb,Ob)
            
            ### derivatives of DMA objects
            dma1=self.__fderiv.copy()
            #dma2=self.__sderiv.copy()
            self.__fdx = []
            for i in range(self.__nModes):
                dma1[i].MAKE_FULL()
                dma1[i].MakeTraceless()
                Ra,qa,Da,Qa,Oa = dma1[i].DMA_FULL
                self.__fdx.append((qa,Da,Qa,Oa))
            #dma2.MAKE_FULL()
            #dma2.MakeTraceless()
            #Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
            #self.__sdx = (qb,Db,Qb,Ob)
            return
        
    def __weight(self):
            """weights L-matrix and derivative tensor elements
            using reduced masses"""
            # L-matrix
            temp = sqrt(self.__redmass*self.AmuToElectronMass)[newaxis,:]
            self.__L = temp * self.__L
            # fderiv & sderiv
            self.__fderiv *= sqrt(self.__redmass*self.AmuToElectronMass)
            #self.__sderiv *=      self.__redmass*self.AmuToElectronMass
            # cubic anharmonic constants
            temp = sqrt(self.__redmass)[:,newaxis,newaxis,]
            self.__gijj = temp * self.__gijj
            temp = sqrt(self.__redmass)[newaxis,:,newaxis,]
            self.__gijj = temp * self.__gijj
            temp = sqrt(self.__redmass)[newaxis,newaxis,:,]
            self.__gijj = temp * self.__gijj