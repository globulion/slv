# --------------------------------------------------------------------- #
#           SOLVATOCHROMIC FREQUENCY SHIFT CORRECTION MODULE            #
# --------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
sys.stdout.flush()

__all__ = ['SLVCOR',]
__version__ = '1.2.2'

class SLVCOR(UNITS):
    """\
It contains the implementation of additional terms arising from 
taking into account potential derivatives and applying local 
approximation to normal coordinate gradient operator:
grad(normal_coord) \approx L \dot \nabla. Here it is assumed,
that solvent DMA object, L-eigenvectors as well as all of the 
derivative tensors are properly rotated and translated!
Reduced masses should be given in AMU, harmonic frequencies
in cm-1 and other quantities in AU."""
    
    def __init__(self,fderiv=0,redmass=0,freq=0,L=[],
                 solute=0,solvent=0,mode_id=0,gijj=[],
                 ua_list=None):
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
                    R    = Ra[i]-Rb[j]
                    Rab=sqrt(sum(R**2,axis=0))
                    
                    ### fi terms
                    for I in xrange(self.__nModes):
                        rf2 = 0; rf3 = 0; rf4 = 0
                        #
                        rf2 -=   qa[i] * qb[j] * dot(self.__Lvec(I,i),R)               / Rab**3
                        #
                        rf3 -= 3*qa[i] * dot(self.__Lvec(I,i),R) * dot(Db[j],R)        / Rab**5
                        rf3 +=   qa[i] * dot(self.__Lvec(I,i),Db[j])                   / Rab**3
                        rf3 += 3*qb[j] * dot(self.__Lvec(I,i),R) * dot(Da[i],R)        / Rab**5
                        rf3 -=   qb[j] * dot(self.__Lvec(I,i),Da[i])                   / Rab**3
                        #
                        rf4 -= 5*qa[i] * dot(self.__Lvec(I,i),R) * dot(dot(Qb[j],R),R) / Rab**7
                        rf4 += 2*qa[i] * dot(dot(Qb[j],self.__Lvec(I,i)),R)            / Rab**5
                        #
                        f1   =   dot(Da[i],R) * dot(self.__Lvec(I,i),Db[j])
                        f1  +=   dot(self.__Lvec(I,i),R) * dot(Da[i],Db[j])
                        f1  +=   dot(Db[j],R) * dot(self.__Lvec(I,i),Da[i])
                        f1  *=-3                                                       / Rab**5
                        rf4 +=   f1
                        rf4 +=15*dot(Da[i],R) * dot(Db[j],R) * dot(self.__Lvec(I,i),R) / Rab**7
                        #
                        rf4 += 2*qb[j]*dot(dot(Qa[i],self.__Lvec(I,i)),R)              / Rab**5
                        rf4 -= 5*qb[j]*dot(dot(Qa[i],R),R) * dot(self.__Lvec(I,i),R)   / Rab**7
                        
                        Rn_fi[I,1]+= rf2
                        Rn_fi[I,2]+= rf3
                        Rn_fi[I,3]+= rf4
                        
                    ### kjj terms
                    rk2 = 0; rk3 = 0; rk4 = 0
                    fqa,fDa,fQa,fOa = self.__fdx
                    
                    rk2 -= 2*fqa[i] * qb[j] * dot(self.__Lvec(J,i),R)                  / Rab**3
                    #
                    rk3 -= 6*fqa[i] * dot(self.__Lvec(J,i),R) * dot(Db[j],R)           / Rab**5
                    rk3 += 2*fqa[i] * dot(self.__Lvec(J,i),Db[j])                      / Rab**3
                    rk3 += 6* qa[i] * qb[j] * dot(self.__Lvec(J,i),R)                  / Rab**5 # error
                    rk3 -=    qa[i] * qb[j] * dot(self.__Lvec(J,i),self.__Lvec(J,i))   / Rab**3 # error
                    rk3 += 6* qb[j] * dot(self.__Lvec(J,i),R) * dot(fDa[i],R)          / Rab**5
                    rk3 -= 2* qb[j] * dot(fDa[i],self.__Lvec(J,i))                     / Rab**3
                    #
                    rk4 += 4*fqa[i] * dot(dot(Qb[j],self.__Lvec(J,i)),R)               / Rab**5
                    rk4 -=10*fqa[i] * dot(self.__Lvec(J,i),R) * dot(dot(Qb[j],R),R)    / Rab**7
                    #
                    k1   = 2*dot(self.__Lvec(J,i),R) * dot(self.__Lvec(J,i),Db[j])
                    k1  +=   dot(self.__Lvec(J,i),self.__Lvec(J,i)) * dot(Db[j],R)
                    k1  *=-3*qa[i]                                                     / Rab**5
                    rk4 +=   k1
                    rk4 +=15*qa[i] * (dot(self.__Lvec(J,i),R))**2 * dot(Db[j],R)       / Rab**7
                    #
                    k1   = 2*dot(self.__Lvec(J,i),R) * dot(self.__Lvec(J,i),Da[i])
                    k1  +=   dot(self.__Lvec(J,i),self.__Lvec(J,i)) * dot(Da[i],R)
                    k1  *= 3*qb[j]                                                     / Rab**5
                    rk4 +=   k1
                    rk4 -=15*qb[j] * (dot(self.__Lvec(J,i),R))**2 * dot(Da[i],R)       / Rab**7
                    #
                    k1   =   dot(fDa[i],R) * dot(self.__Lvec(J,i),Db[j])
                    k1  +=   dot(self.__Lvec(J,i),R) * dot(fDa[i],Db[j])
                    k1  +=   dot(fDa[i],self.__Lvec(J,i)) * dot(Db[j],R)
                    k1  *=-6                                                           / Rab**5
                    rk4 +=   k1
                    #
                    rk4 +=30*dot(fDa[i],R) * dot(self.__Lvec(J,i),R) * dot(Db[j],R)    / Rab**7
                    #
                    rk4 += 4*qb[j] * dot(dot(fQa[i],self.__Lvec(J,i)),R)               / Rab**5
                    rk4 -=10*qb[j] * dot(dot(fQa[i],R),R) * dot(self.__Lvec(J,i),R)    / Rab**7
                    
                    Rn_kjj[1] += rk2
                    Rn_kjj[2] += rk3
                    Rn_kjj[3] += rk4
                    
            return Rn_fi, Rn_kjj
        
    def eval(self):
        """evaluate the corrections to shifts"""
        self.__fi,self.__kjj  = self.fi_kjj()
        self.__ma = self.ma()
        self.__ea = self.ea()
        corr = array(self.__ma) + array(self.__ea)
        # corrections added together, they can be added to self.shift[0] from SLV class
        self.corr = array([0,
                           corr[1],
                           corr[1]+corr[2],
                           corr[1]+corr[2]+corr[3],
                           0]) * self.HartreePerHbarToCmRec
        return
            
    def __repr__(self):
        """print the correction terms to the frequency shift"""
        A,B,C,D,E = self.shift
        log  = " CORRECTION TERMS [cm-1]         : CORRECTED SHIFTS [cm-1]  \n"
        log += " --------------------------------:--------------------------\n"
        a,b = self.get()
        a*=self.HartreePerHbarToCmRec;b*=self.HartreePerHbarToCmRec
        log += " %12s %10s         :  1        %10.2f\n"       % ('MA','EA',      A)
        log += " %3s %10.2f %10.2f       :  1+2      %10.2f\n" % ('R-2',a[1],b[1],B+a[1]+b[1])
        log += " %3s %10.2f %10.2f       :  1+2+3    %10.2f\n" % ('R-3',a[2],b[2],C+a[2]+b[2]+a[1]+b[1])
        log += " %3s %10.2f %10.2f       :  1+2+3+4  %10.2f\n" % ('R-4',a[3],b[3],D+a[3]+b[3]+a[2]+b[2]+a[1]+b[1])
        log += "                                 :  1+2+3+4+5%10s\n" % ('???')
        log += " --------------------------------:--------------------------\n"
        
        return str(log)
        
    def get(self):
        """return the corrections to shifts"""
        return array(self.__ma), array(self.__ea)
        
    def ma(self):
        """evaluates the contributions from mechanical anharmonicity
        to frequency shift"""
        j = self.__mode_id
        sum = zeros(4,dtype=float64)
        for i in range(self.__nModes):
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
        dma1.makeDMAfromFULL()
        dma1.MAKE_FULL()
        #
        Ra,qa,Da,Qa,Oa = dma1.DMA_FULL
        Rb,qb,Db,Qb,Ob = dma2.DMA_FULL
        self.__sx = (Ra,qa,Da,Qa,Oa)
        self.__sy = (Rb,qb,Db,Qb,Ob)
        del dma1
        ### derivatives of DMA objects
        dma1=self.__fderiv.copy()
        ##### ten zahaszowany kod dotyczy przypadku używania wszystkich pochodnych, jednak to nie jest konieczne i wytarczą j-th mode's pochodne
        #for i in range(self.__nModes):                                               
            #dma1[i].set_structure(pos=self.__solute.pos, origin=self.__solute.pos)   
            #print dma1[i].pos                                                        
            #dma1[i].MakeUa(ua_list,change_origin=True)                               
        #    dma1[i].MAKE_FULL()                                                      
        #    dma1[i].MakeTraceless()                                                  
        #    Ra,qa,Da,Qa,Oa = dma1[i].DMA_FULL                                        
        dma1.set_structure(pos=self.__solute.pos, origin=self.__solute.pos)
        if self.__ua_list is not None: dma1.MakeUa(self.__ua_list,change_origin=True)
        dma1.MAKE_FULL()
        dma1.MakeTraceless()
        Ra,qa,Da,Qa,Oa = dma1.DMA_FULL
        self.__fdx = (qa,Da,Qa,Oa)
        #self.__fdx.append((qa,Da,Qa,Oa))
        return
        
    def __weight(self):
        """weights L-matrix and derivative tensor elements
using reduced masses"""
        # L-matrix
        temp = sqrt(self.__redmass*self.AmuToElectronMass)[newaxis,:]
        self.__L = temp * self.__L
        # fderiv
        self.__fderiv *= sqrt(self.__redmass[self.__mode_id]*self.AmuToElectronMass)
        # cubic anharmonic constants
        temp = sqrt(self.__redmass)[:,newaxis,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,:,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,newaxis,:,]
        self.__gijj = temp * self.__gijj
