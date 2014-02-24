# --------------------------------------------------------------------- #
#                           HESSIAN MODULE                              #
# --------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy, pylab, solscf
sys.stdout.flush()

__all__ = ['HESSIAN',]
__version__ = '1.0.2'

def get_dq_scf(B1,G,f,dQ,maxit=10000000,threshold=1E-06):
    """
Calculate dQ vector using iterative SOLSCF theory

Description:

Evaluates dQ from quadratic tensor equation 

F + B·dQ + G:dQ² = 0

using iterative scheme.

Usage:

dQ, niters = get_dq_scf(B1,G,f,maxit,threshold)

Notes:

B1    = inverse ov B matrix; B = M + K
G     = gijk/2 tensor
f     = fi vectors
maxit = maximum number ov iterations
threshold = thershold of iteration convergence
"""
    dq, niters = solscf.solscf(B1,G,f,maxit,threshold,dQ)
    return dq, niters

class HESSIAN(UNITS):
    """
It contains the implementation of Hessian. Here it is assumed,
that solvent DMA object, L-eigenvectors as well as all of the 
derivative tensors are properly rotated and translated!
Reduced masses should be given in AMU, harmonic frequencies
in cm-1 and other quantities in AU.

Usage:

a = HESSIAN(fderiv,redmass,freq,L,solute,solvent,mode_id,gijj,
            ua_list,mol)
            
a.eval()        - evaluate properties: structural distortions,
                  force constant matrix and frequencies
A,B,C = a.get() - obtain structural distortions, 
                  Hessian matrix and frequencies
A = a.get_hess()- obtain Hessian matrix only
A = a.get_freq()- obtain frequencies only
f = a.approx()  - obtain approximated frequencies

Notes:

1) Hessian matrix is returned in A.U.
2) Structural distortions are returned in A.U.
3) Frequencies are returned in [cm-1]
"""
    
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
        # molecule
        self.__mol  = mol
        # status of the computation
        self.__eval = False
        self.__theory = 0
        self.__hess = None
        self.__iter = None
        self.__freq_new = None
        
        ### prepare solute and solvent DMA objects
        self._prepare_dma()
                
        ### weight L-matrix and derivative tensors elements
        self._weight()
        
        ### inline function withdrawing eigenvector for given mode i and atom x
        self._Lvec  = lambda i,x: self.__L[:,i].copy().reshape(self.__nAtoms,3)[x]
    
    # public
    
    def get(self):
        """return all properties: dQ, Hessian and frequencies"""
        return self.__dq_th, self.__hess, self.__freq_new

    def get_hess(self):
        """return the Hessian in normal coordinate space Qsqrt(M)"""
        return self.__hess

    def get_freq(self):
        """return new frequencies"""
        return self.__freq_new

    def eval(self,theory=0,zero=False,threshold=1.0E-6,max_iter=10000):
        """evaluate the Hessian matrix, frequencies and structural distortions"""
        self.__eval = True
        self.__theory = theory
        self.__hess = zeros((self.__nModes,self.__nModes),dtype=float64)
        self.__fi = self.fi()
        self.__dq = self.dq(theory=0)
        self.__dq_th = self.dq(theory=theory,zero=zero,threshold=threshold,max_iter=max_iter)
        if theory: dq = self.__dq_th
        else:      dq = self.__dq
        #self.__fi = where(abs(self.__dq)<0.1,self.__fi,0.0)
        #self.__fi[-1] = 0.0
        SUM = tensordot(self.__gijj,dq,(0,0))
        for j in xrange(self.__nModes):
            for k in xrange(j+1):
                self.__hess[j,k]-= SUM[j,k]/sqrt(self.__redmass[j]*self.__redmass[k])
                if j==k:
                   self.__hess[j,k]+= self.__freq_mcho[j]**2
                else:
                   self.__hess[k,j] = self.__hess[j,k]
        
        self.__freq_new  = self.diag()
        self.__freq_approx = self.approx()
        return

    def diag(self):
        """diagonalize Hessian. Returns new frequencies in [cm-1]"""
        if self.__hess is None: return 0
        else:
            print self.__freq
            n = 0
            L = dot(self._M1(),self.__L)#[:,:-n]
            H = self.__hess#[:-n,:-n]
            self.__hess_mwc = dot(dot(L,H),transpose(L))
            VecAnal = VIB(self.__mol,self.__hess_mwc,weight=False)
            VecAnal.eval()
            self.__freq_new, self.__redmass_new, self.__u = VecAnal.get()
        return self.__freq_new
    
    def approx(self):
        """
Evaluate approximate frequencies without Hessian diagonalization, 
using appropriate level of SOL-X theory"""
        diag_hess = diag(self.__hess)
        diag_hess = where(diag_hess>0,diag_hess,0)
        #freq = where(diag_hess>0,sqrt(diag_hess),sqrt(-diag_hess) )
        freq = sqrt(diag_hess)
        freq*= self.HartreePerHbarToCmRec
        self.__freq_approx = freq
        print freq
        return freq

    def fi(self):
        """calculate fi vector. result in AU"""
        f_i = []
        for i in xrange(self.__nModes):
            f = FrequencyShift(solute=self.__fderiv[i],
                               solvent=self.__solvent,
                               solute_structure=self.__solute.get_origin())[3]
            f_i.append(f)
        f_i = array(f_i,dtype=float64)
        f_i/= self.HartreePerHbarToCmRec
        return f_i
    
    def dq(self,theory=0,zero=False,threshold=1.0E-6,max_iter=10000):
        """evaluate structural distortion using appropriate level of SolX theory"""
        dQ = []
        f_i = self.fi()
        ### --- MCHO theory ---
        ### --- translational approximation 0 theory
        if theory == 0:
           for i in range(self.__nModes):
               dq = f_i[i] / (self.__redmass[i]*self.__freq_mcho[i]**2)
               dQ.append(dq)
           dQ = array(dQ,float64)
        ### --- New theories ---
        ### --- translational approximation 1 theory
        elif theory == 1:
           pass
        ### --- translational approximation 2 theory
        elif theory == 2:
           N = self.__nModes
           G = self.__gijj.copy()/2.
           M = self._m()
           F = self.fi()
           M1= linalg.inv(M)
           I = identity(self.__nModes,float64)
           A = tensordot(M1,G,(1,1))
           B = tensordot(A,M1,(1,0))
           a = tensordot(B,F,(2,0))
           #dQ = dot(linalg.inv(M-dot(M,a)),F)
           dQ = dot(dot(linalg.inv(I-a),M1),F)
        ### --- iterative scheme
        elif theory ==-1:
           # properties
           M = self._m()
           M1= linalg.inv(M)
           G = self.__gijj.copy()/2.
           #G.fill(0)
           F = self.fi()
           I = identity(self.__nModes,float64)
           # starting guess vector and utilities
           if zero: dq_old = zeros(self.__nModes,dtype=float64)
           else:    dq_old = self.dq(theory=2)
           #dQ, iter = get_dq_scf(M1,G,F,dq_old,threshold=threshold,maxit=max_iter)
           #print dQ
           #print iter
           if 1:
              conv = lambda dq_new,dq_old: sqrt(dot(dq_new-dq_old,dq_new-dq_old))
              GdQ  = lambda G,dq : tensordot(G,dq,(0,0))
              #GdQ  = lambda G,dq : tensordot(tensordot(G,dq,(0,0)),dq,(0,0))
              #### first iteration
              ####print shape(G),shape(dq_old)#,shape(dq_new)
              gdq = GdQ(G,dq_old)
              ####C   = dot(M1,dot(gdq,M))
              C = dot(linalg.inv((I + dot(M1,gdq))),M1)
              #####dq_new = dot(linalg.inv(M+C),F)
              dq_new = dot(C,F)
              #dq_new  = dot(M1,F+gdq)
              ##### start iterations
              ####fig = pylab.plt.figure()
              ####ax = fig.add_subplot(111)
              ####ax.plot(dq_new,'r')
              ####fig.show()
              iter = 0
              while conv(dq_new,dq_old) > threshold:
                 dq_old = dq_new.copy()
                 gdq = GdQ(G,dq_new)
                 ####C   = dot(M1,dot(gdq,M))
                 C = dot(linalg.inv((I + dot(M1,gdq))),M1)
                 dq_new = dot(C,F)
                 #dq_new  = dot(M1,F+gdq)
                 ####dq_new = dot(linalg.inv(M+C),F)
                 ####print      conv(dq_new,dq_old)
                 if iter == max_iter: break
                 iter+=1
                 ####print dq_new
                 ####ax.plot(dq_new,'r')
                 if (not iter%10000): print iter,max_iter
              ####fig.show()
              dQ = dq_new
           self.__iter = iter
        return dQ

    def __repr__(self):
        """print the output"""
        log = '\n'
        if self.__eval:
           log+= ' NORMAL MODE DISPLACEMENTS DUE TO SOLVATION\n'
           if self.__theory:
              log+= ' %s %s %s\n'%('Mode'.rjust(4),'dQ [A.U.]'.rjust(10),'dQ [A.U.]'.rjust(10))
              if self.__iter is not None:
                 log+= ' Iterations: %10d\n'%self.__iter
              for i in range(self.__nModes):
                  log+= ' %4i %10.6f %10.6f\n' % ((i+1),self.__dq[i],self.__dq_th[i])
           else:
              log+= ' %s %s\n'%('Mode'.rjust(4),'dQ [A.U.]'.rjust(10))
              for i in range(self.__nModes):
                  log+= ' %4i %10.6f\n' % ((i+1),self.__dq[i])
           log+= ' NEW HARMONIC FREQUENCIES\n'
           log+= ' %s %s\n'%('Mode'.rjust(4),'Freq [cm-1]'.rjust(10))
           for i in xrange(len(self.__freq_new)):
               log+= ' %4i %10.2f\n' % ((i+1),self.__freq_new[i])
           log+= ' APPROXIMATED FREQUENCIES\n'
           log+= str(self.__freq_approx)
        else:
           log+= 'NO EVALUATION PERFORMED\n\n'
        return str(log)
    
    # protected
    
    def _prepare_dma(self):
        fder = []
        for i in range(self.__nModes):
            dma1 = self.__fderiv[i].copy()
            if self.__ua_list is not None:
                dma1.MakeUa(self.__ua_list,change_origin=True)
            dma1.set_structure(pos=self.__solute.pos, origin=self.__solute.pos)
            fder.append(dma1)
        self.__fderiv = fder
        return
    
    def _M(self):
        """M matrix"""
        N = len(self.__mol.atoms)
        M = zeros((N*3,N*3),dtype=float64)
        for i in xrange(N):
            m = 1./sqrt(self.__mol.atoms[i].mass()*self.AmuToElectronMass)
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M
    
    def _M1(self):
        """M matrix"""
        N = len(self.__mol.atoms)
        M = zeros((N*3,N*3),dtype=float64)
        for i in xrange(N):
            m = sqrt(self.__mol.atoms[i].mass()*self.AmuToElectronMass)
            M[3*i+0,3*i+0] = m
            M[3*i+1,3*i+1] = m
            M[3*i+2,3*i+2] = m
        return M
    
    def _m(self):
        """M_iw_i^2 diagonal matrix + K_ij (B matrix)"""
        return diag(self.__redmass*self.__freq_mcho**2)
        
    def _weight(self):
        """weights derivative tensor elements
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
        #self.__gijj[-3:,-3:,-3:] = 0.0
        # reduced masses
        self.__redmass *= self.AmuToElectronMass