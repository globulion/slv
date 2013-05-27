﻿# --------------------------------------------------------------- #
#                       DIPOLE DERIVATIVES                        #
# --------------------------------------------------------------- #

from head import *

class DipDeriv(UNITS,FREQ):
    """dipole derivatives wrt normal mode class in helico representation"""
    
    def __init__(self,file="",L=0,step=0.006):
        FREQ.__init__(self,file)
        # step of differentiation
        self.h = step
        # derive first derivatives (helico!)
        self.fdip = self.FDeriv(Print=0,Debye=0,divide=0)
        # calculate second derivatives (helico)
        self.sdip = self.SDeriv()
        
        #print "pierwsze pochodne z fdipa"       
        #print self.fdip[7]
        #print "drugie pochodne z sdipa"
        #print self.sdip[-1][7]
        #print self.sdip[-2][7,7]
        
        
    def SDeriv(self,symm=0,Print=0):
        """ computes second derivatives of dipole moment wrt normal modes """
        A = os.listdir('.')
        A.sort()
        B = A[:]
        for file in B:
            if (not file.endswith('_.log') ) : A.remove(file)
        #for i in A:print i
        a = len(A)
        b = self.Nmodes ; N = self.Natoms
        D = zeros((a  ,3*N,3  ))  # first derivatives wrt cartesian coordinate
        for i in range(a):
            D[i] = self.DipoleDeriv(A[i])

        # make second and third derivatives wrt cartesian coordinates
        S = zeros((3*N,3*N,3  ))   # second ij      # tensor 9x9x3 dla wody
        for i in range(3*N):
            K = 4*i + 1
            S[i] = (1./12.) * ( 8.*(D[K+1] - D[K+2]) + 
                                   (D[K+3] - D[K+0]) )  / (self.h * self.AngstromToBohr )
        LT = transpose(self.L)
        #print self.h

        # transform to normal mode space
        XS = dot(dot(LT,S[:,:,0]),self.L)
        YS = dot(dot(LT,S[:,:,1]),self.L)
        ZS = dot(dot(LT,S[:,:,2]),self.L)

        E = tensordot(S,self.L,(1,0))
        E = tensordot(LT,S,(1,0))
        if symm:
           for M in [XS,YS,ZS]:
               for i in range(len(XS)):
                   for j in range(len(XS)):
                       if i>j:
                          M[i][j] = 0.5 * ( M[i][j] + M[j][i] )
                          M[j][i] = M[i][j]

        s_diag = zeros((b,3),dtype=float64)
        for i in range(b):
            s_diag[i,0] = XS[i,i]
            s_diag[i,1] = YS[i,i]
            s_diag[i,2] = ZS[i,i]
        
        sderiv = zeros((b,b,3),dtype=float64)
        sderiv[:,:,0] = XS
        sderiv[:,:,0] = YS
        sderiv[:,:,0] = ZS
        
        if Print:
             y = arange(len(XS))+1
             print " \n Second Derivatives (Sij) of Dipole Moment ( 4-point central difference, h = %f Å ) \n" % self.h
             if symm:  print "            --- symmetrized arithmetically --- "
             #if Debye: 
             #          print "                deriv units: Debye "
             else:     print "                deriv units: AU    "
             print           "                frequencies given in [cm-1]" 
             print
             print "         X -component: \n"
             PRINTV(XS,y,self.freq,y)
             print "         Y -component: \n"
             PRINTV(YS,y,self.freq,y)
             print "         Z -component: \n"
             PRINTV(ZS,y,self.freq,y)

        return sderiv,s_diag     
       
    
    
    def __repr__(self):
        log = "\n"
        return str(log)
    
