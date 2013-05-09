# --------------------------------------- #
#             FITTING MODULE              #
# --------------------------------------- #

from numpy import array, float64, zeros, sqrt, sum
from units import *
from dma   import *
from diff  import *
import os

class FIT(UNITS,DIFF):
    """contains usefull procedures for fitting
       of DMA and molecular multipole moments in
       """
       
    def __init__(self,freq=0,step=0,
                 dir=0):

        # number of normal modes
        self.nModes = len(freq)
        self.freq = freq
        if not dir.endswith("/"): dir+="/"
        # step of differentiation and differentiation mode (pointity)
        self.step, self.n_point, self.file_type = self.ReadStep(dir)
        # DMA set from fitting input files
        self.DMA_set, self.Fragments_set = self.ParseDMA_set(dir,self.file_type)
        # number of distributed fragments
        self.nfrags = len(self.Fragments_set[0])
        # --- [!] Calculate the derivatives!
        self.Fder, self.Sder ,self.FDip, self.SDip = self.CalcDer()
        # --- [!] Calculate IR Harmonic intensities
        self.IR_Harm_Int = self.CalcIrInt()
              
    def CalcDer(self):
        """calculates first and second (diagonal) 
           derivatives wrt normal modes
           of DMA distribution as well as
           first derivatives of molecular dipole moment
           wrt normal modes"""
           
        # calculate set of molecular dipole moments in [Bohr*ElectronCharge]
        Dipole_Moments = []
        for i,dma in enumerate(self.DMA_set):
            moment = dma[0].reshape((self.nfrags,1)) * self.Fragments_set[i] + dma[1]
            Dipole_Moments.append(sum(moment,axis=0))
        Dipole_Moments = array(Dipole_Moments)
        
        # store first and diagonal second derivatives in the lists
        Fder = [] 
        Sder = []
        FDip = []
        SDip = []
        
        # make a array of dx
        dx=[x for x in [-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5, 
                         0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 ,4.0]]
        dx=array(dx,dtype=float64)
                                     
        # fit the derivatives
        for i in range(self.nModes):
            k = 1 + 2* i* abs(self.n_point)     
            # calculate the -12A- coefficients
            A=[ DMA(nfrag=self.nfrags) for x in range(12) ]
            Ysr = DMA(nfrag=self.nfrags)
            for p in range(2*abs(self.n_point)):
                df = self.DMA_set[k+p] - self.DMA_set[0]
                x = dx[p]
                if i==0: pass
                  #print (self.DMA_set[k+p]-self.DMA_set[0])[0]
                A[0] += df * x      ### A
                A[1] += df * x**2   ### B
                A[2] += df * x**3   ### C      
                A[3] += x**2        ### D
                A[4] += x**3        ### E
                A[5] += x**4        ### F
                A[6] += x**3 /2.    ### G
                A[7] += x**4 /2.    ### H
                A[8] += x**5 /2.    ### I
                A[9] += x**4 /6.    ### J
                A[10]+= x**5 /6.    ### K
                A[11]+= x**6 /6.    ### L  
                Ysr  += self.DMA_set[k+p] 
            Ysr/=  2*abs(self.n_point)

            # calculate first derivatives
            fder_DMA = A[0]/A[3]

            fdedr_DMA = ( A[0]*( A[8] * A[10] - A[7] * A[11] ) + 
                         A[6]*( A[1] * A[11] - A[2] * A[10] ) +
                         A[9]*( A[2] * A[7]  - A[8] * A[1]  ) )/\
                       ( A[3]*( A[8] * A[10] - A[7] * A[11] ) +
                         A[6]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[5] * A[7]  - A[4] * A[8]  ) )

            # calculate second derivatives
            sder_DMA = DMA(nfrag=self.nfrags)
            sdedr_DMA = ( A[3]*( A[2] * A[10] - A[1] * A[11] ) +
                         A[0]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[1] * A[5]  - A[2] * A[4]  ) )/\
                       ( A[3]*( A[8] * A[10] - A[7] * A[11] ) +
                         A[6]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[5] * A[7]  - A[4] * A[8]  ) )
                         
            # calculate R^2 coefficients
            Serr=DMA(nfrag=self.nfrags)
            Stot=DMA(nfrag=self.nfrags)
            for p in range(2*abs(self.n_point)):
                x = dx[p]
                Serr+= (self.DMA_set[k+p] - self.DMA_set[0] - fder_DMA*x )*(self.DMA_set[k+p] - self.DMA_set[0] - fder_DMA*x )
                Stot+= (Ysr- self.DMA_set[k+p] )*(Ysr- self.DMA_set[k+p]  )
            R2 = -Serr/Stot + 1
            print "R-SQUARE COEFFICIENTS FOR MODE %d" %(i+1)
            print
            print R2    




            
            # Calculate A, B and C coefficients for the dipole moments
            A=zeros((12,3),dtype=float64)
            for p in range(2*abs(self.n_point)):
                df = Dipole_Moments[k+p] - Dipole_Moments[0]
                x = dx[p]
                A[0] += df * x      ### A
                A[1] += df * x**2   ### B
                A[2] += df * x**3   ### C   
                A[3] += x**2        ### D
                A[4] += x**3        ### E
                A[5] += x**4        ### F
                A[6] += x**3 /2.    ### G
                A[7] += x**4 /2.    ### H
                A[8] += x**5 /2.    ### I
                A[9] += x**4 /6.    ### J
                A[10]+= x**5 /6.    ### K
                A[11]+= x**6 /6.    ### L   

            # calculate first derivatives of dipole moments
            fdip     = ( A[0]*( A[8] * A[10] - A[7] * A[11] ) + 
                         A[6]*( A[1] * A[11] - A[2] * A[10] ) +
                         A[9]*( A[2] * A[7]  - A[8] * A[1]  ) )/\
                       ( A[3]*( A[8] * A[10] - A[7] * A[11] ) +
                         A[6]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[5] * A[7]  - A[4] * A[8]  ) )


            # calculate second derivatives of dipole moments
            sdip     = ( A[3]*( A[2] * A[10] - A[1] * A[11] ) +
                         A[0]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[1] * A[5]  - A[2] * A[4]  ) )/\
                       ( A[3]*( A[8] * A[10] - A[7] * A[11] ) +
                         A[6]*( A[4] * A[11] - A[5] * A[10] ) +
                         A[9]*( A[5] * A[7]  - A[4] * A[8]  ) )
                         
                                    
            # accumulate the results in the lists             
            Fder.append( fder_DMA )
            Sder.append( sder_DMA )
            FDip.append( fdip )
            SDip.append( sdip )
            
        return Fder, Sder, array( FDip ), array( SDip )