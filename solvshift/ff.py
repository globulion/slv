# --------------------------------------- #
#         DIFFERENTIATION MODULE          #
# --------------------------------------- #

from numpy import array, float64, zeros, sqrt, sum, shape
from units import *
from dma   import *
from utilities  import *
import os

__all__ = ['FF',]
__version = '2.0.0'

class FF:
   """contains finite field variables for 
      differentiation over normal modes"""
      
   # default step of cartesian differentiation (in Angstroms!!!)
   k = 0.005
   # default step of normal mode differentiation 
   s = 0.010
   # default pointity
   n_point = 5
   # cartesian 5-point displacements in Angstroms!!!
   DisplCart=[#[0.000,0.000,0.000],   # D 0
              [+2*k ,0.000,0.000],   # D 1
              [+k   ,0.000,0.000],   # D 2
              [-k   ,0.000,0.000],   # D 3
              [-2*k, 0.000,0.000],   # D 4
              [0.000,+2*k ,0.000],   # D 5
              [0.000,+k   ,0.000],   # D 6
              [0.000,-k   ,0.000],   # D 7
              [0.000,-2*k ,0.000],   # D 8
              [0.000,0.000,+2*k ],   # D 9
              [0.000,0.000,+k   ],   # D10
              [0.000,0.000,-k   ],   # D11
              [0.000,0.000,-2*k ]]   # D12
   DisplCart = array(DisplCart) 
   
   @staticmethod
   def setStep(step,n_point=5):
       """actualize step and cartesian displacement matrix"""
       FF.k=step
       FF.n_point=n_point
       k=step
       DisplCart=[#[0.000,0.000,0.000],   # D 0
                  [+2*k ,0.000,0.000],   # D 1
                  [+k   ,0.000,0.000],   # D 2
                  [-k   ,0.000,0.000],   # D 3
                  [-2*k, 0.000,0.000],   # D 4
                  [0.000,+2*k ,0.000],   # D 5
                  [0.000,+k   ,0.000],   # D 6
                  [0.000,-k   ,0.000],   # D 7
                  [0.000,-2*k ,0.000],   # D 8
                  [0.000,0.000,+2*k ],   # D 9
                  [0.000,0.000,+k   ],   # D10
                  [0.000,0.000,-k   ],   # D11
                  [0.000,0.000,-2*k ]]   # D12
       FF.DisplCart = array(DisplCart) 
       return
  
   @staticmethod
   def setDispl(L,redmass,mode_id):
       """cartesian displacement matrices. 
       All values in atomic units !"""
       
       displacements = []
       Cart = []
       nAtoms = len(L)/3
       # --- make displacement matrices for each mode
       # add zero-th displacement 
       displacements.append( zeros(len(L),dtype=float64).reshape(nAtoms,3) )
       
       # calculate cartesian norms
       for i in range(len(L[0])):
           disp_mat = zeros(len(L),dtype=float64)
           disp=zeros(len(L),dtype=float64)
           for j in range(3*nAtoms):
               disp[j] = L[j,i]             ### gvib-like mass-weighted matrix 
           Cart.append( disp.reshape(nAtoms,3) )
   
       # Calculate the normalized norms of cartesian displacements
       Norms = sum( sqrt( sum(array(Cart)**2,axis=2) ) , axis=1)
       OverallDisplacements = sqrt( sum(array(Cart)**2,axis=2) )[mode_id] 
       Min_disp = min(Norms)
       Max_disp = max(Norms)
       Norms = Norms/min(Norms)           

       for i in range(len(L[0])):
         if i == mode_id:
           disp_mat = zeros(len(L),dtype=float64)
           disp=zeros(len(L),dtype=float64)
           for j in range(3*nAtoms):
               disp_mat[j] = L[j,i] * FF.k 
               disp[j]     = L[j,i]                   ### gvib-like mass-weighted matrix 
           if   FF.n_point==3 : T = [ 1,-1] 
           elif FF.n_point==5 : T = [ 2, 1,-1,-2] 
           elif FF.n_point==7 : T = [ 3, 2, 1,-1,-2,-3]            
           elif FF.n_point==9 : T = [ 4, 3, 2, 1,-1,-2,-3,-4] 
           elif FF.n_point==11: T = [ 5, 4, 3, 2, 1,-1,-2,-3,-4,-5]
           for d in T:
               displacements.append( d * disp_mat.reshape(nAtoms,3) )

        
       # --- cartesian displacements 
       FF.displacements = array(displacements) 
       FF.Cart = array(Cart)
       FF.Norms = Norms
       FF.Min_disp = Min_disp
       FF.Max_disp = Max_disp
       FF.OverallDisplacements = OverallDisplacements
       print "\n NORMALIZED NORMS OF MOTIONS [DIMENSIONLESS]\n"
       PRINT(FF.Norms)
       print " Minimum total displacement in Å: %f   " % (Min_disp*UNITS.BohrToAngstrom)
       print " Maximum total displacement in Å: %f \n" % (Max_disp*UNITS.BohrToAngstrom)
       
       full_output = 0
       if full_output:
          print "\n STANDARD AMPLITUDES OF MOTIONS [ANGSTROM]\n"
          for i,d in enumerate(FF.Cart): 
            if i==MODE_ID:
              print "MODE %d" % (i+1)
              PRINTL( d*UNITS.BohrToAngstrom )
       FF.redmass = redmass
       FF.L = L
       FF.nAtoms = nAtoms
       FF.nModes = nAtoms * 3 - 6
