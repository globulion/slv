# --------------------------------------------------------------- #
#           MINHAENG CHO's SOLVATOCHROMIC MODEL MODULE            #
# --------------------------------------------------------------- #

from libbbg.units    import *
from numpy    import *
from libbbg.dma      import *
from dipderiv import *
from onsager  import *
import os

__all__ = ['MCHO',]
__version__ = '3.0.0'

class Stark(UNITS):
    """represents Stark tunning rate"""
    
class MCHO(UNITS):
    """contains useful procedures for calculation
    of solvatochromic shifts from M. Cho's model.
    (add citations!)"""
    
    def __init__(self,
                 ### basic variables
                 fderiv=0,sderiv=0,
                 gijj=0,
                 redmass=0,freq=0,
                 ### Cho-Onsager calculations
                 onsager=False,
                 fdip=0,
                 sdip=0,
                 sdip_full=0,
                 dipole=0,
                 giijk=0,
                 fpol=0,
                 spol=0,
                 ### Cho-Onsager including self-polarization
                 ons_pol=False,
                 cavity_radius=7.0,
                 pol_1=[0],
                 pol_2=[0],
                 step_cart=0,
                 step_mode=0,
                 L=0,
                 max_iter=20,
                 threshold=1.E-6,
                 epsilon=0,
                 iterate=False,
                 ### overall solvatochromic multipole moments
                 pol=False,
                 overall_MM=False,
                 fder_MM=0,
                 sder_MM=0,
                 MM_mode_id=0,
                 ### solvatochromic CAMMs
                 camm=False,
                 fder_CAMM=0,
                 sder_CAMM=0,
                 ### Hessian reconstruction
                 Hess=False,
                 ### EDS interaction energy derifatives
                 eds=False,
                 fEDS=0,
                 sEDS=0,
                 nmodes=0,
                 dir='./'):
                    
        self.dir = dir
        if not dir.endswith('/'): self.dir+='/'

        if not eds:self.nModes = len(fderiv)
        else:      self.nModes = nmodes
        self.nAtoms = (self.nModes + 6)/3
        self.__ReadSlvFragFile()  ### ---> fragments,atoms
        self.nfrag=len(self.fragments)
        
        self.step_cart = step_cart
        self.step_mode = step_mode
        self.pol_1     = pol_1
        self.pol_2     = pol_2
        self.L = L
        
        # first and second derivatives
        self.fderiv = fderiv
        self.sderiv = sderiv
        self.fdip = fdip
        self.sdip = sdip
        self.sdip_full = sdip_full
        self.fpol = fpol
        self.spol = spol
        self.fEDS = fEDS
        self.sEDS = sEDS
        
        # permanent dipole moment and polarizability
        self.dipole = dipole
        self.polarizability = pol_1[0]
        
        # cubic anh. constants, red. masses and frequencies
        self.gijj   = gijj
        self.giijk  = giijk
        self.redmass = redmass
        freq_cmrec = zeros(self.nModes,dtype=float64)
        freq_cmrec[:] = array(freq)
        self.freq_mcho = array(freq) * self.CmRecToHz * self.HzToAuAngFreq
        self.freq_cmrec = freq_cmrec
        
        # Cho-Onsager
        self.max_iter = max_iter
        self.threshold = threshold
        self.epsilon = epsilon
        self.iterate = iterate
        # write the report!
        self.report()
        
        # --- Calculate the parameters!!!
        if not eds:
           if not onsager: self.parameters_set = self.Parameters(fderiv=self.fderiv,
                                                              sderiv=self.sderiv)
                                                              
           # --- Calculate solute DMA object of overall molecular solvatochromic CAMMs
           if camm:
              self.parameters_CAMM = self.Parameters(fderiv=fder_CAMM,
                                                     sderiv=sder_CAMM)
                                                          
           # --- Calculate the Stark tunning rates
           if not onsager: self.stark          = self.StarkTunningRates(self.parameters_set)
        
           # --- Calculate the Cho-Onsager coefficients
           if onsager: self.onsager = self.ChoOnsagerShift(cavity_radius)
        
           # --- Calculate the Cho-Onsager+Polarisation frequency shifts
           if ons_pol:
              self.Moja(cavity_radius=cavity_radius,mode_id=MM_mode_id,
                        max_iter=self.max_iter,iterate=self.iterate,
                        threshold=self.threshold,epsilon=self.epsilon)
           
           # --- Calculate solute DMA object of overall molecular solvatochromic MMs
           if overall_MM:
              self.parameters_MM = self.SolvatochromicMultipoles(fder_MM=fder_MM,sder_MM=sder_MM,
                                                                 mode_id=MM_mode_id,fdip=fdip,sdip=sdip)
                                                              
              #p = self.Moja(cavity_radius=7.0298,mode_id=MM_mode_id)
              #p = self.Moja(cavity_radius=7.8611,mode_id=MM_mode_id)
              #p = self.Moja(cavity_radius=7.9711,mode_id=MM_mode_id)

           # --- Molecular Solvatochromic polarizability
           self.solpol = [None]
           if pol:
              self.solpol = self.SolvatochromicPolarizability(MM_mode_id,fpol,spol)
           # --- Reconstruct new Hessian matrix from DMA objects of fder and sder
           if Hess:
              print "not yet written but soon will be"
        # --- Calculate explicitly frequency shifts using interaction energy derivatives
        else:
           self.eds_shifts = self.EDS_FrequencyShifts(MM_mode_id,fEDS,sEDS)
           
    def StructuralChange(self,fderiv=0):
        """estimates structural change due to solvation"""
        
        dQ = []
        for j in range(self.nModes):
            1
        
    def Parameters(self,fderiv=0,sderiv=0):
        """calculates solvatochromic multipole
        transition coefficients"""

        parameters_set = []
        for j in range(self.nModes):
            sum = DMA(nfrag=self.nfrag)
            for i in range(self.nModes):
             #if i==j:
                ### ACCUMULATE THE SUM OVER OTHER MODES
                sum -= fderiv[i]*sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                    (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                    (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
                #
            ### ADD SDERIV CONTRIBUTION
            sum += sderiv *self.redmass[j]*self.AmuToElectronMass
            ###
            sum /= 2.0 * self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j]
            ###
            sum.pos = array(self.fragments)
            sum.origin = array(self.fragments)
            parameters_set.append(sum)

        # multiply by appropriate constants
        #for j in range(self.nModes):
        #    parameters_set[j].DMA[1] *= 1
        #    parameters_set[j].DMA[2] *= 1/3.
        #    parameters_set[j].DMA[3] *= 1/15.

        return parameters_set


    def EDS_FrequencyShifts(self,mode_id,fEDS,sEDS):
        """calculate frequency shifts from interaction energy derivatives"""
        j = mode_id
        n = fEDS.shape[1]
        a_MA = zeros(n,dtype=float64)
        a_EA = zeros(n,dtype=float64)
        ### calculate mechanical anharmonicity component, C_MA
        for i in range(self.nModes):
          #if i==j:
            ### ACCUMULATE THE SUM OVER OTHER MODES
            a_MA -= fEDS  [i]*sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                 (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                 (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
        ### calculate electronic anharmonicity component, C_EA
        a_EA += sEDS *self.redmass[j]*self.AmuToElectronMass
        ###
        a_MA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        a_EA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        
        a_tot = a_MA + a_EA
        
        return array([ a_MA, a_EA, a_tot ]) * self.HartreePerHbarToCmRec
    
    def SolvatochromicPolarizability(self,mode_id,fpol,spol):
        """calculate molecular solvatochromic polarizability for a given mode"""
        
        j = mode_id
        a_MA = zeros((3,3),dtype=float64)
        a_EA = zeros((3,3),dtype=float64)
        print j
        ### calculate mechanical anharmonicity component, C_MA
        for i in range(self.nModes):
          #if i==j:
            ### ACCUMULATE THE SUM OVER OTHER MODES
            a_MA -= fpol  [i]*sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                 (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                 (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
        ### calculate electronic anharmonicity component, C_EA
        a_EA += spol *self.redmass[j]*self.AmuToElectronMass
        ###
        a_MA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        a_EA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        
        a_tot = a_MA + a_EA
        
        return a_MA, a_EA, a_tot
                    
    def SolvatochromicMultipoles(self,fder_MM,sder_MM,mode_id,fdip,sdip):
        """overall molecular solvatochromic multipoles
        from overall molecular derivatives of MMMs.
        It contains also partitioning for:
        - mechanichal anharmonicity term (MA),
        - electronical anharmonicity term (EA).
        """
        
        j = mode_id
        C_MA = DMA(nfrag=1)
        C_EA = DMA(nfrag=1)
        d_MA = zeros(3,dtype=float64)
        d_EA = zeros(3,dtype=float64)
        
        ### calculate mechanical anharmonicity component, C_MA        
        for i in range(self.nModes):
          #if i==j:
            ### ACCUMULATE THE SUM OVER OTHER MODES
            C_MA -= fder_MM[i]*sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                 (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                 (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
            
            d_MA -= fdip  [i]*sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                 (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                 (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
            #
        ### calculate electronic anharmonicity component, C_EA
        C_EA += sder_MM *self.redmass[j]*self.AmuToElectronMass
        d_EA += sdip[mode_id] *self.redmass[j]*self.AmuToElectronMass
        ###
        C_MA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        C_EA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        d_MA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
        d_EA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )

        C = C_MA+C_EA
        d = d_MA+d_EA
        
        ### mechanical solvatochromic moments
        solvatochromic_MM_MA = DMA(nfrag=1)
        solvatochromic_MM_MA.pos = zeros((1,3),dtype=float64)
        solvatochromic_MM_MA.DMA[1] = array([d_MA])
        #solvatochromic_MM_MA.DMA[1] = array(C_MA.DMA[1])
        solvatochromic_MM_MA.DMA[2] = array(C_MA.DMA[2])
        solvatochromic_MM_MA.DMA[3] = array(C_MA.DMA[3])
        
        ### electronic solvatochromic moments
        solvatochromic_MM_EA = DMA(nfrag=1)
        solvatochromic_MM_EA.pos = zeros((1,3),dtype=float64)
        solvatochromic_MM_EA.DMA[1] = array([d_EA])
        solvatochromic_MM_EA.DMA[2] = array(C_EA.DMA[2])
        solvatochromic_MM_EA.DMA[3] = array(C_EA.DMA[3])
        
        ### total solvatochromic moments
        solvatochromic_MM = solvatochromic_MM_MA + solvatochromic_MM_EA
        
        ### label the objects
        solvatochromic_MM   .name = 'TOTAL ANHARMONICITY'
        solvatochromic_MM_MA.name = 'MECHANICAL ANHARMONICITY'
        solvatochromic_MM_EA.name = 'ELECTRONIC ANHARMONICITY'
        
        solvatochromic_MM   .set_structure(pos=self.fragments)
        solvatochromic_MM_MA.set_structure(pos=self.fragments)
        solvatochromic_MM_EA.set_structure(pos=self.fragments)

        return solvatochromic_MM_MA, solvatochromic_MM_EA, solvatochromic_MM
            
    def ChoOnsagerShift(self,cavity_radius):
        """Onsager-Cho solvatochromic model. Returns mechanical and electronic
        coefficients for each mode (add citations!)"""

        ### calculate mechanical anharmonicity component, C_MA
        coeffs = zeros((self.nModes,3),dtype=float64)
        for j in range(self.nModes):
            C_MA = 0
            C_EA = 0
            C_QF = 0
            for i in range(self.nModes):
             #if i==j:
                ### ACCUMULATE THE SUM OVER OTHER MODES
                C_MA -= dot(self.fdip[i],self.dipole)*\
                    sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                    (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                    (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))
                #
                #for k in range(self.nModes):
                #   C_QF += self.giijk[j,j,i,k]*sqrt(self.redmass[i]*self.redmass[k]*self.redmass[j]**2) /\
                #     (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2*\
                #      self.redmass[k]*self.AmuToElectronMass*self.freq_mcho[k]**2)*\
                #      dot(self.fdip[i],self.dipole) * dot(self.fdip[k],self.dipole)*\
                #      sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                #      sqrt(self.redmass[k]*self.AmuToElectronMass)
            #C_QF /= 2
            ### ADD SDERIV CONTRIBUTION
            C_EA += (dot(self.sdip[j],self.dipole) + 
                     dot(self.fdip[j],self.fdip[j]))\
                    *self.redmass[j]*self.AmuToElectronMass
            ###
            C_MA /= (self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
            C_EA /= (self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
            C_QF /= (self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
            ### save the coefficients
            coeffs[j,0] = C_MA
            coeffs[j,1] = C_EA
            coeffs[j,2] = C_QF
  
        #print "\n COEFFICIENTS\n"
        #print coeffs[:,1]
        #print sqrt(sum(self.dipole**2)) * self.BohrElectronToDebye
        #print self.dipole
                           
        return coeffs
    
    def Moja(self,cavity_radius,mode_id,iterate=True,epsilon=0,threshold=0,max_iter=0):
        """Cho-Onsager+Polarisation frequency shifts. Differentiates FI tensor
        (see the reference) wrt normal modes for each dielectric constant value
        separately and calculates frequency shift. Assumes no structural changes
        if <iterate> is set to <False>"""
        
        ChoOnsager = ONSAGER(gijj=self.gijj,redmass=self.redmass,freq=self.freq_cmrec,
                 step_cart=self.step_cart,step_mode=self.step_mode,
                 cavity_radius=cavity_radius,mode_id=mode_id,
                 L=self.L,fdip=self.fdip,sdip=self.sdip,pol_1=self.pol_1,pol_2=self.pol_2,
                 dipole=self.dipole,calculate=True)
                 
        if (iterate and epsilon == 0):
           print " ------------------------------ STARTING SCF CALCULATIONS ------------------------------ "
           print "                                max_iter : %i" % max_iter
           print "                                threshold: %.4E" % threshold
           print
           for i,eps in enumerate(ChoOnsager.epsilons):
               solver = BlasiakIterator(epsilon=eps,
                                        mode_id=mode_id, redmass=self.redmass,
                                        dipole=self.dipole, fdip=self.fdip, sdip=self.sdip_full,
                                        polarizability=self.polarizability,
                                        fpol=ChoOnsager.get_pol_fder(),
                                        spol=ChoOnsager.get_pol_sder(),
                                        fFi=ChoOnsager.get_FI_fder(epsilon),
                                        sFi=ChoOnsager.get_FI_sder(epsilon),
                                        freq=self.freq_cmrec,
                                        cavity_radius=cavity_radius,
                                        Fi=ChoOnsager.get_FI(epsilon),
                                        dipole_moment=ChoOnsager.dipole_moments[i],
                                        solvation_energy=ChoOnsager.solvation_energies[i],
                                        frequency_shift=ChoOnsager.frequency_shifts[i],
                                        gijj=self.gijj)
               solver.iterate(max_iter=max_iter,threshold=threshold)
               # substitute converged values in ONSAGER instance
               ChoOnsager.dipole_moments[i],\
               ChoOnsager.solvation_energies[i],\
               ChoOnsager.frequency_shifts[i],\
               ChoOnsager.Onsager_fields[i] = solver.get_result()
               ChoOnsager.dipole_moments[i]*= self.BohrElectronToDebye
               
               
                 
        print ChoOnsager
    
    def Moja_old(self,cavity_radius,mode_id,iterate=True):
        """Cho-Onsager+Polarisation frequency shifts. Differentiates FI tensor
        (see the reference) wrt normal modes for each dielectric constant value
        separately and calculates frequency shift. Assumes no structural changes"""
        
        s1 = 1./(self.step_cart * self.AngstromToBohr)
        s2 = 1./self.step_mode**2
        
        data = []
        onsager_factors = []
        dipole_moments = []
        solvation_energies = []
        # iterate over dielectric constants
        for epsilon in [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,3.0,
                        4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,
                        18.0,20.0,23.0,26.0,29.0,33.0,37.0,45.0,50.0,
                        58.0,65.0,78.0,90.0,100.0]:
            
            f = 2*(epsilon-1)/(cavity_radius**3 * (2*epsilon+1))
        
            # cartesian first derivatives of FI tensor
            B_cart = zeros((1+4*self.nAtoms*3,3,3),dtype=float64)
            C_cart = zeros((1+4*self.nAtoms*3,3,3),dtype=float64)
            
            for i in range(1+4*self.nAtoms*3):
                b = linalg.inv(identity(3) - f*self.pol_1[i])
                B_cart[i] = b.copy()
                C_cart[i] = dot(transpose(b),dot(self.pol_1[i],b))
            
            # instantaneous dipole moment
            inst_dipole = dot( (identity(3) + f*dot(self.pol_1[0],B_cart[0])) , self.dipole) 
            inst_dipole*= self.BohrElectronToDebye
            dipole_moments.append(inst_dipole)
            
            #
            FI_cart = B_cart*f +  C_cart *f**2
            FI = FI_cart[0]
            
            # instantaneous solvation energy
            inst_e_sol = -0.5 * tensordot(tensordot(FI,self.dipole,(0,0)),self.dipole,(0,0))
            solvation_energies.append(inst_e_sol)

            first_der_FI_cart = []
            for i in range(self.nAtoms*3):
                K = 1 + 4*i 
                ### FI(epsilon,a) tensor
                fd   = (1./12.) *  (\
                           (FI_cart[K+3] - FI_cart[K+0] ) \
                     + 8*  (FI_cart[K+1] - FI_cart[K+2] ))\
                     * s1
                first_der_FI_cart.append( fd ) 
            first_der_FI_cart = array(first_der_FI_cart)
                
            # wrt-mode second derivatives of FI tensor
            B_mode = zeros((5,3,3),dtype=float64)
            C_mode = zeros((5,3,3),dtype=float64)
            
            for i in range(5):
                b = linalg.inv(identity(3) - f*self.pol_2[i])
                B_mode[i] = b.copy()
                C_mode[i] = dot(transpose(b),dot(self.pol_2[i],b))
                                
            FI_mode = B_mode*f + C_mode*f**2
            
            second_der_FI_mode = (1./12.) *  (\
                 -     (FI_mode[1+3] + FI_mode[1+0] ) \
                 +16*  (FI_mode[1+1] + FI_mode[1+2] ) \
                 -30*   FI_mode[0  ] )                \
                 *s2
                 
            ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
            first_der_FI_mode = tensordot(self.L,first_der_FI_cart,(0,0))
            
            # calculate frequency shifts!!!
            j = mode_id
            dw_MA = 0
            dw_EA = 0
        
            ### calculate mechanical anharmonicity component, C_MA    
            #print second_der_FI_mode  
            for i in range(self.nModes):
             #if i==j:
                ### ACCUMULATE THE SUM OVER OTHER MODES
                dw_MA -= sqrt(self.redmass[i]*self.AmuToElectronMass)*\
                     (self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)/
                     (self.redmass[i]*self.AmuToElectronMass*self.freq_mcho[i]**2))*\
                     (
                      tensordot(tensordot(first_der_FI_mode[i],self.dipole,(0,0)),
                                self.dipole,(0,0)) +
                      tensordot(tensordot(FI,self.fdip[i],(0,0)),
                                self.dipole,(0,0)) +     
                      tensordot(tensordot(FI,self.dipole,(0,0)),
                                self.fdip[i],(0,0)) \
                     )
                     

                #
            ### calculate electronic anharmonicity component, C_EA
            dw_EA += self.redmass[j]*self.AmuToElectronMass*\
                  (
                      tensordot(tensordot(second_der_FI_mode,self.dipole,(0,0)),
                                self.dipole,(0,0)) +
                      tensordot(tensordot(FI,self.sdip[j],(0,0)),
                                self.dipole,(0,0)) +                   
                      tensordot(tensordot(FI,self.dipole,(0,0)),
                                self.sdip[j],(0,0)) +
                  2 * tensordot(tensordot(first_der_FI_mode[j],self.fdip[j],(0,0)),
                                self.dipole,(0,0)) +    
                  2 * tensordot(tensordot(first_der_FI_mode[j],self.dipole, (0,0)),
                                self.fdip[j],(0,0)) + 
                  2 * tensordot(tensordot(FI,self.fdip[j],(0,0)),
                                self.fdip[j],(0,0)) \
                   )
            

            ###
            dw_MA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )
            dw_EA /= ( 2*self.redmass[j]*self.AmuToElectronMass * self.freq_mcho[j] )

            dw = -1.*(dw_MA+dw_EA)/2.
            ### switch from au to cm-1
            dw*= self.HartreePerHbarToCmRec
            
            ### collect the result
            data.append( dw )
            onsager_factors.append( f * cavity_radius**3 )
            
        data = array(data)
        onsager_factors = array( onsager_factors )
        dipole_moments = array( dipole_moments )
        solvation_energies = array( solvation_energies )
#Dipole Moment [Debye]
        N = 100
        print
        print " "+"-"*N
        print "CHO-ONSAGER+POLARISATION FREQUENCY SHIFT".center(N)
        print " "+"-"*N
        print " %12s %12s %12s %12s %12s %12s %12s %12s" % ('Ons. factor'.rjust(12),
                                'Shift'.rjust(12),
                                'Frequency'.rjust(12),
                                'X'.rjust(12),
                                'Y'.rjust(12),
                                'Z'.rjust(12),
                                'AVER'.rjust(12),
                                'E_sol'.rjust(12))
        print " %12s %12s %12s %12s %12s %12s %12s %12s" % ('[---]'.rjust(12),
                                '[cm-1]'.rjust(12),
                                '[cm-1]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Hartree]'.rjust(12)) 
        print " "+"-"*N
        d = self.dipole * self.BohrElectronToDebye
        print     " %12.4f %12.2f %12.2f %12.4f %12.4f %12.4f %12.4f %12.6f"\
            % (0, 0, self.freq_cmrec[mode_id],d[0], d[1], d[2], sqrt(sum(d**2)),0)
        for i,j in enumerate(data):
            print " %12.4f %12.2f %12.2f %12.4f %12.4f %12.4f %12.4f %12.6f" % (onsager_factors[i], j,
                                         j+self.freq_cmrec[mode_id],
                                         dipole_moments[i,0],
                                         dipole_moments[i,1],
                                         dipole_moments[i,2],
                                         sqrt(sum(dipole_moments[i]**2)),
                                         solvation_energies[i])
        print " "+"-"*N
            
        return data, onsager_factors

    
    def StarkTunningRates(self,parameters_set):
        """calculates Stark Tunning Rates for each mode.
        Units: cm-1/(MV/cm)"""
        stark = []
        for mode, parameter in enumerate(parameters_set):
            stark_mode = zeros(3,dtype=float64)
            for fragment,transition_charge in enumerate(parameter[0]):
                stark_mode += transition_charge * self.fragments[fragment]
            for fragment,transition_dipole in enumerate(parameter[1]):
                stark_mode += transition_dipole
            stark.append(stark_mode) 
        return array(stark)  * 42.6810
    
    def Update(self,parameters_set):
        """Updates Stark Tunning rates for given parameters object"""
        self.stark = self.StarkTunningRates(parameters_set)
    
    def __ReadSlvFragFile(self):
        """reads slv fragment file.
        The coordinates are in AU!
        """
        
        data = open(self.dir+'slv.frags','r').readlines()
        fragments = []
        atoms = []
        for fragment in data:
            fragments.append(array(fragment.split()[1:-1],dtype=float64))
            atoms.append(Atom(fragment.split()[-1]))
            
        self.atoms =     atoms
        self.fragments = array(fragments)
    
    def __repr__(self):
        """prints the parameters and Stark tunning rates
        in a standard fashion"""
        
        # PARAMETERS FOR EACH MODE (UP TO TRANSITION OCTUPOLES)
        log = "\n" 
        for mode, parameters in enumerate(self.parameters_set):
            log+= " ================ MODE %d:%8.2f [cm-1] ================ \n" %\
                           (mode+1,self.freq_cmrec[mode])
            log+= "\n"
            log+= repr(parameters)
            
        # STARK TUNNING RATES FOR EACH MODE 
        log+= "\n"
        log+= " ====== STARK TUNNING RATES ======\n"
        log+= "\n"                
        log+= "%10s %10s %10s %10s %10s %20s\n" %\
           ("Freq".rjust(10),"[cm-1]".rjust(10),
            "X".rjust(10),"Y".rjust(10),"Z".rjust(10),
            "AVER [cm-1/(MV/cm)]".rjust(20))
        for mode, rate in enumerate(self.stark):
            avg = sqrt(sum(rate**2))
            log+= "%10d %10.2f %10.4f %10.4f %10.4f %10.4f\n" %\
            (mode+1,self.freq_cmrec[mode],
             rate[0],rate[1],rate[2],avg)
             
        return str(log)

    def __getitem__(self,index): 
        """returns the parameters for a selected mode"""
        return self.parameters_set[index]
    
    def __len__(self):
        """size of the parameter set"""
        return len(self.parameters_set)

    def report(self,yes=True):
        """make a report from hippi dipsy. It creates a list of
        reduced masses, cubic force constants. The file name is
        < report.txt > and is written in the working directory"""
                   
        if yes: 
         out = open('report.txt','w')   

         print >> out, "HARMONIC FREQUENCIES [cm-1]\n"        
         for i in range(   self.nModes): 
             print >> out, "%d %14.3f"%(i+1,self.freq_cmrec[i] )
                         
         print >> out, "REDUCED MASSES [AMU]\n"        
         for i in range(   self.nModes): 
             print >> out, "%d %14.4f"%(i+1,self.redmass[i] ) 
             
         print >> out, "CUBIC ANHARMONIC CONSTANTS [Eh/ao^3]\n"
         print >> out, "DIAGONALS"
         for i in range(   self.nModes): 
             print >> out, "%d %d %d %14.6f"%(i+1,i+1,i+1,
                self.gijj[i,i,i]*sqrt(self.redmass[i]*self.redmass[i]**2)) 
                
         print >> out, "SEMI-OFFDIAGONALS"
         print " WARNING: SLV_MODE_ID for report generation is set to %d (Python convention)" % int(os.environ['SLV_MODE_ID'])
         for i in range(   self.nModes): 
             for j in range(   self.nModes):
                if j==int(os.environ['SLV_MODE_ID']): 
                 print >> out, "%d %d %d %14.6f"%(i+1,j+1,j+1,
                self.gijj[i,j,j]*sqrt(self.redmass[i]*self.redmass[j]**2)) 
        return
