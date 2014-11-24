# --------------------------------------------------------------------- #
#           ONSAGER-CHO MODEL FOR FREQUENCY SHIFT PREDICTION            #
# --------------------------------------------------------------------- #

from libbbg.units import *
from numpy import *
from libbbg.dma   import *

__all__ = ['ONSAGER','BlasiakIterator',]
__version__ = '2.0.1'

class BlasiakIterator(object,UNITS):
    """represents iterator for converging
    solvation free energy, Onsager field and dipole moment
    for a given epsilon value"""
    
    def __init__(self,
                 epsilon=50.0, mode_id=0, redmass=0,
                 dipole=0, fdip=0, sdip=0, polarizability=0, fpol=0, spol=0,
                 fFi=0, sFi=0, dipole_moment=0, freq=0,
                 solvation_energy=0, frequency_shift=0, cavity_radius=0, Fi=0,
                 gijj=0):

        self.__epsilon = epsilon
        self.__mode_id = mode_id
        self.__cavity_radius = cavity_radius
        self.__nModes = len(fdip)
        self.__redmass = redmass
        self.__freq_mcho = freq * self.CmRecToHz * self.HzToAuAngFreq
        self.__dipole = dipole
        self.__fdip = fdip.copy()
        self.__sdip = sdip.copy()
        self.__polarizability = polarizability
        self.__fpol = fpol.copy()
        self.__spol = spol.copy()
        self.__Fi = Fi
        self.__fFi = fFi
        self.__sFi = sFi
        self.__gijj = gijj
        self.__onsager_factor = 2*(epsilon-1) / (2*epsilon+1)
        self.__g = self.__onsager_factor /cavity_radius**3
        self.dipole_moment = dipole_moment
        self.solvation_energy = solvation_energy
        self.frequency_shift = frequency_shift
        
        ### [1] re-mass-weight the derivatives
        self.__re_mass_weight_derivatives()
        ### [2] start!
        self.__Start()
        
    def iterate(self,threshold=1.e-14,max_iter=80):
        """perform SCF calculations for a given epsilon
        and Esol-threshold controlling parameter"""
        self.iter = 1
        conv = abs(self.__Esol_new - self.__Esol_old)
        #print self.__Esol_old
        #print self.__dipole_old
        #print self.__Esol_new
        # START THE SCF ITERATIONS!
        #print " ------------- STARTING SCF CALCULATIONS ------------- "
        #print "               max_iter : %i" % max_iter
        #print "               threshold: %.4E" % threshold
        #print
        while conv > threshold:
            ### make the old/new iterates
            self.__dQ_old = self.__dQ_new.copy()
            self.__dipole_old = self.__dipole_new.copy()
            self.__Fons_old = self.__Fons_new.copy()
            self.__fdip_old = self.__fdip_new.copy()
            self.__Esol_old = self.__Esol_new
            
            ### Onsager Field derivatives!
            A = self.__A + tensordot(self.__fpol,self.__dQ_old,(0,0))
            Ainv = linalg.inv(A)
            fA = self.__fpol #+ self.__spol * self.__dQ_old[:,newaxis,newaxis]
            Mu = self.__dipole + tensordot(self.__fdip,self.__dQ_old,(0,0))
            self.__fFons_new = self.__fdip + tensordot(self.__sdip,self.__dQ_old,(0,0))
            #self.__fFons_new-= dot(dot(fA,Ainv),Mu)
            self.__fFons_new-= tensordot(fA,dot(Ainv,Mu),(2,0))
            self.__fFons_new =-tensordot(Ainv,self.__fFons_new,(0,1))
            ### --- neglect all the field derivatives!
            #self.__fFons_new = zeros((3,self.__nModes),dtype=float64)
            
            ### displacements
            self.__dQ_new = tensordot(self.__fdip_old,self.__Fons_old,(1,0))
            #print self.__dipole_old.shape,self.__fFons_new.shape
            self.__dQ_new+= tensordot(self.__dipole_old,self.__fFons_new,(0,0))
            self.__dQ_new/=2.
            for i in xrange(self.__nModes):
                self.__dQ_new[i]/= self.__redmass[i]*self.AmuToElectronMass * self.__freq_mcho[i]**2
            #print self.__dQ_new
            #print self.__dQ_old
            
            ### Onsager Field
            A = self.__A + tensordot(self.__fpol,self.__dQ_new,(0,0))
            A = linalg.inv(A)
            
            Mu = self.__dipole + tensordot(self.__fdip,self.__dQ_new,(0,0))
            self.__Fons_new =-dot(A,Mu)
            
            ### dipole moment
            self.__dipole_new = self.__dipole + dot(self.__polarizability,self.__Fons_new)
            self.__dipole_new+= tensordot(self.__fdip,self.__dQ_new,(0,0))
            self.__dipole_new+= tensordot(tensordot(self.__fpol,self.__dQ_new,(0,0)),self.__Fons_new,(0,0))
            
            self.__dip_aver_new  = sqrt(sum(self.__dipole_new**2))* self.BohrElectronToDebye
            
            ### Solvation Free Energy
            self.__Esol_new = -0.5* dot(self.__dipole_new,self.__Fons_new)
            #print self.__Esol_new
            conv = abs(self.__Esol_new - self.__Esol_old)
            
            ### Frequency Shift
            self.__freq_shift =2* tensordot(self.__gijj,self.__dQ_new,(0,0))[self.__mode_id,self.__mode_id,]
            self.__freq_shift-= tensordot(diagonal(self.__sdip,axis1=0,axis2=1),self.__Fons_new,(0,0))[self.__mode_id]
            self.__freq_shift-=2*dot( dot(self.__fpol[self.__mode_id],self.__fFons_new[:,self.__mode_id]) , self.__Fons_new)
            self.__freq_shift-=2*dot( self.__fdip_old[self.__mode_id], self.__fFons_new[:,self.__mode_id])
            #self.__freq_shift-=  dot( dot(self.__spol[self.__mode_id],self.__Fons_new),self.__Fons_new )
            
            self.__freq_shift/=4*self.__redmass[self.__mode_id]*self.AmuToElectronMass * self.__freq_mcho[self.__mode_id]
            self.__freq_shift *= self.HartreePerHbarToCmRec
            
            ### check if converged and proceed if not
            if self.iter>max_iter: 
               print " * Epsilon %10.2f: Solvation Free Energy Not Converged in --- %10i --- Iterations!" % (self.__epsilon,self.iter,)
               break
            if conv <= threshold: print " * Epsilon %10.2f: Solvation Free Energy Converged in --- %10i --- Iterations!" % (self.__epsilon,self.iter,)
            self.iter+=1
            
            ### dipole moment derivatives
            self.__fdip_new = self.__fdip.copy()
            self.__fdip_new+= tensordot(self.__fpol,self.__Fons_new,(2,0))
            self.__fdip_new+= tensordot(self.__sdip,self.__dQ_new,(0,0))
            ### --- polarizability spol correction
            #k = tensordot(self.__spol,self.__Fons_new,(1,0))
            #for i in xrange(self.__nModes):
            #    k[i]*= self.__dQ_new[i]
            #self.__fdip_new+= k
            ### --- field derivative correction
            self.__fdip_new+= tensordot(self.__fFons_new,self.__polarizability,(0,0))
            ggg = tensordot(self.__fpol,self.__dQ_new,(0,0))
            self.__fdip_new+= tensordot(self.__fFons_new,ggg,(0,0))
            
            
            
        #print self.__dipole_new  * self.BohrElectronToDebye
        #print self.__Fons_new
            
    def get_result(self):
        """withdraw converged data"""
        return self.__dipole_new, self.__Esol_new, self.__freq_shift, self.__Fons_new
            
    def __Start(self):
        """calculate properties as a starting point
        of iterations"""
        ### Onsager Field
        self.__A = self.__polarizability  - identity(3)/ self.__g
        dQ_Start = self.get_dQ_Start()
        A = self.__A + tensordot(self.__fpol,dQ_Start,(0,0))
        A = linalg.inv(A)
        
        Mu = self.__dipole + tensordot(self.__fdip,dQ_Start,(0,0))
        Fons = -dot(A,Mu)
        
        ### Dipole Moment
        dipole = self.__dipole + dot(self.__polarizability,Fons)
        dipole+= tensordot(self.__fdip,dQ_Start,(0,0))
        dipole+= tensordot(tensordot(self.__fpol,dQ_Start,(0,0)),Fons,(0,0))
        dipole_aver_Debye = sqrt(sum(dipole**2))  * self.BohrElectronToDebye
    
        ### Onsager Field derivatives!
        A = self.__A + tensordot(self.__fpol,dQ_Start,(0,0))
        Ainv = linalg.inv(A)
        fA = self.__fpol #+ self.__spol * self.__dQ_old[:,newaxis,newaxis]
        Mu = self.__dipole + tensordot(self.__fdip,dQ_Start,(0,0))
        self.__fFons_new = self.__fdip + tensordot(self.__sdip,dQ_Start,(0,0))
        #self.__fFons_new-= dot(dot(fA,Ainv),Mu)
        self.__fFons_new-= tensordot(fA,dot(Ainv,Mu),(2,0))
        self.__fFons_new =-tensordot(Ainv,self.__fFons_new,(0,1))
        
        ### Solvation Free Energy
        Esol = -0.5* dot(dipole,Fons)
        #print dipole*self.BohrElectronToDebye
        #print dipole_aver_Debye
        #print Fons
        
        ### set the start variables!
        self.__dQ_new   = dQ_Start                                 # [ Bohr ]
        self.__Fons_new = Fons                                     # [ A.U. ]
        self.__dipole_new = dipole                                 # [ A.U. ]
        self.__dip_aver_new = dipole_aver_Debye                    # [ Debye ]
        self.__Esol_new = Esol                                     # [ A.U. ]
        
        ### dipole derivatives !!!
        self.__fdip_new = self.__fdip.copy()                       # [ A.U. ]
        self.__fdip_new+= tensordot(self.__fpol,self.__Fons_new,(2,0))
        self.__fdip_new+= tensordot(self.__sdip,self.__dQ_new,(0,0))
        ### --- polarizability spol correction
        #k = tensordot(self.__spol,self.__Fons_new,(1,0))
        #for i in xrange(self.__nModes):
        #    k[i]*= self.__dQ_new[i]
        #self.__fdip_new+= k
        ### --- field derivative correction
        self.__fdip_new+= tensordot(self.__fFons_new,self.__polarizability,(0,0))
        ggg = tensordot(self.__fpol,self.__dQ_new,(0,0))
        self.__fdip_new+= tensordot(self.__fFons_new,ggg,(0,0))
        
        ### set the old variables!
        self.__dQ_old   = zeros(self.__nModes,dtype=float64)       # [ Bohr ]
        self.__Fons_old = dot(self.__g*(linalg.inv(identity(3) \
                        - self.__g*self.__polarizability)),
                                   self.__dipole)                  # [ A.U. ]
        self.__dipole_old = self.dipole_moment                     # [ A.U. ]
        self.__dip_aver_old = sqrt(sum(self.dipole_moment**2))\
                            * self.BohrElectronToDebye             # [ Debye ]
        self.__Esol_old = self.solvation_energy                    # [ A.U. ]
        
    def get_dQ_Start(self):
        """calculate dQ_i as a starting point"""
        dQ = tensordot(tensordot(self.__fFi,self.__dipole,(2,0) ),self.__dipole,(1,0))
        dQ+= tensordot(tensordot(self.__Fi ,self.__fdip  ,(1,1) ),self.__dipole,(0,0))
        dQ+= tensordot(tensordot(self.__Fi ,self.__dipole,(1,0) ),self.__fdip  ,(0,1))
        dQ/=2.
        for i in xrange(self.__nModes):
            dQ[i]/= self.__redmass[i]*self.AmuToElectronMass * self.__freq_mcho[i]**2
        return dQ
    
    def __re_mass_weight_derivatives(self):
        """turns the dimension of the derivatives from [x/(Bohr*ElMass^1/2)]
        to [x/Bohr]"""
        # first of dipole moments
        temp = sqrt(self.__redmass*self.AmuToElectronMass)[:,newaxis]
        self.__fdip = temp * self.__fdip
        # first of Fi tensor
        temp = sqrt(self.__redmass*self.AmuToElectronMass)[:,newaxis,newaxis]
        self.__fFi = temp * self.__fFi
        # second of polarizability tensor
        self.__spol = temp**2 * self.__spol
        # first of polarizability tensor
        self.__fpol = temp * self.__fpol
        # second of dipole moments
        self.__sdip = temp * self.__sdip
        temp = sqrt(self.__redmass*self.AmuToElectronMass)[newaxis,:,newaxis]
        self.__sdip = temp * self.__sdip
        # cubic anharmonic constants
        temp = sqrt(self.__redmass)[:,newaxis,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,:,newaxis,]
        self.__gijj = temp * self.__gijj
        temp = sqrt(self.__redmass)[newaxis,newaxis,:,]
        self.__gijj = temp * self.__gijj
    
    def __repr__(self):
        """print the iteration process and the final data"""
        pass
        
        
class ONSAGER(UNITS):
    """contains useful procedures for calculation
    of solvatochromic shifts using hybrid Onsager-Cho solvatochromic model.
    (add citations!)"""
    
    def __init__(self,
                 gijj=0, redmass=0, freq=0,
                 step_cart=0, step_mode=0,
                 cavity_radius=0, mode_id=0,
                 L=0, fdip=0, sdip=0, pol_1=0, pol_2=0,
                 dipole=0, calculate=False, epsilon=None): 
                    
        self.cavity_radius = cavity_radius
        self.nModes = len(fdip)
        self.nAtoms = (self.nModes + 6)/3 
        self.fdip = fdip
        self.sdip = sdip
        self.pol_1     = pol_1
        self.pol_2     = pol_2
        self.L = L
        self.dipole = dipole
        self.pol = pol_1[0]
        self.mode_id = mode_id
        self.step_cart = step_cart
        self.step_mode = step_mode
        self.gijj   = gijj
        self.redmass = redmass
        self.freq_mcho = freq * self.CmRecToHz * self.HzToAuAngFreq
        self.freq_cmrec = freq.copy()
        self.epsilons = [ 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0,
                          4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,12.0,14.0,16.0,17.0,
                         18.0,20.0,23.0,26.0,29.0,33.0,37.0,45.0,50.0,
                         58.0,60.0,65.0,78.0,90.0,100.0]
                    
        if epsilon is not None:
           pass
                    
        if calculate:
           self.onsager_factors, self.dipole_moments, \
           self.solvation_energies, self.frequency_shifts,\
           self.Onsager_fields = self.Calculate_Approx()
 
 
    def Calculate_Approx(self):
        """calculate approximate properties"""
        onsager_factors = []
        dipole_moments = []
        solvation_energies = []
        frequency_shifts = []
        Onsager_fields = []
        
        for epsilon in self.epsilons:
            onsager_factors   .append( self.get_Onsager_factors (epsilon) )
            dipole_moments    .append( self.get_Dipole_Approx   (epsilon) )
            solvation_energies.append( self.get_Esol_Approx     (epsilon) )
            frequency_shifts  .append( self.get_Shift_Approx    (epsilon) )
            Onsager_fields    .append( self.get_Fons_Approx     (epsilon) )
            
        onsager_factors     =array( onsager_factors    ,dtype=float64 )
        dipole_moments      =array( dipole_moments     ,dtype=float64 )
        solvation_energies  =array( solvation_energies ,dtype=float64 )
        frequency_shifts    =array( frequency_shifts   ,dtype=float64 )
        Onsager_fields      =array( Onsager_fields     ,dtype=float64 )
        
        return onsager_factors, dipole_moments, solvation_energies, frequency_shifts, Onsager_fields 
    
    def get_Shift_Approx(self,epsilon):
        """calculate approximated frequency shift for a given epsilon 
        by neglecting structural distortions"""
        # calculate frequency shifts!!!
        j = self.mode_id
        dw_MA = 0
        dw_EA = 0
        
        FI = self.get_FI(epsilon)
        first_der_FI_mode = self.get_FI_fder(epsilon)
        second_der_FI_mode= self.get_FI_sder(epsilon)
        
        ### calculate mechanical anharmonicity component, C_MA    
        #print second_der_FI_mode  
        for i in xrange(self.nModes):
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
        
        return dw

    def get_Onsager_factors(self,epsilon):
        """calculate Onsager factors"""
        Onsager_factor = 2*(epsilon-1) / (2*epsilon+1)
        return Onsager_factor
    
    def get_Fons_Approx(self,epsilon):
        """calculate Onsager field for a given epsilon. Returned in AU unit"""
        f = 2*(epsilon-1)/(self.cavity_radius**3 * (2*epsilon+1))
        Fons = dot(f*(linalg.inv(identity(3) \
                        - f*self.pol_1[0])),self.dipole)
        return Fons
    
    def get_Dipole_Approx(self,epsilon):
        """calculate dipole moments for a given epsilon by neglecting
        the structural distortions. Returned in Debye unit"""
        f = 2*(epsilon-1)/(self.cavity_radius**3 * (2*epsilon+1))
        b = linalg.inv(identity(3) - f*self.pol_1[0])
        inst_dipole = dot( (identity(3) + f*dot(self.pol_1[0],b)) , self.dipole) 
        inst_dipole*= self.BohrElectronToDebye 
        return inst_dipole
        
    def get_Esol_Approx(self,epsilon):
        """calculate solvation free energy for a given epsilon by neglecting
        the structural distortions. Returned in Hartree unit"""
        FI = self.get_FI(epsilon)
        inst_e_sol = -0.5 * tensordot(tensordot(FI,self.dipole,(0,0)),self.dipole,(0,0))
        return inst_e_sol
    
    def get_FI(self,epsilon):
        """calculate FI tensor for a given epsilon. Returned in AU unit"""
        f = 2*(epsilon-1)/(self.cavity_radius**3 * (2*epsilon+1))
        b = linalg.inv(identity(3) - f*self.pol_1[0])
        c = dot(transpose(b),dot(self.pol_1[0],b))
        FI = b*f +  c*f**2
        return FI

    def get_pol_fder(self):
        """calculate first derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        first_der_pol_cart = []
        for i in range(self.nAtoms*3):
            K = 1 + 4*i 
            ### FI(epsilon,a) tensor
            fd   = (1./12.) *  (\
                       (self.pol_1[K+3] - self.pol_1[K+0] ) \
                 + 8*  (self.pol_1[K+1] - self.pol_1[K+2] ))\
                 /(self.step_cart* self.AngstromToBohr)
                 
            first_der_pol_cart.append( fd ) 
        first_der_pol_cart = array(first_der_pol_cart)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_pol_mode = tensordot(self.L,first_der_pol_cart,(0,0)) 

        return first_der_pol_mode
    
    def get_pol_sder(self):
        """calculate second diagonal derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        ### FI(epsilon,a) tensor
        second_der_pol_mode = (1./12.) *  (\
            -     (self.pol_2[1+3] + self.pol_2[1+0] ) \
            +16*  (self.pol_2[1+1] + self.pol_2[1+2] ) \
            -30*   self.pol_2[0  ] )                \
            /self.step_mode**2

        return second_der_pol_mode
            
    def get_FI_fder(self,epsilon):
        """calculate first derivatives wrt normal coordinates of FI tensor
        for a given epsilon. Returned in AU unit"""
        f = 2*(epsilon-1)/(self.cavity_radius**3 * (2*epsilon+1))
        
        # cartesian first derivatives of FI tensor
        B_cart = zeros((1+4*self.nAtoms*3,3,3),dtype=float64)
        C_cart = zeros((1+4*self.nAtoms*3,3,3),dtype=float64)
            
        for i in range(1+4*self.nAtoms*3):
            b = linalg.inv(identity(3) - f*self.pol_1[i])
            B_cart[i] = b.copy()
            C_cart[i] = dot(transpose(b),dot(self.pol_1[i],b))        

        FI_cart = B_cart*f +  C_cart *f**2
        FI = FI_cart[0]
        
        first_der_FI_cart = []
        for i in range(self.nAtoms*3):
            K = 1 + 4*i 
            ### FI(epsilon,a) tensor
            fd   = (1./12.) *  (\
                       (FI_cart[K+3] - FI_cart[K+0] ) \
                 + 8*  (FI_cart[K+1] - FI_cart[K+2] ))\
                 /(self.step_cart* self.AngstromToBohr)
                 
            first_der_FI_cart.append( fd ) 
        first_der_FI_cart = array(first_der_FI_cart)
        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_FI_mode = tensordot(self.L,first_der_FI_cart,(0,0)) 
                
        return first_der_FI_mode
    
    def get_FI_sder(self,epsilon):
        """calculate second diagonal derivatives wrt normal coordinates of FI tensor
        for a given epsilon as well as normal coordinate specified in constructor 
        by mode_id key"""
        f = 2*(epsilon-1)/(self.cavity_radius**3 * (2*epsilon+1))
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
            /self.step_mode**2
             
        return second_der_FI_mode
    
    def __repr__(self):
        """print me!"""
        N = 150
        log = "\n"
        log+= " "+"-"*N+"\n"
        log+= "CHO-ONSAGER+POLARISATION FREQUENCY SHIFT".center(N)+"\n"
        log+= " "+"-"*N+"\n"
        log+= " %12s %12s %12s %12s %12s %12s %12s %12s %14s %14s %14s\n" \
                             % ('Ons. factor'.rjust(12),
                                'Shift'.rjust(12),
                                'Frequency'.rjust(12),
                                'Dip X'.rjust(12),
                                'Dip Y'.rjust(12),
                                'Dip Z'.rjust(12),
                                'AVER'.rjust(12),
                                'E_sol'.rjust(12),
                                'Field X'.rjust(12),
                                'Field Y'.rjust(12),
                                'Field Z'.rjust(12))
        log+= " %12s %12s %12s %12s %12s %12s %12s %12s %14s %14s %14s\n" % ('[---]'.rjust(12),
                                '[cm-1]'.rjust(12),
                                '[cm-1]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Debye]'.rjust(12),
                                '[Hartree]'.rjust(12),
                                '[AU]'.rjust(12),
                                '[AU]'.rjust(12),
                                '[AU]'.rjust(12),) 
        log+= " "+"-"*N+"\n"
        d = self.dipole * self.BohrElectronToDebye
        log+=     " %12.4f %12.2f %12.2f %12.4f %12.4f %12.4f %12.4f %12.6f %14.6f %14.6f %14.6f\n"\
            % (0, 0, self.freq_cmrec[self.mode_id],d[0], d[1], d[2], sqrt(sum(d**2)),0,0,0,0)
        for i,j in enumerate(self.frequency_shifts):
            log+= " %12.4f %12.2f %12.2f %12.4f %12.4f %12.4f %12.4f %12.6f %14.6f %14.6f %14.6f\n" % (self.onsager_factors[i], j,
                                         j+self.freq_cmrec[self.mode_id],
                                         self.dipole_moments[i,0],
                                         self.dipole_moments[i,1],
                                         self.dipole_moments[i,2],
                                         sqrt(sum(self.dipole_moments[i]**2)),
                                         self.solvation_energies[i],
                                         self.Onsager_fields[i,0],
                                         self.Onsager_fields[i,1],
                                         self.Onsager_fields[i,2])
        log+= " "+"-"*N+"\n"
        
        return str(log)
