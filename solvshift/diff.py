# ------------------------------------------------------------ #
#         D I F F E R E N T I A T I O N   M O D U L E          #
# ------------------------------------------------------------ #

from numpy      import array, float64, zeros, sqrt,\
                       sum, shape, ceil, transpose,\
                       tensordot
from units      import *
from dma        import *
from utilities  import *
from gaussfreq  import *
import os, glob

class DIFF(UNITS,FREQ):
    """contains usefull procedures for differentiation
       of DMA and molecular multipole moments in FFXpt-diag scheme
       where X=5,9"""
       
    def __init__(self,freq=0,step=0,
                 dir="",cartesian=False,L=0,
                 camm=False,pol=False,eds=False):

        self.camm = camm
        # number of normal modes
        self.L = L
        self.nModes = len(freq)
        self.nAtoms = (len(freq)+6)/3
        self.freq = freq
        # fragment names
        self.frag_names = self.ReadAtoms(dir)
        self.fEDS  = None
        self.sEDS  = None
        
        if not dir.endswith("/"): dir+="/"
        
        # step of differentiation and differentiation mode (pointity)
        self.step, self.n_point, self.file_type, self.mode_id, \
                   self.sderiv_wrk, self.sderiv_step = self.ReadStep(dir)
        if self.mode_id>0: self.calculate_sder=True
        else:              self.calculate_sder=False
        
        # DMA set from FFxpt input files
        if camm:
           self.DMA_set, smiec1, smiec2 = self.ParseDMA_set(dir,'coulomb')
           smiec1, self.Fragments_set, self.Overall_MM_set = self.ParseDMA_set(dir,self.file_type)
           for i,dma in enumerate(self.DMA_set):
               dma.pos = array(self.Fragments_set[i])
               dma.origin = array(self.Fragments_set[i])
           del smiec1, smiec2
        else:
           self.DMA_set, self.Fragments_set, self.Overall_MM_set = self.ParseDMA_set(dir,self.file_type)
           
        # number of distributed fragments
        self.nfrags = len(self.Fragments_set[0])
        
        # parse polarizability sets
        self.polarizability_set_1 = self.Parse_Polarizability(dir)
        self.polarizability_set_2 = self.Parse_Polarizability(self.sderiv_wrk)
        #if pol:
        self.fpol = 0#self.get_pol_fder()
        self.spol = 0#self.get_pol_sder()
        if eds:
           self.EDS_set_1 = self.Parse_EDS(dir+"w1c_nl/")
           self.EDS_set_2 = self.Parse_EDS(dir+"w1c_nl/sder/20/")
           #self.EDS_set_1 = self.Parse_EDS(dir+"el40/")
           #self.EDS_set_2 = self.Parse_EDS(dir+"el40/sder/20/")
           self.fEDS = self.get_EDS_fder()
           self.sEDS = self.get_EDS_sder()
        
        # --- [!] Calculate the derivatives!
        self.Fder, self.Sder ,self.FDip, self.SDip = self._CalcDerCartesian()
        #if cartesian: self.Fder, self.Sder ,self.FDip, self.SDip = self._CalcDerCartesian() <--- obsolete
        #else: self.Fder, self.Sder ,self.FDip, self.SDip = self._CalcDer() <--- obsolete
        if self.calculate_sder:
           if camm: 
              self.sder_DMA_set, smiec1, smiec2\
              =  self.ParseDMA_set(self.sderiv_wrk,'coulomb')
              smiec1, self.sder_DMA_frags, self.Overall_MM_set\
              =  self.ParseDMA_set(self.sderiv_wrk,self.file_type)
              del smiec1, smiec2
           else:
              self.sder_DMA_set, self.sder_DMA_frags, self.Overall_MM_set\
              =  self.ParseDMA_set(self.sderiv_wrk,self.file_type)
           
           for i,dma in enumerate(self.sder_DMA_set):
               dma.pos = array(self.sder_DMA_frags[i])
               dma.origin = array(self.sder_DMA_frags[i])
               
           self.Sder[self.mode_id],self.sder_MM = \
                self._Sderivatives(self.sderiv_wrk,self.mode_id,self.sderiv_step) 
                
        # --- [!] Calculate IR Harmonic intensities
        self.IR_Harm_Int = self.CalcIrInt()
        
        # --- [!] Calculate CAMM derivatives!
        if 0:#camm:
           self.fder_CAMM_set = self.Fderivatives_CAMM(dir=dir,
                                                       basis='',
                                                       step=self.step)
           if self.calculate_sder:
              self.sder_CAMM_set = self.Sderivatives_CAMM(dir=self.sderiv_wrk,
                                                          basis='',
                                                          step=self.sderiv_step)


    def get_EDS_fder(self):
        """calculate first derivatives wrt normal coordinates
        of interaction energy components from EDS"""
        
        first_der_EDS_cart = []
        for i in range(self.nAtoms*3):
            K = 1 + 4*i
            fd   = (1./12.) *  (\
                       (self.EDS_set_1[K+3] - self.EDS_set_1[K+0] ) \
                 + 8*  (self.EDS_set_1[K+1] - self.EDS_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            first_der_EDS_cart.append( fd )
        first_der_EDS_cart = array(first_der_EDS_cart)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_EDS_mode = tensordot(self.L,first_der_EDS_cart,(0,0))

        return first_der_EDS_mode

    def get_EDS_sder(self):
        """calculate second diagonal derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        second_der_EDS_mode = (1./12.) *  (\
            -     (self.EDS_set_2[1+3] + self.EDS_set_2[1+0] ) \
            +16*  (self.EDS_set_2[1+1] + self.EDS_set_2[1+2] ) \
            -30*   self.EDS_set_2[0  ] )                \
            /self.sderiv_step**2

        return second_der_EDS_mode
                
    def get_pol_fder(self):
        """calculate first derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        first_der_pol_cart = []
        for i in range(self.nAtoms*3):
            K = 1 + 4*i
            fd   = (1./12.) *  (\
                       (self.polarizability_set_1[K+3] - self.polarizability_set_1[K+0] ) \
                 + 8*  (self.polarizability_set_1[K+1] - self.polarizability_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            first_der_pol_cart.append( fd )
        first_der_pol_cart = array(first_der_pol_cart)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_pol_mode = tensordot(self.L,first_der_pol_cart,(0,0))

        return first_der_pol_mode
    
    def get_pol_sder(self):
        """calculate second diagonal derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        second_der_pol_mode = (1./12.) *  (\
            -     (self.polarizability_set_2[1+3] + self.polarizability_set_2[1+0] ) \
            +16*  (self.polarizability_set_2[1+1] + self.polarizability_set_2[1+2] ) \
            -30*   self.polarizability_set_2[0  ] )                \
            /self.sderiv_step**2

        return second_der_pol_mode

    def ParseDmatFromFchk(self,file,basis_size):
        """parses density matrix from Gaussian fchk file"""
        
        data = open(file)
        line = data.readline()
        querry = "Total SCF Density"
        while 1:
            if querry in line: break
            line = data.readline()
        N = int(line.split()[-1])
        line = data.readline()
        dmat = []
        for i in range(int(ceil(N/5.))): 
            dmat+=[x for x in line.split()] 
            line = data.readline()
            
        # construct explicit 2D density matrix
        P = zeros((basis_size,basis_size),dtype=float64)
        #I = 0
        for i in range(basis_size):
            for j in range(i+1):
                P[i,j] = float64(dmat.pop(0))#dmat[I]
                P[j,i] = P[i,j] #dmat[I]
                #I += 1
       
        return array(P)
 

    if os.environ['__IMPORT__COULOMB__']:
        
       def Fderivatives_CAMM(self,dir,basis,step):
           """calculates CAMM_set from fchk files of FF directory"""
           #files2= os.listdir(dir)
           #files = [ ]
           #for i in files2:
           #    if i.endswith('_.fchk'):
           #       files.append(i)
           #del files2
           ## sort the input files (necessary!)
           #files.sort()  
           ##print files
           ##print len(files)
       
           #CAMM_set = []
           #basis_size = len(Ints.getbasis(molecule,'6-311++G**'))
           ## calculate CAMMs from respective density matrices
           #for i,file in enumerate(files):
           #    dmat = self.ParseDmatFromFchk(file,basis_size)
           #    ### create Molecule object
           #    structure = []
           #    for j in range(self.nfrags):
           #        structure.append( (self.atomic_numbers[self.frag_names[j]],
           #                           self.Fragments_set[i,j]) ) 
           #    molecule = Molecule('mol',
           #                         structure,
           #                         multiplicity=1,
           #                         charge=0,
           #                         units='Bohr')
           #    ### calculate CAMMs                    
           #    CAMM = multip.MULTIP(molecule=molecule,
           #                         basis='6-311++G**',
           #                        #basis='sto-3g',
           #                         method='b3lyp',
           #                         matrix=dmat,
           #                         transition=False)
           #    CAMM.camms()
           #    del dmat
           #    ### construct DMA object
           #    camm = DMA(nfrag=self.nfrags)
           #    camm.DMA[0] = CAMM.Mon
           #    camm.DMA[1] = CAMM.Dip # DODAJ resztę!!!!
           #
           #   CAMM_set.append( camm )
           
           # calculate first derivatives!
           s1 = 1./(step * self.AngstromToBohr)
           for i in range(self.nAtoms*3):
               K = 1 + 4*i 
               fder_CAMM = (1./12.) *  (\
                          (CAMM_set[K+3] - CAMM_set[K+0] ) \
                    + 8*  (CAMM_set[K+1] - CAMM_set[K+2] ))\
                    *s1
                
           ### transform derivatives from cart to normal mode space         
           fder_CAMM = DMAMadrixMultiply(transpose(self.L),fder_CAMM) 
          
           return fder_CAMM
       
           
       def Sderivatives_CAMM(self,dir,basis_size,step):
           """calculates CAMM_set from fchk files of FF directory"""
           files2= os.listdir(dir)
           files = [ ]
           for i in files2:
               if i.endswith('_.fchk'):
                  files.append(i)
           del files2
           # sort the input files (necessary!)
           files.sort()  
           #print files
           #print len(files)
       
           CAMM_set = []
           basis_size = len(Ints.getbasis(molecule,'6-311++G**'))
           # calculate CAMMs from respective density matrices
           for i,file in enumerate(files):
               dmat = self.ParseDmatFromFchk(file,basis_size)
               ### create Molecule object
               structure = []
               for j in range(self.nfrags):
                   structure.append( (self.atomic_numbers[self.frag_names[j]],
                                      self.Fragments_set[i,j]) ) 
               molecule = Molecule('mol',
                                    structure,
                                    multiplicity=1,
                                    charge=0,
                                    units='Bohr')
               ### calculate CAMMs                    
               CAMM = multip.MULTIP(molecule=molecule,
                                    basis='6-311++G**',
                                   #basis='sto-3g',
                                    method='b3lyp',
                                    matrix=dmat,
                                    transition=False)
               CAMM.camms()
               del dmat
               ### construct DMA object
               camm = DMA(nfrag=self.nfrags)
               camm.DMA[0] = CAMM.Mon
               camm.DMA[1] = CAMM.Dip # DODAJ resztę!!!!
           
               CAMM_set.append( camm )
           
           # calculate second derivatives!
           s2 = 1./step**2
           for i in range(1):
               K = 1 + 4*i 
               sder_DMA = (1./12.) *  (\
                    -     (CAMM_set[K+3] + CAMM_set[K+0] ) \
                    +16*  (CAMM_set[K+1] + CAMM_set[K+2] ) \
                    -30*   CAMM_set[0  ] )                 \
                    *s2
           return sder_CAMM
       
    def _CalcDerCartesian(self):
        """calculates FIRST (and diagonal second - but ONLY diagonal!!!)
           derivatives wrt cartesian coordinates of:
           - DMA distribution
           - MMM distribution
           """
           
        # calculate set of molecular dipole moments in [Bohr*ElectronCharge]
        #Dipole_Moments = []
        #for i,dma in enumerate(self.DMA_set):
        #    moment = dma[0].reshape((self.nfrags,1)) * self.Fragments_set[i] + dma[1]
        #    Dipole_Moments.append(sum(moment,axis=0))
        #Dipole_Moments = array(Dipole_Moments)
        #print "Dipole moments [au] \n"
        #print Dipole_Moments
        # store first derivatives in the lists
        Fder = []
        Sder = []
        Fmmm = []
        Smmm = []
        # change origin in the case of CAMM
        if self.camm:
           reference_structure = array(self.Fragments_set[0])
           for i,dma in enumerate(self.DMA_set):
               #print i,"fder"
               dma.MAKE_FULL()
               dma.ChangeOrigin(new_origin_set=reference_structure)
        #print Dipole_Moments
        s1 = 1./self.step
        s2 = 1./(self.step**2)
        ### 5-point formulae
        if self.n_point == 5:
          for i in range(self.nAtoms*3):
            #print "calculation for %d mode"%(i+1)
            K = 1 + 4*i 
            
            # first derivatves of DMA
            fder_DMA = (1./12.) *  (\
                       (self.DMA_set[K+3] - self.DMA_set[K+0] ) \
                 + 8*  (self.DMA_set[K+1] - self.DMA_set[K+2] ))\
                 * s1

            # second derivatves of DMA
            sder_DMA = (1./12.) *  (\
                 -     (self.DMA_set[K+3] + self.DMA_set[K+0] ) \
                 +16*  (self.DMA_set[K+1] + self.DMA_set[K+2] ) \
                 -30*   self.DMA_set[0  ] )                     \
                 * s2

            # first derivatves of molecular multipole moments     
            fder_MMM = (1./12.) *  (\
                       (self.Overall_MM_set[K+3] - self.Overall_MM_set[K+0] ) \
                 + 8*  (self.Overall_MM_set[K+1] - self.Overall_MM_set[K+2] ))\
                 * s1
                 
            # second derivatves of molecular multipole moments   
            sder_MMM = (1./12.) *  (\
                 -     (self.Overall_MM_set[K+3] + self.Overall_MM_set[K+0] ) \
                 +16*  (self.Overall_MM_set[K+1] + self.Overall_MM_set[K+2] ) \
                 -30*   self.Overall_MM_set[0  ] )                     \
                 * s2
            
                             
            # first derivatves of molecular dipole moment     
            #fdip     = (1./12.) *  (\
            #           (Dipole_Moments[K+3] - Dipole_Moments[K+0] ) \
            #    + 8*  (Dipole_Moments[K+1] - Dipole_Moments[K+2] ))\
            #     * s1
                 
            # second derivatves of molecular dipole moment     
            #sdip     = (1./12.) *  (\
            #     -     (Dipole_Moments[K+3] + Dipole_Moments[K+0] ) \
            #     +16*  (Dipole_Moments[K+1] + Dipole_Moments[K+2] ) \
            #     -30*   Dipole_Moments[0  ] )                     \
            #     * s2

            Fder.append( fder_DMA/self.AngstromToBohr )
            Sder.append( sder_DMA/self.AngstromToBohr**2 )
            Fmmm.append( fder_MMM/self.AngstromToBohr )
            Smmm.append( sder_MMM/self.AngstromToBohr**2 )

        ### 9-point formulae
        if self.n_point == 9:
            
          for i in range(self.nAtoms*3):
            K = 1 + 8*i 

            # first derivatves of DMA
            fder_DMA = (\
                (4./5.)   * (self.DMA_set[K+3] - self.DMA_set[K+4] ) \
               -(1./5.)   * (self.DMA_set[K+2] - self.DMA_set[K+5] ) \
               +(4./105.) * (self.DMA_set[K+1] - self.DMA_set[K+6] ) \
               -(1./280.) * (self.DMA_set[K+0] - self.DMA_set[K+7] ))\
               * s1

            # second derivatves of DMA
            sder_DMA = (\
                (8./5.)   * (self.DMA_set[K+3] + self.DMA_set[K+4] ) \
               -(1./5.)   * (self.DMA_set[K+2] + self.DMA_set[K+5] ) \
               +(8./315.) * (self.DMA_set[K+1] + self.DMA_set[K+6] ) \
               -(1./560.) * (self.DMA_set[K+0] + self.DMA_set[K+7] ) \
               -(205./72.)*  self.DMA_set[0  ] )\
               * s2
                 
            # first derivatves of molecular dipole moment
            fdip = (\
                (4./5.)   * (Dipole_Moments[K+3] - Dipole_Moments[K+4] ) \
               -(1./5.)   * (Dipole_Moments[K+2] - Dipole_Moments[K+5] ) \
               +(4./105.) * (Dipole_Moments[K+1] - Dipole_Moments[K+6] ) \
               -(1./280.) * (Dipole_Moments[K+0] - Dipole_Moments[K+7] ))\
               * s1
                 
            # second derivatves of molecular dipole moment
            sdip = (\
                (8./5.)   * (Dipole_Moments[K+3] + Dipole_Moments[K+4] ) \
               -(1./5.)   * (Dipole_Moments[K+2] + Dipole_Moments[K+5] ) \
               +(8./315.) * (Dipole_Moments[K+1] + Dipole_Moments[K+6] ) \
               -(1./560.) * (Dipole_Moments[K+0] + Dipole_Moments[K+7] ) \
               -(205./72.)*  Dipole_Moments[0  ] )\
               * s2

            Fder.append( fder_DMA/self.AngstromToBohr )
            Sder.append( sder_DMA/self.AngstromToBohr**2 )
            FDip.append( fdip/self.AngstromToBohr )
            SDip.append( sdip/self.AngstromToBohr**2 )



        # transform to normal mode space  
        Fder_DMA = DMAMadrixMultiply(transpose(self.L),Fder)
        Sder_DMA = [ DMA(nfrag=self.nfrags) for x in range(self.nModes) ]  
        #FDip = dot( transpose(self.L) ,array( FDip ) )
        Fder_MMM = DMAMadrixMultiply(transpose(self.L),Fmmm)
        Sder_MMM = zeros((self.nModes,3))
        
        return Fder_DMA, Sder_DMA, Fder_MMM, Sder_MMM

    def _Sderivatives(self,sderiv_wrk,mode_id,step):
        """calculates second derivatives wrt NORMAL (!) coordinates
        for selected mode. It requires input files for normal mode displacements!!!
        """
        # store first and diagonal second derivatives in the lists
        s1 = 1./step
        s2 = 1./(step**2)

        # change origin in the case of CAMM
        if self.camm:
           reference_structure = array(self.Fragments_set[0])
           for i,dma in enumerate(self.sder_DMA_set):
               #print i, "sder"
               dma.MAKE_FULL()
               dma.ChangeOrigin(new_origin_set=reference_structure)
               
        ### 5-point formulae
        if self.n_point == 5:
          for i in range(1):
            K = 1 + 4*i 
            # first derivatves of DMA
            fder_DMA = (1./12.) *  (\
                       (self.sder_DMA_set[K+3] - self.sder_DMA_set[K+0] ) \
                 + 8*  (self.sder_DMA_set[K+1] - self.sder_DMA_set[K+2] ))\
                 *s1
            # second derivatives of DMA
            sder_DMA = (1./12.) *  (\
                 -     (self.sder_DMA_set[K+3] + self.sder_DMA_set[K+0] ) \
                 +16*  (self.sder_DMA_set[K+1] + self.sder_DMA_set[K+2] ) \
                 -30*   self.sder_DMA_set[0  ] )                     \
                 *s2

            if self.file_type.lower() == "gaussian":
             # first derivatves of MM
             fder_MM  = (1./12.) *  (\
                        (self.Overall_MM_set[K+3] - self.Overall_MM_set[K+0] ) \
                  + 8*  (self.Overall_MM_set[K+1] - self.Overall_MM_set[K+2] ))\
                  *s1
             # second derivatives of MM
             sder_MM  = (1./12.) *  (\
                  -     (self.Overall_MM_set[K+3] + self.Overall_MM_set[K+0] ) \
                  +16*  (self.Overall_MM_set[K+1] + self.Overall_MM_set[K+2] ) \
                  -30*   self.Overall_MM_set[0  ] )                     \
                  *s2

        #print "INPUT FOR FDER"
        #print fder_DMA
        #print "INPUT FOR SDER"
        #print sder_DMA
        #print "INPUT FOR FDER"
        #print fder_MM
        #print "INPUT FOR SDER"
        #print sder_MM.DMA[1][-1]
        #print sder_MM
        #print "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        #for i in self.Overall_MM_set:
        #    print i.DMA[1]
        return sder_DMA, sder_MM
            
            
    def ParseDMA_set(self,dir,type):
        """withdraws DMA_set from log files of FF directory"""
        files2= os.listdir(dir)
        files = [ ]
        if type == 'coulomb':        
           files= glob.glob('%s*_.camm' % dir)
        else:
           files= glob.glob('%s*_.log' % dir)
        del files2
        # sort the input files (necessary!)
        files.sort()  
        #print files
        #print len(files)
        
        # construct DMA_set, Fragments_set and Overall_MM_set
        Fragments_set = []
        DMA_set = []
        Overall_MM_set = []
        for file in files:
            dma, Fragment = ParseDMA ( file, type ) # zmieniono z ( dir+file, type )
            DMA_set.append( dma )
            Fragments_set.append( Fragment )
            
            if type.lower() == 'gaussian':
               ### --- TIP ---
               ### w przyszłości dopisz w klasie FREQ argument w konstruktorze 'type'
               ### aby umożliwić sczytywanie momentów z plików GAMESS lub innych.
               ### Wtedy usuń warunek if w tej funkcji!
               # construct overall multipoles set DMA object
               dipole     = self.Dipole(file)# zmieniono z ( dir+file)
               quadrupole = self.Quadrupole(file)# zmieniono z ( dir+file )
               octupole   = self.Octupole(file)# zmieniono z ( dir+file )
            
               MM = DMA(nfrag=1)
               MM.pos = zeros(3,dtype=float64)
               MM.DMA[1][0] = array(dipole)
               MM.DMA[2][0] = array(quadrupole)
               MM.DMA[3][0] = array(octupole)

               Overall_MM_set.append( MM )
        
        return DMA_set, array(Fragments_set), Overall_MM_set

    def Parse_EDS(self,dir):
        """parses files to collect polarisabilities"""
        
        files2= os.listdir(dir)
        files = [ ]
        for i in files2:
            if i.endswith('_.log'):
               files.append(i)
        del files2
        # sort the input files (necessary!)
        files.sort()  

        set = []
        for file in files:
            set.append(Parse_EDS_InteractionEnergies(dir+file))
            
        return array(set)

    def Parse_Polarizability(self,dir):
        """parses files to collect polarisabilities"""
        
        files2= os.listdir(dir)
        files = [ ]
        for i in files2:
            if i.endswith('_.log'):
               files.append(i)
        del files2
        # sort the input files (necessary!)
        files.sort()  

        set = []
        for file in files:
            set.append(self.Polarizability(dir+file))
            
        return array(set)
    
    def CalcIrInt(self):
        """calculate IR Harmonic intensities in [kM/mole]"""
        FDip = zeros((self.nModes,3),dtype=float64)
        #print self.FDip[7]
        for mode in range(self.nModes): FDip[mode] = self.FDip[mode].DMA[1]
        return sum(FDip**2,axis=1) * self.BohrElectronToDebye**2\
                /(self.BohrToAngstrom**2 * self.ElectronMassToAmu) * self.IrIntToKmM

    
    def ReadStep(self,dir):
        """reads the differentiation step, pointity and package
           from a setup file slv.step"""
        
        step_file = open('%s/slv.step' % dir,'r')
        # step
        step = float64(step_file.readline().split()[-1])
        # pointity
        n_point = int(step_file.readline().split()[-1])
        # file type
        file_type = step_file.readline().split()[-1]
        # mode id for sderiv calculation (can be negative then no sderiv calc)
        mode_id  = int(step_file.readline().split()[-1])
        # sderiv working directory
        sderiv_wrk= step_file.readline().split()[-1]
        # step of differentiation wrt normal modes for sderv routine
        sderiv_step = float64(step_file.readline().split()[-1])

        return step , n_point, file_type, mode_id, sderiv_wrk, sderiv_step
    
    def ReadAtoms(self,dir):
        """reads the differentiation step, pointity and package
           from a setup file slv.step"""
        
        # atomic/fragment symbols 
        frag_file = open('%s/slv.frags' % dir,'r')
        frag_names = []
        line = frag_file.readline()
        while 1:
           if not line: break
           frag_names.append( line.split()[-1])
           line = frag_file.readline()
        return frag_names
        
    def __repr__(self):
        """prints:
        Harmonic IR intensities [kM/mole] 
        First and Second derivatives of DMA
        """
        
        # Harmonic IR intensities
        log = "\n"
        log+= " HARMONIC IR INTENSITIES [kM/mole]\n"
        log+= "\n"
        log+= "%10s %10s %10s\n" %\
           ("Freq".rjust(10),"[cm-1]".rjust(10),
            "Int".rjust(10))
        for mode, intensity in enumerate(self.IR_Harm_Int):
            log+= "%10d %10.2f %10.4f\n" %\
            (mode+1,self.freq[mode],intensity)
        log+= "\n"
        # Derivatives of DMA distribution
        log+= " ---  F I R S T    D E R I V A T I V E S  ---\n"
        for mode, derivatives in enumerate(self.Fder):
            log+= " ================ MODE %d:%8.2f [cm-1] ================ \n" %\
                           (mode+1,self.freq[mode])
            log+= "\n"
            log+= repr(derivatives)
        log+= " ---  S E C O N D     D E R I V A T I V E S  ---\n"
        for mode, derivatives in enumerate(self.Sder):
            log+= " ================ MODE %d:%8.2f [cm-1] ================ \n" %\
                           (mode+1,self.freq[mode])
            log+= "\n"
            log+= repr(derivatives)
            
        return str(log)

    def _CalcDer_old(self):
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
        #print "Dipole moments [au] \n"
        #print Dipole_Moments
        # store first and diagonal second derivatives in the lists
        Fder = [] 
        Sder = []
        FDip = []
        SDip = []
        #print Dipole_Moments
        s1 = 1./self.step
        s2 = 1./(self.step**2)
        ### 5-point formulae
        if self.n_point == 5:
          for i in range(self.nModes):
            #print "calculation for %d mode"%(i+1)
            K = 1 + 4*i 
            
            # first derivatves of DMA
            fder_DMA = (1./12.) *  (\
                       (self.DMA_set[K+3] - self.DMA_set[K+0] ) \
                 + 8*  (self.DMA_set[K+1] - self.DMA_set[K+2] ))\
                 * s1

            # second derivatves of DMA
            sder_DMA = (1./12.) *  (\
                 -     (self.DMA_set[K+3] + self.DMA_set[K+0] ) \
                 +16*  (self.DMA_set[K+1] + self.DMA_set[K+2] ) \
                 -30*   self.DMA_set[0  ] )                     \
                 * s2
                 
            # first derivatves of molecular dipole moment     
            fdip     = (1./12.) *  (\
                       (Dipole_Moments[K+3] - Dipole_Moments[K+0] ) \
                 + 8*  (Dipole_Moments[K+1] - Dipole_Moments[K+2] ))\
                 * s1
                 
            # second derivatves of molecular dipole moment     
            sdip     = (1./12.) *  (\
                 -     (Dipole_Moments[K+3] + Dipole_Moments[K+0] ) \
                 +16*  (Dipole_Moments[K+1] + Dipole_Moments[K+2] ) \
                 -30*   Dipole_Moments[0  ] )                     \
                 * s2

            Fder.append( fder_DMA )
            Sder.append( sder_DMA )
            FDip.append( fdip )
            SDip.append( sdip )

        ### 9-point formulae
        if self.n_point == 9:
          for i in range(self.nModes):
            K = 1 + 8*i 

            # first derivatves of DMA
            fder_DMA = (\
                (4./5.)   * (self.DMA_set[K+3] - self.DMA_set[K+4] ) \
               -(1./5.)   * (self.DMA_set[K+2] - self.DMA_set[K+5] ) \
               +(4./105.) * (self.DMA_set[K+1] - self.DMA_set[K+6] ) \
               -(1./280.) * (self.DMA_set[K+0] - self.DMA_set[K+7] ))\
               * s1

            # second derivatves of DMA
            sder_DMA = (\
                (8./5.)   * (self.DMA_set[K+3] + self.DMA_set[K+4] ) \
               -(1./5.)   * (self.DMA_set[K+2] + self.DMA_set[K+5] ) \
               +(8./315.) * (self.DMA_set[K+1] + self.DMA_set[K+6] ) \
               -(1./560.) * (self.DMA_set[K+0] + self.DMA_set[K+7] ) \
               -(205./72.)*  self.DMA_set[0  ] )\
               * s2
                 
            # first derivatves of molecular dipole moment     
            fdip = (\
                (4./5.)   * (Dipole_Moments[K+3] - Dipole_Moments[K+4] ) \
               -(1./5.)   * (Dipole_Moments[K+2] - Dipole_Moments[K+5] ) \
               +(4./105.) * (Dipole_Moments[K+1] - Dipole_Moments[K+6] ) \
               -(1./280.) * (Dipole_Moments[K+0] - Dipole_Moments[K+7] ))\
               * s1
                 
            # second derivatves of molecular dipole moment     
            sdip = (\
                (8./5.)   * (Dipole_Moments[K+3] + Dipole_Moments[K+4] ) \
               -(1./5.)   * (Dipole_Moments[K+2] + Dipole_Moments[K+5] ) \
               +(8./315.) * (Dipole_Moments[K+1] + Dipole_Moments[K+6] ) \
               -(1./560.) * (Dipole_Moments[K+0] + Dipole_Moments[K+7] ) \
               -(205./72.)*  Dipole_Moments[0  ] )\
               * s2
            #sdip     = 1 *  (\
            #            (Dipole_Moments[K+1] + Dipole_Moments[K+2] - 
            #             2* Dipole_Moments[0]) * s2 )          
            #     
            Fder.append( fder_DMA )
            Sder.append( sder_DMA )
            FDip.append( fdip )
            SDip.append( sdip )
                        
        return Fder, Sder, array( FDip ), array( SDip )
    

if __name__ == '__main__':
   def ParseDmatFromFchk(file,basis_size):
        """parses density matrix from Gaussian fchk file"""
        
        data = open(file)
        line = data.readline()
        querry = "Total SCF Density"
        while 1:
            if querry in line: break
            line = data.readline()
        N = int(line.split()[-1])

        line = data.readline()
        dmat = []
        for i in range(int(ceil(N/5.))): 
            dmat+=[x for x in line.split()] 
            line = data.readline()
        print len(dmat)
        #dmat = array(dmat,dtype=float64)
        
        # construct explicit 2D density matrix
        P = zeros((basis_size,basis_size),dtype=float64)
        #I = 0
        for i in range(basis_size):
            for j in range(i+1):
                P[i,j] = float64(dmat.pop(0))#dmat[I]
                P[j,i] = P[i,j] #dmat[I]
                #I += 1

        return array(P)    
    
   def test1(): 
       """calculates CAMMs from density matrux from GAUSSIAN09
       using COULOMB.py routines"""
       
       from sys import argv
       
       file_fchk = argv[1]
       file_log  = argv[2]
       dma, fragment = ParseDMA( file_log, 'gaussian' )
       
       frag_file = open('slv.frags','r')
       frag_names = []
       line = frag_file.readline()
       while 1:
          if not line: break
          frag_names.append( line.split()[-1])
          line = frag_file.readline()  

       ### create Molecule object
       structure = []
       for j in range(len(fragment)):
           structure.append( (UNITS.atomic_numbers[frag_names[j]],
                              fragment[j]) ) 
       molecule = Molecule('mol',
                            structure,
                            multiplicity=1,
                            charge=0,
                            units='Bohr')
                            
       basis_size = len(Ints.getbasis(molecule,'6-311++G**'))
       dmat = ParseDmatFromFchk(file_fchk,basis_size)
       
       ### calculate CAMMs                    
       CAMM = multip.MULTIP(molecule=molecule,
                     basis='6-311++G**',
                     #basis='sto-3g',
                     method='b3lyp',
                     matrix=dmat,
                     transition=False)
       CAMM.camms()
       #CAMM.__printCAMMs__()
       
       result = DMA(nfrag=12)
       result.DMA[0][:] = CAMM.Mon
       #
       result.DMA[1][:] = CAMM.Dip
       #
       result.DMA[2][:,0] = array(CAMM.Quad)[:,0,0]
       result.DMA[2][:,1] = array(CAMM.Quad)[:,1,1]
       result.DMA[2][:,2] = array(CAMM.Quad)[:,2,2]
       result.DMA[2][:,3] = array(CAMM.Quad)[:,0,1]
       result.DMA[2][:,4] = array(CAMM.Quad)[:,0,2]
       result.DMA[2][:,5] = array(CAMM.Quad)[:,1,2]
       #
       result.DMA[3][:,0] = array(CAMM.Oct)[:,0,0,0]
       result.DMA[3][:,1] = array(CAMM.Oct)[:,1,1,1]
       result.DMA[3][:,2] = array(CAMM.Oct)[:,2,2,2]
       result.DMA[3][:,3] = array(CAMM.Oct)[:,0,0,1]
       result.DMA[3][:,4] = array(CAMM.Oct)[:,0,0,2]
       result.DMA[3][:,5] = array(CAMM.Oct)[:,0,1,1]
       result.DMA[3][:,6] = array(CAMM.Oct)[:,1,1,2]
       result.DMA[3][:,7] = array(CAMM.Oct)[:,0,2,2]
       result.DMA[3][:,8] = array(CAMM.Oct)[:,1,2,2]
       result.DMA[3][:,9] = array(CAMM.Oct)[:,0,1,2]
       #
       print result
       out = open('test.dat','w')
       out.write(str(result))
       out.close()
   
   ### first test
   #test1()
   def CalculateCAMM(basis='6-311++G**'): 
       """calculates CAMMs from density matrix from GAUSSIAN09
       using COULOMB.py routines"""
       from head import *
       from sys import argv
       import os, glob
       
       pliki_fchk  = glob.glob('./*_.fchk')
       pliki_fchk.sort()
       pliki_log   = glob.glob('./*_.log')
       pliki_log .sort()    
       print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
       for i in range(len(pliki_log)):
           print pliki_log[i], pliki_fchk[i]
       print
       
       for i,file_log in enumerate(pliki_log):
           dma, fragment = ParseDMA( file_log, 'gaussian' )
       
           frag_file = open('slv.frags','r')
           frag_names = []
           line = frag_file.readline()
           while 1:
              if not line: break
              frag_names.append( line.split()[-1])
              line = frag_file.readline()  

           ### create Molecule object
           structure = []
           for j in range(len(fragment)):
               #structure.append( (UNITS.atomic_numbers[frag_names[j]],
               structure.append( (Atom(frag_names[j]).atno,
                                  fragment[j]) ) 
           molecule = Molecule('mol',
                                structure,
                                multiplicity=1,
                                charge=0,
                                units='Bohr')
                            
           basis_size = len(Ints.getbasis(molecule,basis))
           print " - basis size= ",basis_size
           dmat = ParseDmatFromFchk(pliki_fchk[i],basis_size)
       
           ### calculate CAMMs                    
           CAMM = multip.MULTIP(molecule=molecule,
                         basis=basis,
                         #basis='sto-3g',
                         method='b3lyp',
                         matrix=dmat,
                         transition=False)
           CAMM.camms()
           CAMM.mmms()
           CAMM.__printMMMs__()
           #CAMM.__printCAMMs__()
       
           result = DMA(nfrag=len(structure))
           result.DMA[0][:] = CAMM.Mon
           #
           result.DMA[1][:] = CAMM.Dip
           #
           result.DMA[2][:,0] = array(CAMM.Quad)[:,0,0]
           result.DMA[2][:,1] = array(CAMM.Quad)[:,1,1]
           result.DMA[2][:,2] = array(CAMM.Quad)[:,2,2]
           result.DMA[2][:,3] = array(CAMM.Quad)[:,0,1]
           result.DMA[2][:,4] = array(CAMM.Quad)[:,0,2]
           result.DMA[2][:,5] = array(CAMM.Quad)[:,1,2]
           #
           result.DMA[3][:,0] = array(CAMM.Oct)[:,0,0,0]
           result.DMA[3][:,1] = array(CAMM.Oct)[:,1,1,1]
           result.DMA[3][:,2] = array(CAMM.Oct)[:,2,2,2]
           result.DMA[3][:,3] = array(CAMM.Oct)[:,0,0,1]
           result.DMA[3][:,4] = array(CAMM.Oct)[:,0,0,2]
           result.DMA[3][:,5] = array(CAMM.Oct)[:,0,1,1]
           result.DMA[3][:,6] = array(CAMM.Oct)[:,1,1,2]
           result.DMA[3][:,7] = array(CAMM.Oct)[:,0,2,2]
           result.DMA[3][:,8] = array(CAMM.Oct)[:,1,2,2]
           result.DMA[3][:,9] = array(CAMM.Oct)[:,0,1,2]
           #
           #print result
           out = open(file_log[:-4]+'.camm','w')
           out.write(str(result))
           out.close()
           print " Writing file:  :", file_log[:-4]+'.camm'
       print
           
   ### calculate camms!
   from sys  import argv
   CalculateCAMM(basis=argv[1])

   
