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
import os, glob, PyQuante.Ints, coulomb.multip
os.environ['__IMPORT__COULOMB__']='1'

__all__ = ['DIFF',]
__version__ = '5.2.1'
                
class DIFF(UNITS,FREQ):
    """contains usefull procedures for differentiation
of DMA and molecular multipole moments in FFXpt-diag scheme
where X=5,9"""
       
    def __init__(self,freq=0,step=0,
                 dir="",cartesian=False,L=0,
                 camm=False,pol=False,eds=False,efp=False,
                 solefp=False,nae=0,basis='6-311++G**',
                 dpol=False,lprint=True):

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
        self.__lprint = lprint
        
        if not dir.endswith("/"): dir+="/"
        
        # step of differentiation and differentiation mode (pointity)
        self.step, self.n_point, self.file_type, self.mode_id, \
        self.sderiv_wrk, self.sderiv_step = self.ReadStep(dir)
        if self.mode_id>0: self.calculate_sder=True
        else:              self.calculate_sder=False
        
        if (not eds  and not efp and not solefp and not dpol):
          # DMA set from FFxpt input files
          if camm:
             self.DMA_set, smiec1 = self.ParseDMA_set(dir,'coulomb')
             #smiec1, self.Overall_MM_set = self.ParseDMA_set(dir,self.file_type)
             del smiec1
          else:
             self.DMA_set, self.Overall_MM_set = self.ParseDMA_set(dir,self.file_type)
        
          self.reference_str   = self.DMA_set[0].get_pos()
          self.reference_origin= self.DMA_set[0].get_origin()
          # number of distributed fragments
          ###self.nfrags = len(self.reference_str)
          self.nfrags = len(self.reference_origin)
        
          if pol:
             self.polarizability_set_1 = self.Parse_Polarizability(dir)
             self.polarizability_set_2 = self.Parse_Polarizability(self.sderiv_wrk)
             self.fpol = self.get_pol_fder()
             self.spol = 0#self.get_pol_sder()
        
          # --- [!] Calculate the derivatives!
          #if self.calculate_sder:
          #   self.Fder, self.Sder ,self.FDip, self.SDip = self._CalcDerCartesian()
          #else:
          if self.calculate_sder:
              self.sder_DMA_set, smiec1\
                =  self.ParseDMA_set(self.sderiv_wrk,'coulomb')
          self.Fder, self.Sder = self._NEW(sder=self.calculate_sder)
        
          if 0:#self.calculate_sder:
             if camm:
                self.sder_DMA_set, smiec1\
                =  self.ParseDMA_set(self.sderiv_wrk,'coulomb')
                #smiec1, self.Overall_MM_set\
                #=  self.ParseDMA_set(self.sderiv_wrk,self.file_type)
                del smiec1
             else:
                self.sder_DMA_set, self.Overall_MM_set\
                =  self.ParseDMA_set(self.sderiv_wrk,self.file_type)
               
             self.Sder[self.mode_id],self.sder_MM = \
                  self._Sderivatives(self.sderiv_wrk,self.mode_id,self.sderiv_step) 
                
          # --- [!] Calculate IR Harmonic intensities
          #self.IR_Harm_Int = self.CalcIrInt()
        
        # solvatochromic EFP
        elif solefp:
           fdir = solefp
           print " PIPEK-MEZEY LOCALIZATION AND LMTP EVALUATION..."
           self.tran_set_1, self.lmoc_set_1, self.vecl_set_1 = self.Parse_tran(dir+fdir,nae,basis)
           self.U = self.vecl_set_1[0]
           self.rLMO = self.lmoc_set_1[0]
           print " PARSING FOCK MATRICES..."
           self.fock_set_1 = self.Parse_fock(dir+fdir,self.vecl_set_1)
           self.Fock = self.fock_set_1[0]
           print " CALCULATION OF FIRST DERIVATIVES..."
           self.fock1, self.lmoc1 ,self.vecl1 = self.get_efp_fder()
        # Distributed polarizabilities
        elif dpol:
           self.dpos_set_1, self.dpol_set_1 = self.Parse_dpol(dir+dpol)
           self.dpol1 = self.get_dpol_fder()
        # EDS
        elif eds:
           fdir,sdir = eds.split(':')
           self.EDS_set_1 = self.Parse_EDS(dir+fdir)
           self.EDS_set_2 = self.Parse_EDS(dir+sdir)
           self.fEDS = self.get_EDS_fder()
           self.sEDS = self.get_EDS_sder()
        # EFP
        elif efp:
           fdir,sdir = efp.split(':')
           self.EFP_set_1 = self.Parse_EFP(dir+fdir)
           self.EFP_set_2 = self.Parse_EFP(dir+sdir)
           self.fEFP = self.get_EFP_fder()
           self.sEFP = self.get_EFP_sder()
           

    def Parse_dpol(self,dir):
        """"""
        files = glob.glob('%s/*_.efp' % dir)
        files.sort()
        dPos = []
        dPol = []
        r_ref = None
        for file in files:
            r,a=ParseDistributedPolarizabilitiesFromGamessEfpFile(file)
            if r_ref is None: r_ref = r
            #
            r, sim = order(r_ref,r,start=0,lprint=self.__lprint)
            if self.__lprint: print sim
            a = reorder(a,sim)
            dPos.append(r)
            dPol.append(a)
        dPos = array(dPos,float64)
        dPol = array(dPol,float64)
        return dPos, dPol
    
    def get_dpol_fder(self):
        """calculate first derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        first_der_dpol_cart = []
        N = self.dpol_set_1.shape[1]
        for i in range(self.nAtoms*3):
            K = 1 + 4*i
            fd   = (1./12.) *  (\
                       (self.dpol_set_1[K+3] - self.dpol_set_1[K+0] ) \
                 + 8*  (self.dpol_set_1[K+1] - self.dpol_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            first_der_dpol_cart.append( fd )
        first_der_dpol_cart = array(first_der_dpol_cart,float64).reshape(self.nAtoms*3,N*3*3)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_dpol_mode = tensordot(self.L,first_der_dpol_cart,(0,0))

        return first_der_dpol_mode.reshape(self.nAtoms*3-6,N,3,3)
        #return first_der_dpol_cart.reshape(self.nAtoms*3,N,3,3)
    
    def Parse_tran(self,dir,nae,basis):
        """parses CMO transformation matrices from fchk files and transforms
them to LMO space"""
        files = glob.glob('%s/*_.fchk' % dir)
        files.sort()
        Tran = []
        Lmoc = []
        VecL = []
        vec_ref = None
        
        # evaluate transformation matrices and LMO centroids
        for file in files:
            mol  = Read_xyz_file(file,mol=True,
                                 mult=1,charge=0,
                                 name='happy dummy molecule')
            bfs        = PyQuante.Ints.getbasis(mol,basis)
            basis_size = len(bfs)
            natoms= len(mol.atoms)
            SAO   = PyQuante.Ints.getS(bfs)
            dmat = ParseDmatFromFchk(file,basis_size)
            veccmo= ParseVecFromFchk(file)[:nae,:]
            tran, veclmo = get_pmloca(natoms,mapi=bfs.LIST1,sao=SAO,
                                      vecin=veccmo,nae=nae,
                                      maxit=100000,conv=1.0E-10,
                                      lprint=True,
                                      freeze=None)
            if vec_ref is None: vec_ref = veclmo.copy()
            veclmo, sim = order(vec_ref,veclmo,start=0)
            print sim
            check_sim(sim)
            # calculate LMTPs
            camm = coulomb.multip.MULTIP(molecule=mol,
                                         basis=basis,
                                         method='b3lyp',
                                         matrix=dmat,
                                         transition=False,
                                         bonds=None,vec=veclmo)
            camm.camms()
            dma = camm.get()[0]
            lmoc= dma.get_origin()[natoms:]
            #
            Tran.append(tran)
            Lmoc.append(lmoc)
            VecL.append(veclmo)
            #
            out = open(file[:-4]+'vecl','w')
            log = '\n'
            n=1
            for i in range(len(veclmo)):
                log+= ' ***** LMO NO %d ***** \n\n' % (i+1)
                for j in xrange(len(veclmo[0])):
                    log += "%20.10E"  % veclmo[i,j]
                    if not n%5: log+= '\n'
                    n+=1
            log += '\n'
            out.write(log)
            out.close()
        Tran = array(Tran,float64)
        Lmoc = array(Lmoc,float64)
        VecL = array(VecL,float64)
        return Tran, Lmoc, VecL
    
    def Parse_fock(self,dir,Tran):
        """parse directory to look for Fock matrices for ff procedure"""
        files = glob.glob('%s/*_.log' % dir)
        files.sort()
        Fock = []
        for file in files:
            Fock.append(ParseFockFromGamessLog(file,interpol=False))
        Fock = array(Fock,float64)
        # transform from AO to LMO space
        I = Fock.shape[ 0]
        N = Tran.shape[-2]
        Fock_t = zeros((I,N,N),float64)
        for i in xrange(I):
            T = Tran[i]
            F = Fock[i]
            Fock_t[i] = tensordot(T,tensordot(T,F,(1,0)),(1,1))
        return Fock_t
    
    def get_efp_fder(self):
        """calculate first derivatives wrt normal coordinates
of Fock matrix elements given in LMO space. 
Returns ndarray of dimension (nmodes,nmos,nmos)

calculate first derivatives wrt normal coordinates
of LMO centroids. 
Returns ndarray of dimension (nmodes,nmos,3)

calculate first derivatives wrt normal coordinates
of AO to LMO transformation matrix. 
Returns ndarray of dimension (nmodes,nmos,nbasis)"""
        
        first_der_fock_cart = []
        first_der_lmoc_cart = []
        first_der_vecl_cart = []
        for i in xrange(self.nAtoms*3):
            K = 1 + 4*i
            fd1  = (1./12.) *  (\
                       (self.fock_set_1[K+3] - self.fock_set_1[K+0] ) \
                 + 8*  (self.fock_set_1[K+1] - self.fock_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            fd2  = (1./12.) *  (\
                       (self.lmoc_set_1[K+3] - self.lmoc_set_1[K+0] ) \
                 + 8*  (self.lmoc_set_1[K+1] - self.lmoc_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            fd3  = (1./12.) *  (\
                       (self.vecl_set_1[K+3] - self.vecl_set_1[K+0] ) \
                 + 8*  (self.vecl_set_1[K+1] - self.vecl_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
            #
            first_der_fock_cart.append( fd1 )
            first_der_lmoc_cart.append( fd2 )
            first_der_vecl_cart.append( fd3 )
            #
        first_der_fock_cart = array(first_der_fock_cart,float64)
        first_der_lmoc_cart = array(first_der_lmoc_cart,float64)
        first_der_vecl_cart = array(first_der_vecl_cart,float64)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_fock_mode = tensordot(self.L,first_der_fock_cart,(0,0))
        first_der_lmoc_mode = tensordot(self.L,first_der_lmoc_cart,(0,0))
        first_der_vecl_mode = tensordot(self.L,first_der_vecl_cart,(0,0))

        return first_der_fock_mode, first_der_lmoc_mode, first_der_vecl_mode
    
    def get(self,type=0):
        """Get the resulting derivatives according to the type:
type 10: first cartesian derivatives of DMTP
type 20: second normal mode derivatives of DMTP
type 11: first derivatives of EDS interaction energies
type 21: second derivatives of EDS interaction energies
type 12: first mode derivatives of Fock matrix
type 13: first mode derivatives of LMO centroids
type 14: first mode derivatives of AO to LMO transformation matrix
type -1: Fock matrix in LMO (nmos,nmos)
type -2: LMO centroids (nmos)
type -3: transformation matrix from AO to LMO (nmos, nbasis)
"""
        if   type == 12: return self.fock1
        elif type == 13: return self.lmoc1
        elif type == 14: return self.vecl1
        elif type == -1: return self.Fock
        elif type == -2: return self.rLMO
        elif type == -3: return self.U
    
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
        first_der_EDS_cart = array(first_der_EDS_cart,float64)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_EDS_mode = tensordot(self.L,first_der_EDS_cart,(0,0))

        return first_der_EDS_mode
    
    def get_EFP_fder(self):
        """calculate first derivatives wrt normal coordinates
        of interaction energy components from EFP"""
        
        first_der_EFP_cart = []
        for i in range(self.nAtoms*3):
            K = 1 + 4*i
            fd   = (1./12.) *  (\
                       (self.EFP_set_1[K+3] - self.EFP_set_1[K+0] ) \
                 + 8*  (self.EFP_set_1[K+1] - self.EFP_set_1[K+2] ))\
                 /(self.step* self.AngstromToBohr)
                 
            first_der_EFP_cart.append( fd )
        first_der_EFP_cart = array(first_der_EFP_cart,float64)

        ### TRANSFORM FIRST DERIVATIVES TO NORMAL MODE SPACE
        first_der_EFP_mode = tensordot(self.L,first_der_EFP_cart,(0,0))

        return first_der_EFP_mode
    
    def get_EDS_sder(self):
        """calculate second diagonal derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        second_der_EDS_mode = (1./12.) *  (\
            -     (self.EDS_set_2[1+3] + self.EDS_set_2[1+0] ) \
            +16*  (self.EDS_set_2[1+1] + self.EDS_set_2[1+2] ) \
            -30*   self.EDS_set_2[0  ] )                \
            /self.sderiv_step**2

        return second_der_EDS_mode

    def get_EFP_sder(self):
        """calculate second diagonal derivatives wrt normal coordinates 
        of polarizability tensor. Returned in AU unit"""
        
        second_der_EFP_mode = (1./12.) *  (\
            -     (self.EFP_set_2[1+3] + self.EFP_set_2[1+0] ) \
            +16*  (self.EFP_set_2[1+1] + self.EFP_set_2[1+2] ) \
            -30*   self.EFP_set_2[0  ] )                \
            /self.sderiv_step**2

        return second_der_EFP_mode
                
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
        first_der_pol_cart = array(first_der_pol_cart,float64)

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
           for i,dma in enumerate(self.DMA_set):
               #print i,"fder"
               dma.MAKE_FULL()
               ###dma.ChangeOrigin(new_origin_set=self.reference_str)
               dma.ChangeOrigin(new_origin_set=self.reference_origin)
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
        Sder_MMM = zeros((self.nModes,3),float64)
        
        return Fder_DMA, Sder_DMA, Fder_MMM, Sder_MMM

    def _NEW(self,sder=False):
        """calculates FIRST derivatives wrt cartesian coordinates of DMA distribution"""
        Fder = []
        sder_DMA = None
        print "STEP", self.step, self.sderiv_step
        # change origin in the case of CAMM
        if self.camm:
           for i,dma in enumerate(self.DMA_set):
               dma.MAKE_FULL()
               dma.ChangeOrigin(new_origin_set=self.reference_origin)
        s1 = 1./self.step
        ### 5-point formulae
        if self.n_point == 5:
          for i in range(self.nAtoms*3):
            K = 1 + 4*i
            # first derivatves of DMA
            fder_DMA = (1./12.) *  (\
                       (self.DMA_set[K+3] - self.DMA_set[K+0] ) \
                 + 8*  (self.DMA_set[K+1] - self.DMA_set[K+2] ))\
                 * s1
            Fder.append( fder_DMA/self.AngstromToBohr )
            
        # transform to normal mode space  
        Fder_DMA = DMAMadrixMultiply(transpose(self.L),Fder)
        #
        if sder:
           s2 = 1./(self.sderiv_step**2)
           # change origin in the case of CAMM
           if self.camm:
              for i,dma in enumerate(self.sder_DMA_set):
                  dma.MAKE_FULL()
                  dma.ChangeOrigin(new_origin_set=self.reference_origin)
               
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
        #
        return Fder_DMA, sder_DMA
    
    def _Sderivatives(self,sderiv_wrk,mode_id,step):
        """calculates second derivatives wrt NORMAL (!) coordinates
        for selected mode. It requires input files for normal mode displacements!!!
        """
        sder_MM = 0
        # store first and diagonal second derivatives in the lists
        s1 = 1./step
        s2 = 1./(step**2)

        # change origin in the case of CAMM
        if self.camm:
           for i,dma in enumerate(self.sder_DMA_set):
               dma.MAKE_FULL()
               dma.ChangeOrigin(new_origin_set=self.reference_origin)
               
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
        
        # construct DMA_set, Fragments_set and Overall_MM_set
        ###Fragments_set = []
        DMA_set = []
        Overall_MM_set = []
        for file in files:
            dma = ParseDMA ( file, type )
            Fragment = dma.get_pos()
            DMA_set.append( dma )
            ###Fragments_set.append( Fragment )
            
            if type.lower() == 'gaussian':
               ### --- TIP ---
               ### w przyszłości dopisz w klasie FREQ argument w konstruktorze 'type'
               ### aby umożliwić sczytywanie momentów z plików GAMESS lub innych.
               ### Wtedy usuń warunek if w tej funkcji!
               # construct overall multipoles set DMA object
               dipole     = self.Dipole(file)
               quadrupole = self.Quadrupole(file)
               octupole   = self.Octupole(file)
            
               MM = DMA(nfrag=1)
               MM.pos = zeros(3,dtype=float64)
               MM.DMA[1][0] = array(dipole)
               MM.DMA[2][0] = array(quadrupole)
               MM.DMA[3][0] = array(octupole)

               Overall_MM_set.append( MM )
        
        return DMA_set, Overall_MM_set ###array(Fragments_set), as a second thing

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

    def Parse_EFP(self,dir):
        """parses EFP energies"""
        
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
            set.append(ParseEFPInteractionEnergies(dir+file))
            
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
        """calculate IR Harmonic intensities in [kM/mole] 
using L matrix from Gaussfreq and numerical FDip"""
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
        print "MODE", mode_id
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

##############################################################################
if __name__ == '__main__':
   from head import *
   from sys import argv
   import os, glob, pp, PyQuante, utilities, coulomb.multip
   
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
   def CalculateCAMM(basis='6-311++G**',bonds=[],ncpus=1): 
       """calculates CAMMs from density matrix from GAUSSIAN09
using COULOMB.py routines"""
       job_server = pp.Server(ncpus)
       
       pliki_fchk  = glob.glob('./*_.fchk')
       pliki_fchk.sort()
       pliki_log   = glob.glob('./*_.log')
       pliki_log .sort()    
       print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
       for i in range(len(pliki_log)):
           print pliki_log[i], pliki_fchk[i]
       print
       
       for i,file_log in enumerate(pliki_log):
           dma = ParseDMA( file_log, 'gaussian' )
           fragment = dma.get_pos()
           
           ### read atomic symbols
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
           camm = multip.MULTIP(molecule=molecule,
                         basis=basis,
                         #basis='sto-3g',
                         method='b3lyp',
                         matrix=dmat,
                         transition=False,
                         bonds=bonds)
           camm.camms()
           camm.mmms()
           camm.__printMMMs__()
           #CAMM.__printCAMMs__()
       
           dma = camm.get()[0]
           dma.write(file_log[:-4]+'.camm')
           print " Writing file:  :", file_log[:-4]+'.camm'
       print

   def bua(file_fchk,basis,bonds,vec,vec_ref):
       molecule = utilities.Read_xyz_file(file_fchk,mol=True,
                                          mult=1,charge=0,
                                          name='happy dummy molecule')
       
       bfs        = PyQuante.Ints.getbasis(molecule,basis)
       basis_size = len(bfs)
       #print " - basis size= ", basis_size
       dmat = utilities.ParseDmatFromFchk(file_fchk,basis_size)
       def check_sim(l):
           """check the sim list"""
           for x,y in l:
               i=0;j=0
               for a,b in l:
                   if a==x: i+=1
                   if b==y: j+=1
               if (i>1 or j>1): 
                   print " --- !ERROR! --- "
                   break

       ### parse vectors and make Pipek-Mezey transformation
       if vec is not None:
          natoms= len(molecule.atoms)
          SAO   = PyQuante.Ints.getS(bfs)
          print " - ovelrap AO matrix evaluation..."
          nae = vec
          vec = utilities.ParseVecFromFchk(file_fchk)[:nae,:]
          
          print " - Pipek-Mezey localization of %i orbitals..." %nae
          tran, vec = utilities.get_pmloca(natoms,mapi=bfs.LIST1,sao=SAO,
                                           vecin=vec,nae=nae,
                                           maxit=100000,conv=1.0E-10,
                                           lprint=False,
                                           freeze=None)
          vec, sim = utilities.order(vec_ref,vec,start=0)
          print sim
          check_sim(sim)
       ### calculate CAMMs
       print " - multipole integrals in AO basis evaluation..."
       camm = coulomb.multip.MULTIP(molecule=molecule,
                                    basis=basis,
                                    method='b3lyp',
                                    matrix=dmat,
                                    transition=False,
                                    bonds=bonds,vec=vec)
       print " - calculation of %s"%camm.operation
       camm.camms()
       #camm.mmms()
       #camm.__printMMMs__()
       #CAMM.__printCAMMs__()
       
       dma = camm.get()[0]
       dma.write(file_fchk[:-5]+'.camm')
       print " --- Writing file:  :", file_fchk[:-5]+'.camm'    
       return

   def CalculateCAMM_(basis='6-311++G**',bonds=[],ncpus=4,
                      vec=None,vec_ref=None,natoms=11): 
       """calculates CAMMs from density matrix from GAUSSIAN09
using COULOMB.py routines"""

       job_server = pp.Server(ncpus)
       jobs = []
       
       pliki_fchk  = glob.glob('./*_.fchk')
       pliki_fchk.sort()
       print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
       for i in range(len(pliki_fchk)):
           print pliki_fchk[i]
       print
       
       # compute reference vectors
       if vec is not None:
          ref_mol = utilities.Read_xyz_file(pliki_fchk[0],mol=True,
                                            mult=1,charge=0,
                                            name='happy dummy molecule')
          bfs_ref    = PyQuante.Ints.getbasis(ref_mol,basis)
          basis_size = len(bfs_ref)
          sao_ref    = PyQuante.Ints.getS(bfs_ref)
          print " - basis size= ", basis_size
          nae = vec
          vec_ref = utilities.ParseVecFromFchk(pliki_fchk[0])[:nae,:]
          t, vec_ref = utilities.get_pmloca(natoms,mapi=bfs_ref.LIST1,
                                            sao=sao_ref,
                                            vecin=vec_ref,nae=nae,
                                            maxit=100000,conv=1.0E-19,
                                            lprint=False,
                                            freeze=None)
       # submit the jobs!
       i=0
       for file_fchk in pliki_fchk:
           jobs.append( job_server.submit(bua, (file_fchk,basis,bonds,vec,vec_ref),
                                          (Read_xyz_file,ParseDmatFromFchk,ParseDMA,ParseVecFromFchk,get_pmloca,),
                                          ("coulomb.multip","PyQuante",
                                           "PyQuante.Ints","utilities","units","orbloc") ) )
           if (i%4==0 and i!=0): job_server.wait()
           i+=1
       print
       
       for job in jobs:
           result = job()
           
       job_server.print_stats()
       return
                 
   ### calculate camms!
   from sys  import argv
   #bonds = [map(int,(x.split(','))) for x in argv[-1].split('-')]
   #bonds=[(1,0),(2,1),(3,2),(4,2)] # NMA
   bonds = None
   vec = None
   #for i in range(len(bonds)): bonds[i]=tuple(bonds[i])
   #print bonds
   #os.environ['__IMPORT__COULOMB__']=1
   CalculateCAMM_(basis=argv[1],bonds=bonds,vec=vec)

   
