# --------------------------------------- #
#             FITTING MODULE              #
# --------------------------------------- #

from numpy import array, float64, zeros, sqrt, \
                  sum, transpose, linalg, ones, \
                  dot
from units import *
from dma   import *
from diff  import *
from utilities import QMFile, ModifyStruct
import os, glob

__all__ = ['FIT','LxFit',]
__version__ = '3.2.56'


class LxFit(UNITS):
    """
Empirical Distributed Solvatochromic Charge Method
--------------------------------------------------

Usage:
"""
    def __init__(self,atoms,dat_path='./dat',ext='xyz',
                      unique_idx=None,rings=None,exclude_atoms=0):
        self.__atoms = atoms
        self.__dat_path = dat_path
        self.__ext = ext
        self.__unique_idx = unique_idx
        self.__exclude_atoms = exclude_atoms
        if rings is not None:
           for ring in rings:
               ring.insert(2,0)
        #       self.__atoms += ring[3]
        self.__rings = rings

    def get_par(self):
        """return solvatochromic parameters"""
        return self.__par

    def makeRings(self):
        """create obvolutkis in the benchmark set of clusters"""
        # copy files
        path = self.__dat_path
        os.system('mkdir ./new')
        #os.system('cp -r %s/* ./new ' % path)
        #os.system('mv %s ./old' % path)
        #os.system('mv ./new %s' % path)
        # create rings
        self._rings(self.__dat_path,self.__ext,self.__atoms,self.__rings,self.__exclude_atoms)
        return

    def _rings(self,path,ext,atoms,rings,exclude_atoms):
        """add rings"""
        files = glob.glob(path+'/'+'*.'+ext)
        f = QMFile()
        for file in files:
            f.open(file,format='xyz')
            pos = f.get_pos()

            S = ModifyStruct(pos[:,:3])
            for ring in self.__rings:
                S.makeRing(*ring)
            f.insert(S.ring[1:],id=atoms)
            new_name = f.get_file_name()+'_ring.xyz'
            new_name = new_name.split('/')[-1]
            f.write(new_name)
            os.system('mv %s ./new/%s' % (new_name,new_name))

        self.__exclude_atoms += len(S.ring)-1
        self.__atoms += len(S.ring)-1
        return

    def check(self):
        """check"""
        print " ---------------------------------------- check -----------------------------------------------"
        for i,j in enumerate(self.__shifts):
            print "%50s %10.2f %10.2f" % (self.__file_names[i].ljust(50),j,self.__ref_shifts[i])
        print " R²        : %10.4f" % self.get_r2() 
        print " Sum of l_x: %10.4f" % sum(self.__par)
        print " ---------------------------------------- check -----------------------------------------------"
        return 

    def get_r2(self):
        """return r² of fitting"""
        data_av = average(self.__shifts)
        sse = sum((self.__ref_shifts-self.__shifts)**2)
        sst = sum((self.__shifts-data_av)**2)
        return 1 - sse/sst 

    def eval(self):
        """perform fitting"""
        # parse benchmark data
        self._parse(self.__dat_path,self.__ext)
        # construct potential matrix
        self._phi(self.__atoms,self.__exclude_atoms)
        # symmetry constraints
        if self.__unique_idx is not None:
           self.__phi = self._contract_phi(self.__phi_start, self.__unique_idx)
        else: self.__phi = self.__phi_start.copy()
        # evaluate lx parameters
        pt = transpose(self.__phi)
        C = dot(pt,self.__phi)
        
        B = ones((C.shape[0]+1,C.shape[1]+1),float64)
        B[:-1,:-1] = C; B[-1,-1] = 0.0
        V = dot(pt,self.__freq-self.__ref_freq)
        VV = zeros(len(V)+1,float64)
        VV[:-1] = V

        self.__par = dot(linalg.inv(B), VV )[:-1]
        self.__shifts = dot(self.__phi,self.__par)
        return

    # protected

    def _phi(self,atoms,exclude_atoms):
        """
construct electrostatic potential matrix
phi_ij = phi_i(r_j)"""
        phi = zeros((self.__k,self.__n),float64)
        for i in xrange(self.__k):
            pos = self.__pos[i]
            chg = self.__chg[i]
            # atomic sites
            for j in xrange(self.__n):
                v = 0.0
                rx = pos[j]
                # sum over solvent sites
                for k in xrange(exclude_atoms,pos.shape[0]):
                    qk = chg[k]
                    rk = pos[k]
                    r  = sqrt(sum( (rk-rx)**2 ))
                    v += qk/r
                phi[i,j] = v
        self.__phi_start = phi
        return

    def _parse(self,path,ext):
        """parse benchmark frequencies, structures and charges"""
        files = glob.glob(path+'/'+'*.'+ext)
        f = QMFile()
        freq = list(); pos = list(); chg = list(); file_names = list()
        for file in files:
            f.open(file,format='xyz')
            w = float64(f.get_misc().split(':')[-1])
            q = f.get_pos()

            freq.append(w)
            pos.append(q[:,:3])
            chg.append(q[:,3])
 
            file_names.append(f.get_file_name())

        f.open(path+'/'+'__ref__',format='xyz')

        ref_freq = float64(f.get_misc().split(':')[-1])
        ref_pos = f.get_pos()[:,:3]
        ref_chg = f.get_pos()[:,3]

        self.__freq = array(freq,float64)
        self.__pos  = pos
        self.__chg  = chg

        self.__ref_freq = ref_freq
        self.__ref_pos = ref_pos
        self.__ref_chg = ref_chg

        self.__k = len(freq)
 
        # construct solvatochromic model
        self.__model = ref_pos[:self.__atoms]
        self.__n = len(self.__model)
        #self.__m = 0
        self.__ref_shifts = self.__freq - self.__ref_freq
        self.__file_names = file_names
        return

    def _contract_phi(self,Phi, unique_list):
        """contract electrostatic potential matrix to constrain symmetry of l_x parameters"""
        G = []; n_gr = len(unique_list)                                               
        for i in unique_list: G+=i
        N = Phi.shape[1]    # total number of parameters
        N_= N-len(G)+n_gr   # number of unique parameters
        Phi_new = []
        Phi_contr = zeros((Phi.shape[0],N_),float64)
        for i in range(N):
            if i not in G: Phi_new.append(Phi[:,i])
        Phi_new= transpose(array(Phi_new,float64))
                                                                                     
        for g in range(n_gr):     Phi_contr[:,g] = sum(Phi[:,unique_list[g]],axis=1)
        for i in range(n_gr,N_):  Phi_contr[:,i] = Phi_new[:,i-n_gr]
        return Phi_contr



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
            print




            
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
