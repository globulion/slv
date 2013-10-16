# ------------------------------------------------------------------------- #
#           MINHAENG CHO's SOLVATOCHROMIC FREQUENCY SHIFT MODULE            #
# ------------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
from slvcor    import SLVCOR
from hessian   import HESSIAN
import sys, copy
sys.stdout.flush()


__all__ = ['SLV',]
__version__ = '7.1.1'

class SLV(UNITS):
    """\
Frequency shift calculator
----------------------------------------------------------------------
Frequency shift of solute immersed in solvent. Solute is a set of
solvatochromic parameters including:
- distributed charges
- overall multipole moment
- overall multipole moment from partial charges
- multi-site (distributed) multipole moments
Data to be loaded:
- file with solute
- file with solvent
- benchmark solvent molecule file for CAMM/DMA-based computations
Note:
Solvent has to be set in the same order for each molecule!
Eg. : O1H1H1 O2H2H2 O3H3H3 ...
---------------------------------------------------------------------
"""
    
    def __init__(self,pkg="gaussian",nsatoms=('O','H','H'),nstatoms=(),L=[],L_=[],mode_id=0,
                 solute="",solvent="",solute_origin="com",
                 solute_structure="",mixed=False,
                 overall_MM=False,camm=False,chelpg=False,bsm=None,bsm_file='',
                 fderiv=0,sderiv=0,redmass=0,freq=0,structural_change=False,
                 ref_structure=[],solpol=False,gijj=0,lprint=True,suplist=[0,1,2,3],
                 sol_suplist=[0,1,2]):
        self.log = 'None'
        ### modes of action
        # mixed 
        self.mixed = mixed
        # overall multipoles
        self.overall_MM = overall_MM
        # CAMM
        self.camm = camm
        # CHELPG
        self.chelpg = chelpg
        # SolPOL
        self.solpol = solpol
        ### memorials
        # package 
        self.pkg = pkg
        # tuple of solvent atomic numbers
        self.solventAtoms = [ Atom(x) for x in nsatoms ] 
        # tuple of solute atomic numbers
        self.soluteAtoms = [ Atom(x) for x in nstatoms ] 
        # L matrix
        self.L = L
        self.L_= L_
        # mode of interest for frequency shift
        self.mode_id = mode_id
        # solute parameters
        self.solute = solute
        # initial solvent DMA object
        self.solvent = solvent
        # solute origin to be transformed
        self.solute_origin = solute_origin
        # target solute structure (in solvent)
        self.solute_structure = solute_structure
        # first and second DMA derivatives wrt normal mode
        self.fderiv = array(fderiv)
        self.sderiv = array(sderiv)
        # reduced masses in helico
        self.redmass = array(redmass)
        # harmonic frequencies in helico
        self.freq = array(freq)
        # number of modes
        self.nModes = len(freq)
        # cubic anharmonic constants
        self.gijj = array(gijj)
        # benchmark solvent molecule
        self.bsm = bsm
        self.bsm_file = bsm_file
        # number of atoms in solvent molecule
        self.sat = len(nsatoms)
        # number of atoms in solute molecule
        self.stat= len(solute.pos)
        # number of solvent molecules
        self.nmol = len(solvent.pos)/self.sat
        # solvent position array for overall MMs for solvent
        self.solvent_pos = zeros((self.nmol,3),dtype=float64)
        # reference structure
        self.ref_structure = array(ref_structure)
        # solvatochromic polarizability
        self.solpol = solpol
        # print the calculation progress
        self.lprint = lprint
        # superimposition list
        self.__suplist = suplist
        self.__sol_suplist = sol_suplist
        
        ### add non-atomic sites if necessary
        self.addSites()        
             
        ### BSM
        if self.bsm_file:
           self.benchmark_solvent_molecule  = self._bsm()
        if self.bsm is not None:
           self.benchmark_solvent_molecule  = self._bsm()
        ### construct SOLVENT environment
        self._Solvent()
        ### construct SOLUTE environment
        self._Solute()
        a = self.SOLVENT.copy()
        
    # public methods
    
    def eval(self):
        """calculate shifts"""
        self.shift=self.__shift(self.SOLUTE.copy(),
                                 self.SOLVENT.copy(),
                                 self.solute_structure,
                                 self.solpol)
        ##### estimate structural shifts!!!
        #if structural_change:
        #   self.st=self.StructuralChange(self.fderiv,
        #                         self.SOLVENT,
        #                         self.solute_structure)


    def eval_hessian(self,mol,ua_list=None):
        """evaluates Hessian"""
        fderiv_rot, L_rot, transformed_gas_phase_str, rot =\
         self.get_fder_rotated(self.fderiv, self.L, self.solute_structure, self.ref_structure, )
        hess =HESSIAN(fderiv=fderiv_rot,
                      redmass=self.redmass.copy(),
                      freq=self.freq.copy(),
                      solute=self.SOLUTE.copy(),
                      solvent=self.SOLVENT.copy(),
                      mode_id=self.mode_id,
                      L=L_rot.copy(),
                      #L=self.L.copy(),
                      gijj=self.gijj,
                      ua_list=ua_list,
                      mol=mol)
        hess.eval(theory=0,max_iter=500000000,threshold=1e-6)
        return hess
    
    def get_fder_rotated(self,fderiv, L, solute_structure, ref_structure):
        """rotates first derivatives of DMA and the ref structure to target orientation"""
        ### superimpose structures
        sup = SVDSuperimposer()
        sup.set(solute_structure[self.__suplist],ref_structure[self.__suplist])
        sup.run()
        rms = sup.get_rms()
        rot, transl = sup.get_rotran()
        transformed_gas_phase_str = sup.get_transformed()

        ### rotate the fderiv[i]
        fderiv_copy = copy.deepcopy(fderiv)
        for i in fderiv_copy:
            i.pos  =array(solute_structure)
            i.origin  =array(solute_structure)
            i.MAKE_FULL()
            i.Rotate(rot)
        
        ### rotate the eigenvectors
        Lrot = L.reshape(self.stat,3,self.nModes)
        #Lrot = L.reshape(7,3,self.nModes)
        Lrot = tensordot(Lrot,rot,(1,0))              # dimension: nstat,nmodes,3
        Lrot = transpose(Lrot,(0,2,1))                # dimension: nstat,3,nmodes
        Lrot = Lrot.reshape(self.stat*3,self.nModes) # dimension: nstat*3,nmodes
        #Lrot = Lrot.reshape(7*3,self.nModes)
        
        return fderiv_copy, Lrot, transformed_gas_phase_str, rot
    
    def eval_shiftcorr(self,solute_DMA,ua_list=None):
        """evaluates corrections to frequency shifts"""
        
        ### get rotated quantities: L, fderiv, solute DMA
        fderiv_rot, transformed_gas_phase_str, L_rot, rot = self.get_rotated(self.fderiv, self.L, self.solute_structure, self.ref_structure)
        sol = ParseDMA(solute_DMA,'coulomb')
        self.SOLVENT.makeDMAfromFULL()
        if ua_list is not None: sol.MakeUa(ua_list,change_origin=True)
        sol.MAKE_FULL()
        sol.Rotate(rot)
        sol.set_structure(pos=self.solute_structure, origin=self.solute_structure)

        ### calculate the corrections!!!
        corr = SLVCOR(fderiv=fderiv_rot,redmass=self.redmass.copy(),
                      freq=self.freq.copy(),L=L_rot.copy(),solute=sol.copy(),
                      solvent=self.SOLVENT.copy(),
                      mode_id=self.mode_id,gijj=self.gijj,ua_list=ua_list)
        corr.eval()
        corr.shift = self.shift[0]
        # set the corrected shift
        self.shift_corr = self.shift[0] + corr.corr
        self.log = self.log[:-1]+str(corr)
        return
    
    def _bsm(self):
        """benchmark_solvent_molecule extractor"""
        if self.bsm is not None: benchmark_solvent_molecule = self.bsm.copy()
        else:                    benchmark_solvent_molecule = ParseDMA(self.bsm_file,'coulomb')
        return benchmark_solvent_molecule
    
    def addSites(self):
        """add non-atomic sites if requested by the input"""
        if len(self.solute_structure) != len(self.ref_structure):
           N_atoms = len(self.solute_structure)
           # parameters DMA
           X = self.solute.pos[self.__suplist]
           # target orientation
           Y = self.solute_structure[self.__suplist]
           #rot, rms = RotationMatrix(initial=X,final=Y)
           sup = SVDSuperimposer()
           sup.set(Y,X)
           sup.run()
           rms = sup.get_rms()
           rot, transl = sup.get_rotran()
           transformed = dot(self.ref_structure,rot) + transl
           self.solute_structure = concatenate((self.solute_structure,
                                             transformed[N_atoms:]),axis=0)
           #print  "AAAAA"
           #print self.solute_structure * self.BohrToAngstrom
           #print self.solute_structure * self.BohrToAngstrom
        return
       
    def __shift(self,solute=0,solvent=0,solute_structure=0,solpol=0):
        """calculate frequency shift!!!"""
        
        #print solvent.pos*self.BohrToAngstrom
        #print solvent.origin*self.BohrToAngstrom
        solvent.origin = array(solvent.pos)
        #print "SOLVENT"
        #print solvent.pos * self.BohrToAngstrom
        shift = FrequencyShift(solute=solute,
                               solvent=solvent,
                               solute_structure=solute_structure)
        
        solvent.MAKE_FULL()
        solvent.MakeTraceless()
        if self.solpol:
           solpol = dot(transpose(self.rot),dot(solpol,self.rot))
           shift_pol = FrequencyShiftPol(solvent,solpol,self.r_origin_target)
           print " SHIFT-POL [cm-1]: ",shift_pol
        out = "\n INTERACTION ENERGY TERMS [cm-1]"
        out+= '\n'
        self.log = Emtp.log
        return shift, solute

    def get_rotated(self,fderiv, L, solute_structure, ref_structure, suplist):
        """rotates first derivatives of DMA and the ref structure to target orientation"""
        ### superimpose structures
        sup = SVDSuperimposer()
        sup.set(solute_structure[suplist],ref_structure[suplist])
        sup.run()
        rms = sup.get_rms()
        rot, transl = sup.get_rotran()
        transformed_gas_phase_str = sup.get_transformed()

        ### rotate the fderiv[i]
        fderiv_copy = copy.deepcopy(fderiv)
        #for i in fderiv_copy:
        #    i.pos  =array(ref_structure)
        #    i.origin  =array(ref_structure)
        #    i.MAKE_FULL()
        #    i.Rotate(rot)
        fderiv_copy.pos  =array(ref_structure)
        fderiv_copy.origin  =array(ref_structure)
        fderiv_copy.MAKE_FULL()
        fderiv_copy.Rotate(rot)
            
        ### rotate the eigenvectors
        Lrot = L.reshape(self.stat,3,self.nModes)
        #Lrot = L.reshape(7,3,self.nModes)
        Lrot = tensordot(Lrot,rot,(1,0))              # dimension: nstat,nmodes,3
        Lrot = transpose(Lrot,(0,2,1))                # dimension: nstat,3,nmodes
        Lrot = Lrot.reshape(self.stat*3,self.nModes) # dimension: nstat*3,nmodes
        #Lrot = Lrot.reshape(7*3,self.nModes)
        
        return fderiv_copy, transformed_gas_phase_str, Lrot, rot
    
    def get_StructuralChange(self,fderiv,solvent,solute_structure,L,ref_structure, suplist):
        """estimates the deviations from gas-phase structure"""

        ### superimpose structures and rotate first derivatives
        fderiv, transformed_gas_phase_str, Lrot, rot = self.get_rotated(fderiv, L, solute_structure, ref_structure, suplist)

        ### calculate displacement vector dQ in GC normal mode space
        dQ = []
        for i in range(len(fderiv)):
            dq =  FrequencyShift(solute=fderiv[i],
                                 solvent=solvent,
                                 solute_structure=transformed_gas_phase_str)\
                                 /self.HartreePerHbarToCmRec*sqrt(self.redmass[i]*self.AmuToElectronMass)
            dq /= - self.redmass[i]*self.AmuToElectronMass*\
                  (self.freq[i]* self.CmRecToHz * self.HzToAuAngFreq)**2
            # switch to normal mode unit [Bohr*me^1/2]
            dq *= sqrt(self.redmass[i]*self.AmuToElectronMass)
            dQ.append(dq)
            
        dQ = array(dQ)
        
        ### delete wrongly described vibrations
        dQ[-1] = zeros(5,dtype=float64)
        dQ[-2] = zeros(5,dtype=float64)
        #dQ[-3] = zeros(5,dtype=float64)
        #dQ[-4] = zeros(5,dtype=float64)
        #dQ[-5] = zeros(5,dtype=float64)
        
        ### calculate cartesian structure [in Bohr]
        dX = []
        for i in range(5):
            dx = (dot(L,dQ[:,i])).reshape(12,3)
            dX.append(dx)
        print "Displacements in Angstrom"
        #print array(dX)*self.BohrToAngstrom
        # calculate structure
        pred_structure = array(dX) + ref_structure#transformed_gas_phase_str#ref_structure
        print "Predicted structure in Angstrom"
        #print pred_structure*self.BohrToAngstrom
        for i in range(5):
            rot, rms = RotationMatrix(initial=pred_structure[i],final=solute_structure)
            print " RMS CALC/SOL : %16.4f"% rms
            
        for i in range(5):
            rot, rms = RotationMatrix(initial=pred_structure[i],final=ref_structure)
            print " RMS CALC/REF: %16.4f"% rms
            
        rot, rms = RotationMatrix(initial=solute_structure,final=ref_structure)
        print " RMS  SOL/REF: %16.4f"% rms
        return pred_structure
    
    def __repr__(self):
        """print nice table of results"""
        return str(self.log)
    
    # private methods
    
    ## ------------------------------ ##
    ##         S O L V E N T          ##
    ## ------------------------------ ##
           
    def _Solvent(self):
        """constructs solvent environment for frequency shift calculations"""
        if self.pkg.lower()=="gaussian":
           ### construct DMA for solvent environment from charges
           if (self.overall_MM and self.mixed):
              SOLVENT = DMA(nfrag=self.nmol)
              SOLVENT.pos= zeros((self.nmol,3),dtype=float64)
              for i in range(self.nmol):
                  ### calculate the center of mass
                  r_com = zeros(3,dtype=float64)
                  mass_sum = 0
                  for atom in range(self.sat):
                      #r_com   += self.mass[ self.nsatoms[atom] ] * self.solvent.pos[i*(self.sat)+atom]
                      #mass_sum+= self.mass[ self.nsatoms[atom] ]
                      r_com   += self.solventAtoms[atom].mass * self.solvent.pos[i*(self.sat)+atom]
                      mass_sum+= self.solventAtoms[atom].mass
                  r_com/=mass_sum
                  SOLVENT.pos[i] = r_com

                  ### calculate molecular moments
                  mu   = zeros((3),dtype=float64)
                  quad = zeros((3,3),dtype=float64)
                  oct  = zeros((3,3,3),dtype=float64)
                  for atom in range(self.sat):
                      r     = self.solvent.pos[i*(self.sat)+atom] - r_com
                      qatom = self.solvent[0][i*(self.sat)+atom]
                      ### calculate dipole moment
                      mu   += qatom * r
                      ### calculate quadrupole moment
                      quad += qatom * outer (r,r)
                      ### calculate octupole moment
                      oct  += qatom * outer( r, outer (r,r) ).reshape(3,3,3)
           
                  SOLVENT.origin = array(SOLVENT.pos)
                  ### set the molecular moments into the DMA solvent object
                  SOLVENT.DMA[1][i] = mu
           
                  SOLVENT.DMA[2][i,0] = quad[0,0]
                  SOLVENT.DMA[2][i,1] = quad[1,1]
                  SOLVENT.DMA[2][i,2] = quad[2,2]
                  SOLVENT.DMA[2][i,3] = quad[0,1]
                  SOLVENT.DMA[2][i,4] = quad[0,2]
                  SOLVENT.DMA[2][i,5] = quad[1,2]
           
                  SOLVENT.DMA[3][i,0] = oct[0,0,0]
                  SOLVENT.DMA[3][i,1] = oct[1,1,1]
                  SOLVENT.DMA[3][i,2] = oct[2,2,2]
                  SOLVENT.DMA[3][i,3] = oct[0,0,1]
                  SOLVENT.DMA[3][i,4] = oct[0,0,2]
                  SOLVENT.DMA[3][i,5] = oct[0,1,1]
                  SOLVENT.DMA[3][i,6] = oct[1,1,2]
                  SOLVENT.DMA[3][i,7] = oct[0,2,2]
                  SOLVENT.DMA[3][i,8] = oct[1,2,2]
                  SOLVENT.DMA[3][i,9] = oct[0,1,2]
                  
           elif (self.mixed and not self.overall_MM):
              SOLVENT = self.solvent.copy()
              
           elif (not self.mixed and self.overall_MM):
              ### initialize SOLVENT object 
              SOLVENT = DMA(nfrag=self.nmol)
              SOLVENT.pos= zeros((self.nmol,3),dtype=float64)
              SOLVENT.origin = zeros((self.nmol,3),dtype=float64)  
                          
              ### calculate the target origins and positions of SOLVENT multipoles
              for i in range(self.nmol):
                  r_com = zeros(3,dtype=float64)
                  mass_sum = 0
                  for atom in range(self.sat):
                      r_com   += self.solventAtoms[atom].mass * self.solvent.pos[i*(self.sat)+atom]
                      mass_sum+= self.solventAtoms[atom].mass
                  r_com/=mass_sum
                  SOLVENT.pos[i] = r_com
                  SOLVENT.origin[i] = r_com

              ### calculate molecular moments
              X = array(self.benchmark_solvent_molecule.pos.copy())
              if self.lprint:
                 print " -------------------------------"
                 print "          RMS analysis"
                 print " -------------------------------"
              for water in range(len(self.solvent.pos)/self.sat):
                      # copying of BSM
                      bsm_copy = self.benchmark_solvent_molecule.copy()
                      # superimposition
                      Y = self.solvent.pos[self.sat*water:self.sat*water+self.sat]
                      rot, rms = RotationMatrix(final=Y,initial=X)
                      # rotation
                      bsm_copy.MAKE_FULL()
                      bsm_copy.Rotate(rot)
                      # building DMA object for solvent environment 
                      SOLVENT.DMA[0][water] = array(bsm_copy[0])
                      SOLVENT.DMA[1][water] = array(bsm_copy[1])
                      SOLVENT.DMA[2][water] = array(bsm_copy[2])
                      SOLVENT.DMA[3][water] = array(bsm_copy[3])
                      #
                      if self.lprint: print "  - solvent rms: %10.6f" % rms
              
        
           elif (self.camm and not self.mixed):
                ### --- transform DMA for each water molecule
                ### --- using solvent_DMA object for one GAS-PHASE
                ### --- water molecule
                nmol = len(self.solvent.get_pos())/self.sat
                X = self.benchmark_solvent_molecule.get_origin()
                nscnt = len(X)
                X = X[self.__sol_suplist]
                # make the solvent DMA object
                SOLVENT = DMA(nfrag=nscnt*nmol)
                if self.lprint:
                   print " -------------------------------"
                   print "          RMS analysis"
                   print " -------------------------------"
                Ys = self.solvent.get_pos()
                orig = None
                for mol in xrange(nmol):
                    # copying of BSM
                    bsm_copy = self.benchmark_solvent_molecule.copy()
                    # superimposition
                    Y = Ys[self.sat*mol:self.sat*mol+self.sat]
                    sup = SVDSuperimposer()
                    sup.set(Y,X)
                    sup.run()
                    rms = sup.get_rms()
                    rot, transl = sup.get_rotran()
                    transformed = dot(self.benchmark_solvent_molecule.get_origin(),rot) + transl
                    #Yp = concatenate((Y,transformed[self.sat:]),axis=0)
                    Yp = transformed.copy()
                    if orig is not None:
                       orig = concatenate((orig,Yp),axis=0)
                    else: orig = Yp.copy()
                    # rotation
                    bsm_copy.MAKE_FULL()
                    bsm_copy.Rotate(rot)
                    # building DMA object for solvent environment 
                    SOLVENT.DMA[0][nscnt*mol:nscnt*mol+nscnt] = array(bsm_copy[0])
                    SOLVENT.DMA[1][nscnt*mol:nscnt*mol+nscnt] = array(bsm_copy[1])
                    SOLVENT.DMA[2][nscnt*mol:nscnt*mol+nscnt] = array(bsm_copy[2])
                    SOLVENT.DMA[3][nscnt*mol:nscnt*mol+nscnt] = array(bsm_copy[3])
                    #
                    if self.lprint: print "  - solvent rms: %10.6f" % rms
                #
                SOLVENT.set_structure(pos=orig,equal=True)

           else:
                  SOLVENT = self.solvent.copy()
                  
           ### memorize the solvent environment 
           self.SOLVENT = SOLVENT

    ## ---------------------------- ##
    ##         S O L U T E          ##
    ## ---------------------------- ##

    def _Solute(self):
      """constructs solvent environment for frequency shift calculations"""
      if self.pkg.lower()=="gaussian":
        ### from SVD analysis withdraw rotation matrix
        #sup = SVDSuperimposer()
        # parameters DMA
        X = self.solute.pos[self.__suplist]
        # target orientation
        Y = self.solute_structure[self.__suplist]
        #sup.set(Y,X)
        #sup.run()
        #rot, transl = sup.get_rotran()
        rot, rms = RotationMatrix(initial=X,final=Y)
        if self.lprint: print "  - solute  rms: %10.6f" % rms  #% sup.get_rms()
       
        if (not self.camm and not self.chelpg):
        #if 1:
           ### calculate origin for SOLUTE!
           r_origin = zeros(3,dtype=float64)
           r_origin_target = zeros(3,dtype=float64)
           mass_sum = 0
           if   self.solute_origin.lower() == 'com':
                for atom in range(self.stat):
                    r_origin+= self.soluteAtoms[atom].mass * self.ref_structure[atom]
                    r_origin_target+= self.soluteAtoms[atom].mass * self.solute_structure[atom]
                    mass_sum+= self.soluteAtoms[atom].mass
                r_origin/=mass_sum
                r_origin_target/=mass_sum
                            
           elif self.solute_origin.lower() == 'cos':
                for atom in range(self.stat):
                    r_origin+= self.solute.DMA[0][atom] ** 2 * self.ref_structure[atom]
                    r_origin_target+= self.solute.DMA[0][atom] ** 2 * self.solute_structure[atom]
                    mass_sum+= self.solute.DMA[0][atom] ** 2
                r_origin/=mass_sum
                r_origin_target/=mass_sum

           elif self.solute_origin.lower() == 'coe':
                vec = self.L[:,self.mode_id].reshape(self.stat,3)
                for atom in range(self.stat):
                    r_origin+= sum(vec[atom]**2)  * self.ref_structure[atom] / sum(vec**2)
                    vec_rot = dot(vec,rot)
                    r_origin_target+= sum(vec_rot[atom]**2)  * self.solute_structure[atom] / sum(vec_rot**2)
                
           elif self.solute_origin.lower().startswith('at'):
                at,atom1,atom2,degree = self.solute_origin.lower().split(',')
                atom1 = int(atom1) - 1
                atom2 = int(atom2) - 1
                degree= float64(degree)
                r_origin = degree* self.ref_structure[atom2] -\
                    (degree - 1) * self.ref_structure[atom1]
                r_origin_target = degree* self.solute_structure[atom2] -\
                           (degree - 1) * self.solute_structure[atom1]
           else:
                # provide your own origin in Angstroms for gas-phase structure!
                r_origin = array(self.solute_origin.split(','),dtype=float64) * self.AngstromToBohr
                
           self.r_origin = r_origin
           self.r_origin_target = r_origin_target
           self.rot = rot

        ### make solute DMA object
        #if (not self.overall_MM and not self.camm and not self.chelpg):
        if (self.chelpg and self.mixed):
           SOLUTE = DMA(nfrag=1)
           SOLUTE.pos = zeros((1,3),dtype=float64)
       
           ### set the SOLUTE origin     
           SOLUTE.pos[0] = array([r_origin_target])
          
           ### compute molecular solvatochromic moments
           mu   = zeros((3),dtype=float64)
           quad = zeros((3,3),dtype=float64)
           oct  = zeros((3,3,3),dtype=float64)
           
           for atom in range(self.stat):
               r     = self.solute_structure[atom] - r_origin_target
               qatom = self.solute[0][atom]
               ### calculate dipole moment
               mu   += qatom * r
               ### calculate quadrupole moment
               quad += qatom * outer (r,r)
               ### calculate octupole moment
               oct  += qatom * outer( r, outer (r,r) ).reshape(3,3,3)

           ### set the molecular moments into the DMA solvent object
           SOLUTE.DMA[1][0] = mu
           
           SOLUTE.DMA[2][0,0] = quad[0,0]
           SOLUTE.DMA[2][0,1] = quad[1,1]
           SOLUTE.DMA[2][0,2] = quad[2,2]
           SOLUTE.DMA[2][0,3] = quad[0,1]
           SOLUTE.DMA[2][0,4] = quad[0,2]
           SOLUTE.DMA[2][0,5] = quad[1,2]
           
           SOLUTE.DMA[3][0,0] = oct[0,0,0]
           SOLUTE.DMA[3][0,1] = oct[1,1,1]
           SOLUTE.DMA[3][0,2] = oct[2,2,2]
           SOLUTE.DMA[3][0,3] = oct[0,0,1]
           SOLUTE.DMA[3][0,4] = oct[0,0,2]
           SOLUTE.DMA[3][0,5] = oct[0,1,1]
           SOLUTE.DMA[3][0,6] = oct[1,1,2]
           SOLUTE.DMA[3][0,7] = oct[0,2,2]
           SOLUTE.DMA[3][0,8] = oct[1,2,2]
           SOLUTE.DMA[3][0,9] = oct[0,1,2]

           self.solute_structure = array([r_origin_target])
           
        elif self.overall_MM:
           SOLUTE = self.overall_MM.copy()
           SOLUTE.MAKE_FULL()
           SOLUTE.ChangeOrigin( array([r_origin]) )
           #r_origin = dot(r_origin,rot)
           SOLUTE.pos = array([r_origin_target])
           SOLUTE.origin = array([r_origin_target])
           SOLUTE.MAKE_FULL()
           SOLUTE.Rotate(rot)
           self.solute_structure = array([r_origin_target])
          
        elif (self.chelpg and not self.mixed):
           SOLUTE = self.chelpg.copy()
           SOLUTE.pos = array(self.solute_structure)
           SOLUTE.origin = zeros((1,len(SOLUTE.pos)))
           #print SOLUTE
              
        elif (self.camm):
           SOLUTE = self.camm.copy()
           SOLUTE.MAKE_FULL()
           SOLUTE.Rotate(rot)
           SOLUTE.set_structure(pos=self.solute_structure,equal=True)
           self.update_structure = 0
           #print "BBBB"
           #print self.solute_structure * self.BohrToAngstrom
           if self.update_structure:
              sup = SVDSuperimposer()
              sup.set(self.solute_structure,self.ref_structure)
              sup.run()
              rms = sup.get_rms()
              rot, transl = sup.get_rotran()
              transformed_gas_phase_str = sup.get_transformed()
              SOLUTE.origin = array(transformed_gas_phase_str)
              SOLUTE.pos = array(transformed_gas_phase_str)
              #SOLUTE.MAKE_FULL()
              #SOLUTE.ChangeOrigin(new_origin_set=self.solute_structure)
          
        #print "solute pos ", SOLUTE.pos

      ### memorize the solute
      self.SOLUTE = SOLUTE


### tests
if __name__=="__main__": pass
     #ShiftFromMD(pkg="amber",charges="z.prmtop",trajectory="ZMD6.mdcrd")
