# ------------------------------------------------------------------------- #
#           MINHAENG CHO's SOLVATOCHROMIC FREQUENCY SHIFT MODULE            #
# ------------------------------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
sys.stdout.flush()

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
    
    def __init__(self,pkg="gaussian",nsatoms=('O','H','H'),nstatoms=(),L=[],mode_id=0,
                 solute="",solvent="",solute_origin="com",
                 solute_structure="",mixed=False,
                 overall_MM=False,camm=False,chelpg=False,bsm_file='',
                 fderiv=0,redmass=0,freq=0,structural_change=False,
                 ref_structure=[],solpol=False):
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
        # first DMA derivatives wrt normal mode
        self.fderiv = array(fderiv)
        # reduced masses in helico
        self.redmass = redmass
        # harmonic frequencies in helico
        self.freq = freq
        # benchmark solvent molecule file
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
        # add non-atomic sites if necessary
        self.addSites()
             
        ### construct SOLVENT environment
        self._Solvent()
        ### construct SOLUTE environment
        self._Solute()
        a = self.SOLVENT.copy()
        ### calculate frequency shift!!!
        self.shift=self.FrequencyShift(self.SOLUTE,
                                 self.SOLVENT,
                                 self.solute_structure,
                                 self.solpol)
        ##### estimate structural shifts!!!
        #if structural_change:
        #   self.st=self.StructuralChange(self.fderiv,
        #                         self.SOLVENT,
        #                         self.solute_structure)

    # public methods

    def addSites(self):
        """add non-atomic sites if requested by the input"""
        if len(self.solute_structure) != len(self.ref_structure):
           N_atoms = len(self.solute_structure)
           # parameters DMA 
           X = self.solute.pos[:5]
           # target orientation
           Y = self.solute_structure[:5]
           #rot, rms = RotationMatrix(initial=X,final=Y)
           sup = SVDSuperimposer()
           sup.set(Y,X)
           sup.run()
           rms = sup.get_rms()
           rot, transl = sup.get_rotran()
           transformed = dot(self.ref_structure,rot) + transl
           self.solute_structure = concatenate((self.solute_structure,
                                             transformed[N_atoms:]),axis=0)
           #print self.solute_structure * self.BohrToAngstrom
        return
       
    def FrequencyShift(self,solute=0,solvent=0,solute_structure=0,solpol=0):
        """calculate frequency shift!!!"""
        
        #print solvent.pos*self.BohrToAngstrom
        #print solvent.origin*self.BohrToAngstrom
        solvent.origin = array(solvent.pos)
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
        return shift, solute

    def get_StructuralChange(self,fderiv,solvent,solute_structure,L,ref_structure):
        """estimates the deviations from gas-phase structure"""

        ### superimpose structures
        sup = SVDSuperimposer()
        sup.set(solute_structure,ref_structure)
        sup.run()
        rms = sup.get_rms()
        rot, transl = sup.get_rotran()
        transformed_gas_phase_str = sup.get_transformed()

        ### rotate the fderiv[i]
        for i in fderiv:
            i.pos  =array(ref_structure)
            i.origin  =array(ref_structure)
            i.MAKE_FULL()
            i.Rotate(rot)

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
              # benchmark_solvent_molecule
              benchmark_solvent_molecule, smiec = ParseDMA(self.bsm_file,'coulomb')
              smiec, bsm_pos = ParseDMA(self.bsm_file[:-4]+'log','gaussian')
              del smiec
              benchmark_solvent_molecule.pos = zeros((1,3),dtype=float64)
              benchmark_solvent_molecule.origin = zeros((1,3),dtype=float64)

              X = array(bsm_pos)
              print " -------------------------------"
              print "          RMS analysis"
              print " -------------------------------"
              for water in range(len(self.solvent.pos)/self.sat):
                      # copying of BSM
                      bsm_copy = benchmark_solvent_molecule.copy()
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
                      print "  - solvent rms: %10.6f" % rms
              
        
           elif (self.camm and not self.mixed):
                ### --- transform DMA for each water molecule
                ### --- using solvent_DMA object for one GAS-PHASE
                ### --- water molecule
                # benchmark_solvent_molecule
                benchmark_solvent_molecule, smiec = ParseDMA(self.bsm_file,'coulomb')
                smiec, struct = ParseDMA(self.bsm_file[:-4]+'log','gaussian')
                del smiec
                benchmark_solvent_molecule.pos = array(struct)
                benchmark_solvent_molecule.origin = array(struct)
                # make the solvent DMA object
                SOLVENT = DMA(nfrag=len(self.solvent.pos))
                SOLVENT.pos = array(self.solvent.pos)
                SOLVENT.origin = array(self.solvent.pos)
            
                X = array(benchmark_solvent_molecule.pos)
                print " -------------------------------"
                print "          RMS analysis"
                print " -------------------------------"
                for water in range(len(self.solvent.pos)/self.sat):
                    # copying of BSM
                    bsm_copy = benchmark_solvent_molecule.copy()
                    # superimposition
                    Y = self.solvent.pos[self.sat*water:self.sat*water+self.sat]
                    rot, rms = RotationMatrix(final=Y,initial=X)
                    # rotation
                    bsm_copy.MAKE_FULL()
                    bsm_copy.Rotate(rot)
                    # building DMA object for solvent environment 
                    SOLVENT.DMA[0][self.sat*water:self.sat*water+self.sat] = array(bsm_copy[0])
                    SOLVENT.DMA[1][self.sat*water:self.sat*water+self.sat] = array(bsm_copy[1])
                    SOLVENT.DMA[2][self.sat*water:self.sat*water+self.sat] = array(bsm_copy[2])
                    SOLVENT.DMA[3][self.sat*water:self.sat*water+self.sat] = array(bsm_copy[3])
                    #
                    print "  - solvent rms: %10.6f" % rms

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
        X = self.solute.pos[:4]
        # target orientation
        Y = self.solute_structure[:4]
        #sup.set(Y,X)
        #sup.run()
        #rot, transl = sup.get_rotran()
        rot, rms = RotationMatrix(initial=X,final=Y)
        print "  - solute  rms: %10.6f" % rms  #% sup.get_rms()
       
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
           SOLUTE.pos = array(self.solute_structure)
           SOLUTE.origin = array(self.solute_structure)
           self.update_structure = 0
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
