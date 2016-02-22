#!*-* coding=utf-8 *-*
# -------------------------------------------------------------------------- #
#            SOLVATOCHROMIC EFFECTIVE FRAGMENT POTENTIAL MODULE              #
# -------------------------------------------------------------------------- #

import sys, copy, os, re, math, glob, PyQuante.Ints, coulomb.multip, libbbg.qm.clemtp,\
       numpy, libbbg.units, libbbg.utilities, diff, solvshift.shftex, solvshift.shftce,\
       solvshift.solpol, solvshift.solpol2, solvshift.efprot, solvshift.exrep,\
       solvshift.slvpar
sys.stdout.flush()

__all__     = [ 'EFP', 'FragFactory', ]
__version__ = '1.0.3'
__author__  = 'Bartosz Błasiak'

class EFP(object, libbbg.units.UNITS):
    """
 =============================================================================================
                            EFFECTIVE FRAGMENT POTENTIAL METHOD                            
 =============================================================================================
 

 Usage example:

 from solvshift.solefp import EFP 
 efp = EFP(elect=elect, pol=pol, rep=rep, corr=corr, disp=disp, ea=ea, freq=True, mode=mode,
           ccut=ccut, pcut=pcut, ecut=ecut,
           cunit=True)
 efp.set(xyz,ind,nmol,bsm,supl_copy,reord_copy)
 efp.eval(num=True, lwrite=False, step=0.006, theory=2)

 ---------------------------------------------------------------------------------------------
 
 
 The global settings for a computation:
 ccut         - Coulomb cutoff
 pcut         - Polarization cutoff
 pairwise_all - count all pairwise interactions (otherwise count
                only these that include central molecule, False 
                when frequency shift mode is switched on in freq and mode)
 elect        - evaluate electrostatics
 pol          - evaluate polarization
 rep          - evaluate exchange-repulsion
 ct           - evaluate charge-transfer
 disp         - evaluate dispersion
 all          - evaluate all these interactions
 =============================================================================================
                                                                 Last Revision: 7 Jun 2015
"""
    def __init__(self,
                      elect=True , pol=False , rep=False, ct=False, disp=False, all=False,
                      ccut=None  , pcut=None , ecut=None,
                      ea=True    , corr=False,
                      nlo=False  , freq=False,
                      cunit=False, mode=None ,                                           ):
        self._create()
        self.__cunit = cunit
        #
        self.__ccut = 1.E+10 if ccut is None else ccut
        self.__pcut = 1.E+10 if pcut is None else pcut
        self.__ecut = 1.E+10 if ecut is None else ecut
        #
        if all:
           self.__eval_elect = True; self.__eval_pol   = True
           self.__eval_rep   = True; self.__eval_ct    = True
           self.__eval_disp  = True; self.__eval_corr  = False
           self.__eval_elect_ea = True
        else:
           self.__eval_elect = elect; self.__eval_pol   = pol
           self.__eval_rep   = rep  ; self.__eval_ct    = ct
           self.__eval_disp  = disp ; self.__eval_corr  = corr
           self.__eval_elect_ea = ea
        #
        if freq: 
           self.__pairwise_all = False; self.__eval_freq = True
           self.__mode = mode
        else: self.__pairwise_all = True
        #
        if nlo: 
           self.__eval_nlo = True
        return
    
    def __call__(self, lwrite=False, num=False, step=0.006, theory=0):
        """Evaluate the property"""
        self.eval(lwrite, num, step, theory)
        return
    
    # P U B L I C
    # OPERATIONAL METHODS

    def eval(self, lwrite=False, num=False, step=0.006, theory=0, remove_clashes=False, rcl_algorithm='remove_by_name'):
        """
 ----------------------------------------------------------------------------------------------------------------------
 Evaluate the properties. Currently available ways of usage:

 Task                                                               Relevant keywords                                   
 ----------------------------------------------------------------------------------------------------------------------
 Interaction energies                                               lwrite
 Frequency shifts                                                   lwrite
  a)  analytically                                                  remove_clashes, rcl_algorithm
  b)  analytically with removal of EFP clashes for polarization
  c)  numerically (A,B)                                             num, step, theory
 ----------------------------------------------------------------------------------------------------------------------
 Comments:

 (A) removal of EFP clashes is not implemented
 (B) dispersion is not added yet
 ----------------------------------------------------------------------------------------------------------------------
 Usage:

 efp.eval(lwrite=False, num=False, step=0.006, theory=0, remove_clashes=False, rcl_algorithm='remove_by_name')
 
 Options:
 lwrite                           Print additional information (superimposition data, molecules in each layer)
 num                              Evaluate frequency shifts numercially
 step                             Step for numerical first-derivatives. Given in Angstroms
 theory                           The level of SolX theory (0 or 2). Relevant if num=True
 remove_clashes                   Remove clashes between EFP residues (only affects polarization frequency shifts!)
 rcl_algorithm                    Algorithm for removing EFP residues when remove_clashes=True
 ----------------------------------------------------------------------------------------------------------------------
                                                                                            Last Revision: 13 Aug 2015
 """
        self._eval(lwrite, num, step, theory, remove_clashes, rcl_algorithm)
        return

    def eval_dma(self, dma, lwrite=False):
        """
 ----------------------------------------------------------------------------------------------------------------------
 Evaluate the properties due to electrostatic environment (relevant only in central mode)
 ----------------------------------------------------------------------------------------------------------------------

 This is designed only to frequency shift evaluation. 

 Usage:

 efp.eval_dma(dma, lwrite=False)

 Options:

 dma                              DMA distribution of the environment. Can be given in two forms:
                                     a) libbbg.dma.DMA object with primitive moments (not traceless!)
                                     b) numpy.ndarray of shape (ndma,4) where [:,:3] slice contains XYZ in Bohr
                                        whereas [:,3] slice are the corresponding charges (A.U.)
 lwrite                           Print additional information
 ----------------------------------------------------------------------------------------------------------------------
 Comments:

 Cut-offs are not implemented for this layer yet. Therefore whole DMA distribution will be taken into account.
 ----------------------------------------------------------------------------------------------------------------------
                                                                                           Last Revision: 13 Aug 2015
 """
        self._eval_dma(dma, lwrite)
        return

    def close(self):
        """
 Close the debug files if any (if lwrite in eval was True it has to be called 
 at the end of working of EFP object). Otherwise self.__debug file won't be closed.
"""
        if self.__debug is not None: 
           self.__debug.close()
           os.system('mv __debug_file__ __73nfbod783748x__')
           debg = open('__debug_file__','w')
           temp = open('__73nfbod783748x__'); d = temp.readlines() ; temp.close()
           os.system('rm __73nfbod783748x__')
           debg.write('%d\n\n' % len(d))
           debg.write(''.join(d))
           debg.close()
        return

    def transform_vec_to_bsm_orientation(self, vec, t=False):
        """
 Transforms the vectors to the BSM molecule orientation (if __bsm[0] is copied in eval).
 vec can be a single array of size 3 or a set of N vectors in array of size N x 3.
 Usage:
       vec = transform_vec_to_bsm_orientation(vec, t=False)
  if t=True: then perform also translation; otherwise only rotation is performed.
 See also: get_rot_transl, get_pos_c, get_pos_ref
"""
        rot, transl = self.get_rot_transl()
        vec = numpy.dot(vec, rot)
        if t: vec += transl
        return vec

    def __repr__(self):
        """print the status"""
        log = '\n'
        log+= " STORED MOLECULES\n"
        if self.__bsm is None:
           log+= " --- no molecules added ---"
        else:
           for i,bsm in enumerate(self.__bsm):
               log+= " %5i %15s\n" % (i,bsm.get_name())
        return str(log)

    #   SET METHODS

    def set(self, pos, ind, nmol, bsm=None, supl=None, reord=None):
        """
 Set the initial molecular coordinates and the moltype index list. 
 Also set the BSM parameters if not done in set_bsm. 
 <pos>     -  an 2D array of atomic coordinates of dimension natoms, 3
 <ind>     -  a list of moltype ids in bsm list. Has length equal to the 
           -  total number of molecules
 <nmol>    - a list of number of atoms (integers) in each molecule
 <bsm>     - a list of a form: bsm = list(<Frag object>)
 <supl>    - a list of superimposition lists
 <reord>   - a list of reordering lists
"""
        self.__ind = numpy.array(ind,int)
        self.__nmol= numpy.array(nmol,int)
        if bsm is not None: self.__bsm = bsm
        if supl is not None: 
           self.__suplist = supl
           self.__suplist_c = supl[0]
        if reord is not None:
           self.__reordlist = reord
           self.__reordlist_c = reord[0]
        self._update(pos)
        self.__flds = None
        # store the information about central molecule
        if self.__eval_freq and bsm is not None:
           parc = self.__bsm[0].get()
           self.__pos_ref= parc['pos']
           self.__nmodes = parc['nmodes']
           self.__nata   = parc['natoms']
           self.__atno   = parc['atno']
           self.__atms   = parc['atms']
           self.__gijk   = parc['gijk']
           self.__redmss = parc['redmass']
           self.__freq   = parc['freq']
           self.__mode_sder_index = list(parc['mode']).index(self.__mode-1)
           self.__mol    = self.__bsm[0].get_mol()
           self.__avec   = None
           # calculate M matrix                                                  
           N = len(self.__mol.atoms)
           M = numpy.zeros((N*3,N*3), dtype=numpy.float64)
           for i in xrange(N):
               m = numpy.sqrt(self.__mol.atoms[i].mass()*self.AmuToElectronMass)
               M[3*i+0,3*i+0] = m
               M[3*i+1,3*i+1] = m
               M[3*i+2,3*i+2] = m
           self.__M = M
        # initialize vibrational forces and structural distortions
        if self.__eval_freq:
           self.__fi_el  = numpy.zeros(self.__nmodes, dtype=numpy.float64)   
           self.__fi_pol = numpy.zeros(self.__nmodes, dtype=numpy.float64)
           self.__fi_rep = numpy.zeros(self.__nmodes, dtype=numpy.float64)
           self.__dq_el  = numpy.zeros((self.__nata,3), dtype=numpy.float64)
           self.__dq_pol = numpy.zeros((self.__nata,3), dtype=numpy.float64)
           self.__dq_rep = numpy.zeros((self.__nata,3), dtype=numpy.float64)
        return    
    
    def set_cut(self,ccut,pcut,ecut):
        """Set the cutoff radii for Coulomb, polarization and exchange-repulsion forces"""
        self.__ccut = ccut
        self.__pcut = pcut
        self.__ecut = ecut
        return
    
    def set_bsm(self,bsm):
        """Set the BSM parameters. <bsm> is a list of a form:
  bsm = list(<Frag object>)"""
        self.__bsm = bsm
        return
        
    def set_pos(self,pos):
        """Update position array"""
        self._update(pos)
        return

    #   GET METHODS: SOLVATOCHROMISM

    def get_shifts(self):
        """Return SolEFP0 frequency shift data"""
        return self.__shift.copy()

    def get_fields(self):
        """Return the electric fields at all distributed polarizability sites as a tuple (Fields,r)"""
        return self.__flds

    def get_dipind(self):
        """Return the induced dipoles at all distributed polarizability sites as a tuple (Dip,r)"""
        return self.__dipind

    def get_avec(self):
        """Return the solvatochromic induced dipoles at all distributed polarizability sites as a tuple (Dip,r)"""
        return self.__avec

    def get_xyz_string(self):
        """Returns the fragment of XYZ file containing all EFP fragments within Coulomb cutoff. Implemented only for central mode!"""
        if not self.__eval_freq: 
           raise NotImplementedError, "Implemented only for central (frequency shift) mode!"
        else:
           log = ""
   
           frg = self.__bsm[0].copy()
           rms = frg.sup( self._reorder_xyz( self.__rc, self.__reordlist_c),
                                                          self.__suplist_c, )
           parc = frg.get() 
           pars = [parc,]

           N = len(self.__ntc)
           for i in xrange(N):
               im = self.__mtc[i]
               nm_prev = sum(self.__ntc[:i])
               nm_curr = sum(self.__ntc[:i+1])
               STR = self._reorder_xyz( self.__rcoordc[nm_prev:nm_curr]    , self.__reordlist[im] )  # reorder to BSM-FRG atom order
               frg = self.__bsm[im].copy()
               rms = frg.sup( STR , suplist= self.__suplist[im] )
               pars.append(frg.get())
           for par in pars:
               log += self._xyz(par['pos'], par['atno'], units='Angs')
           return log

    def get_forces(self, tot=False):
        """
Return solvation forces associated with solute's normal coordinates. 
Returns:
  tot=False (default)    tuple of three numpy.ndarray's of size (nmodes)
                         corresponds to electrostatic, repulsion and polarization 
                         forces, respectively
  tot=True               numpy.ndarray of size (nmodes)
                         which is the sum of all the three contributions to the 
                         total forces listed above

Notes: 
  o Dispersion and charge-gransfer forces are not implemented yet!
  o Values are in A.U. units of [Hartree/Bohr]
"""
        if   tot: return self.__fi_el+ self.__fi_rep+ self.__fi_pol
        else    : return self.__fi_el, self.__fi_rep, self.__fi_pol

    def get_kmatrix(self, tot=False):
        """
Return K-matrix from numerical calculation of derivatives. 
Returns:
  tot=False (default)    tuple of three numpy.ndarray's of size (nmodes x nmodes)
                         corresponds to electrostatic, repulsion and polarization 
                         K-matrix, respectively
  tot=True               numpy.ndarray of size (nmodes x nmodes)
                         which is the sum of all the three contributions to the 
                         total K-matrix listed above

Notes: 
  o Dispersion and charge-gransfer forces are not implemented yet!
  o Values are in A.U. units of [Hartree/Bohr**2]
"""
        if   tot: return self.__kij_el+ self.__kij_rep+ self.__kij_pol
        else    : return self.__kij_el, self.__kij_rep, self.__kij_pol

    def get_hess(self, tot=False):
        """
Return Hessian matrix from numerical calculation of derivatives. 
Returns:
  tot=False (default)    tuple of three numpy.ndarray's of size (nmodes x nmodes)
                         corresponds to electrostatic, repulsion and polarization 
                         Hessian, respectively
  tot=True               numpy.ndarray of size (nmodes x nmodes)
                         which is the sum of all the three contributions to the 
                         total Hessian listed above

Notes: 
  o Dispersion and charge-gransfer forces are not implemented yet!
  o Values are in A.U. units of [Hartree/Bohr**2]
"""
        if   tot: return self.__hess_el+ self.__hess_rep+ self.__hess_pol
        else    : return self.__hess_el, self.__hess_rep, self.__hess_pol

    def get_dq(self, theory=0, tot=False, cart=True):
        """
Return structural distortions in a given Tn-th level of theory=n. 
If tot=0, returns tuple of size 3 with electrostatics, repulsion and polarization
structural distortions numpy.ndarrays (nata x 3). If tot=1, returns total sum
of the structural distortion contributions. Returns everything in A.U. (Bohr)!
Uses suspicious eval_dq method! Should use _calc_dq instead"""
        dq_el = self._eval_dq(self.__fi_el , theory, cart=cart)
        dq_rep= self._eval_dq(self.__fi_rep, theory, cart=cart)
        dq_pol= self._eval_dq(self.__fi_pol, theory, cart=cart)
        if   tot: return dq_el+ dq_rep+ dq_pol
        else    : return dq_el, dq_rep, dq_pol


    def get_hessian(self, theory=0, mwc=True):
        """return the Hessian matrix in gas-phase normal coordinate space"""
        dq_el, dq_rep, dq_pol = self.get_dq(theory=theory, tot=0, cart=False)
        dq_tot = dq_el+ dq_rep+ dq_pol
        # calculate the Hessian in gas-phase normal coordinate space
        m = numpy.diag( self.__redmss * self.__freq**2.0 )
        s = numpy.tensordot(self.__gijk, dq_tot,(0,0))
        t = numpy.outer(numpy.sqrt(self.__redmss), numpy.sqrt(self.__redmss))
        h = (m + s)/ t
        print self.__freq * self.HartreePerHbarToCmRec
        # Calculate mass-weighted Hessian in Cartesian space (MWC)
        if mwc:
           t = numpy.sqrt(self.__redmss)[:,numpy.newaxis]
           l = self.__lvec.reshape(self.__nmodes, self.__nata * 3) / t
           l = numpy.dot(self.__M, l.transpose() )
           h = numpy.dot(l,numpy.dot(h,l.transpose()))
        return h


    def get_freq_approx(self,kij_el=False,diag=False):
        """Evaluate approximate frequencies without Hessian diagonalization, 
using appropriate level of SOL-X theory"""
        t = numpy.sqrt(numpy.outer(self.__redmss, self.__redmss))
        hess_tot_app = self.__hess_tot_app / t
        hess_tot     = self.__hess_tot     / t
        if diag:
           if kij_el: diag_hess,u = numpy.linalg.eig(hess_tot_app)
           else:      diag_hess,u = numpy.linalg.eig(hess_tot)
        else:
           if kij_el: diag_hess = numpy.diag(hess_tot_app)
           else:      diag_hess = numpy.diag(hess_tot)
        #diag_hess = numpy.where(diag_hess>0,diag_hess,0)
        #freq = where(diag_hess>0,sqrt(diag_hess),sqrt(-diag_hess) )
        freq = numpy.sqrt(diag_hess)
        if self.__cunit:
           freq *= self.HartreePerHbarToCmRec
        return freq

    def get_shift_approx(self,kij_el=False):
        """Evaluate approximate frequencies without Hessian diagonalization, 
using appropriate level of SOL-X theory"""
        #t = numpy.sqrt(numpy.outer(self.__redmss, self.__redmss))
        hess_tot_app = self.__hess_tot_app #/ t
        hess_tot     = self.__hess_tot     #/ t
        if kij_el: diag_hess = numpy.diag(hess_tot_app)
        else:      diag_hess = numpy.diag(hess_tot)
        diag_hess-= self.__redmss*self.__freq**2
        shift = diag_hess / (2.*self.__freq*self.__redmss)
        if self.__cunit:
           shift *= self.HartreePerHbarToCmRec
        return shift

    def get_slv(self, theory=0):
        """Calculate the frequencies, reduced masses and transformation matrix"""
        h = self.get_hessian(theory=theory, mwc=True)
        vib = libbbg.utilities.VIB(self.__mol, h, weight=False)
        vib.eval()
        freq, redmss, U = vib.get()
        return freq, redmss, U

    def get_shift_from_fi(self, fi, kij=None, iter=False, K4=None):
        """Calculate frequency shift from forces (in normal coordinate space)"""
        # approximate frequency shift based on (W^2 - w^2) \approx 2*w*(W-w)
        mid = self.__mode-1
        gijj = self.__gijk[:,mid,mid]                                   
        if kij is None:
           s = -numpy.sum(fi * gijj /(self.__redmss*self.__freq**2) ) / \
                (2.*self.__redmss[mid]*self.__freq[mid]) 
           if K4 is not None:
              q = fi/(self.__redmss*self.__freq**2)
              s+=  numpy.dot( numpy.dot(K4, q), q) [mid,mid]/ (4.*self.__redmss[mid]*self.__freq[mid]) 
        elif not iter:
           m = self.__redmss*self.__freq*self.__freq
           N = self.__nmodes
           G = self.__gijk.copy()/2.
           #G.fill(0.0)
           B = numpy.diag(m) + kij
           B1= numpy.linalg.inv(B)
           I = numpy.identity(self.__nmodes, numpy.float64)
           a = numpy.tensordot(B1,G,(1,1))
           b = numpy.tensordot(a,B1,(2,0))
           A = numpy.tensordot(b,fi,(2,0))

           D = numpy.dot( numpy.linalg.inv(I - A), B1)
           dq= numpy.dot(D,fi)

           T = numpy.tensordot(a, dq, (1,0))
           #print T
           #S = numpy.dot(self.__gijk, numpy.dot(D,fi))
           
           s =-numpy.sum(gijj*dq) / (2.*self.__redmss[mid]*self.__freq[mid])
           if K4 is not None:
              #q = fi/(self.__redmss*self.__freq**2)
              q = dq
              s+= numpy.dot( numpy.dot(K4, q), q) [mid,mid] / (4.*self.__redmss[mid]*self.__freq[mid]) 

        else:
           m = self.__redmss*self.__freq*self.__freq
           N = self.__nmodes
           G = self.__gijk.copy()/2.
           B = numpy.diag(m) + kij
           B1= numpy.linalg.inv(B)
           I = numpy.identity(self.__nmodes, numpy.float64)

           H = numpy.tensordot(B1,G,(1,0))
           fip = numpy.dot(B1, fi)
           dq_new =-fi/m
           dq_old = numpy.zeros(15,numpy.float64)

           while numpy.sum(numpy.abs(dq_old-dq_new)) > 0.001 :
              dq_old = dq_new.copy()
              dq_new =-fip - numpy.tensordot(H, numpy.outer(dq_old, dq_old), ((1,2),(0,1)) ) 
              print dq_new
              #A   = I + numpy.dot(B1, numpy.dot(G,dq_old))
              #A   = numpy.linalg.inv(A)
              #dq_new = -numpy.dot(A,fip)
           s = numpy.sum(gijj*dq_new) / (2.*self.__redmss[mid]*self.__freq[mid])

        #if self.__cunit:
        #   s*= self.HartreePerHbarToCmRec
        return s

    def get_shift_from_kii(self, kii):
        """Calculate frequency shift from diagonal hessian (in normal coordinate space)"""
        mid = self.__mode-1
        s =  kii[mid,mid] / \
             (2.*self.__redmss[mid]*self.__freq[mid])
        #if self.__cunit:
        #   s*= self.HartreePerHbarToCmRec
        return s

    def get_shift_from_fi_kii(self, fi, kii):
        """"""
        mid = self.__mode-1
        gijj = self.__gijk[:,mid,mid]
        fii = self.__redmss[mid] * self.__freq[mid]**2
        s = ( fii + kii[mid,mid] - numpy.sum(fi * gijj /(self.__redmss*self.__freq**2) ) )/ self.__redmss[mid]
        s = numpy.sqrt(s) - self.__freq[mid]
        if self.__cunit:
           s*= self.HartreePerHbarToCmRec
        return s


    # STRUCTURAL FUNCTIONALITIES
 
    def get_pos_calc(self, theory=0, units='bohr'):
        """Return the structure of the fragment calculated from SolX theory"""
        c = 1.0000
        if units.lower().startswith('ang'): c = self.BohrToAngstrom
        pos = self.__pos_c + self.get_dq(theory,tot=1)
        return pos * c

    def get_pos_c(self, units='bohr'):
        """Return solute superimposed structure (relevant if central molecule is ON). 
This structure has identical internal coordinates as fragment BSM gas-phase structure"""
        c = 1.0000
        if units.lower().startswith('ang'): c = self.BohrToAngstrom
        return self.__pos_c.copy() * c
  
    def get_pos_ref(self, units='bohr'):
        """Return solute reference structure (relevant if central molecule is ON). 
This structure has identical internal coordinates as fragment BSM gas-phase structure.
Note that this structure is the one BEFORE performing superimposition (while get_pos_c is AFTER), 
and if there is no copying of __bsm[0] object the get_pos_ref orientation and position will be different than the original
BSM structure from SLV parameter file library."""
        c = 1.0000
        if units.lower().startswith('ang'): c = self.BohrToAngstrom
        return self.__pos_ref.copy() * c

    def get_rot_transl(self):
        """Return the rotation matrix and translation vector defining the translation:
   pos_ref = numpy.dot(pos_c,rot) + transl
Note that the RMS of superimposition is always exactly ZERO in this case.
See also: get_pos_c, get_pos_ref"""
        s = libbbg.utilities.SVDSuperimposer()
        s.set(self.__pos_ref, self.__pos_c)
        s.run()
        rms         = s.get_rms()
        rot, transl = s.get_rotran()
        return rot, transl

    def get_center_c(self,r1,r2=None):
        """Return a center [x,y,z]
Usage:
 get_center_c(r1)    - [r1+1]-th atomic centre in __pos_c
 get_center_c(r1,r2) - midbond point between [r1+1]-th and [r2+1]-th centres in __pos_c.
See also: get_pos_c, get_pos_ref"""
        if r2 is None: return self.__pos_c[r1]
        else: return 0.500*(self.__pos_c[r1] + self.__pos_c[r2])

    def get_center_ref(self,r1,r2=None):
        """Return a center [x,y,z]
Usage:
 get_center_c(r1)    - [r1+1]-th atomic centre in __pos_ref
 get_center_c(r1,r2) - midbond point between [r1+1]-th and [r2+1]-th centres in __pos_ref.
See also: get_pos_c, get_pos_ref"""
        if r2 is None: return self.__pos_ref[r1]
        else: return 0.500*(self.__pos_ref[r1] + self.__pos_ref[r2])

    def get_rms(self):
        """
Return RMS of superimposition of central molecule 
with its parameters (relevant for central molecule mode)"""
        return self.__rms_central

    def get_rms_sol(self):
        """
Return maximal RMS of superimposition of (solvent) molecules with its parameters.
If central molecule mode is ON, rms_central is not included. Otherwise, all superimpositions
are counted."""
        return self.__rms_solvent_max

    def get_rms_ave(self):
        """
Return average RMS over all molecules (including central molecule if any)
"""
        return self.__rms_ave

    def get_rms_max(self):
        """Return maximal RMS"""
        rms_c = self.__rms_central
        rms_s = self.__rms_solvent_max
        if self.__rms_central is not None:
           return max(self.__rms_central, self.__rms_solvent_max)
        else:
           return self.__rms_solvent_max

    def get_pos_sol(self, units='bohr'):
        """Return solvent superimposed structure (relevant if central molecule is ON).
Now, only for exchange-repulsion layer"""
        c = 1.0000
        if units.lower().startswith('ang'): c = self.BohrToAngstrom
        return self.__solvent_exrep.copy() * c

    def get_rc(self):
        """Return solute target structure (relevant if central molecule mode is ON)"""
        return self.__rc.copy()

    def get_u(self):
        """Returns unitary transformation matrix diagonalizing the Hessian"""
        return self.__u


    
    # P R O T E C T E D
    
    def _create(self):
        """namespace of objects"""
        self.__debug             = None
        self.__mode              = None
        self.__bsm               = None
        self.__eval_freq         = False
        self.__eval_nlo          = None
        self.__shift             = dict() # numpy.zeros(9,dtype=numpy.float64)
        self.__rms_central       = None
        self.__rms_ave           = None
        self.__suplist           = None
        self.__reordlist         = None
        self.__rms_solvent_max   = None
        self.__suplist_c         = None
        self.__reordlist_c       = None
        # current position of solute superimposed
        self.__pos_c             = None
        # current position of solvent superimposed
        self.__solvent_elect     = None
        self.__solvent_polar     = None
        self.__solvent_exrep     = None
        return
  
    def _update(self,pos):
        """update the neighbour/in-sphere lists based on actual <pos> coordinate array"""
        if not self.__pairwise_all and len(self.__nmol)-1:
           nm = len(self.__nmol)-1
           nac= self.__nmol[0]
           natoms = len(pos)
           rc = pos[:nac ]
           rm = pos[ nac:].reshape((natoms-nac)*3)
           #
           ic  = numpy.zeros(natoms-nac, dtype=bool)
           ip  = numpy.zeros(natoms-nac, dtype=bool)
           ie  = numpy.zeros(natoms-nac, dtype=bool)
           icm = numpy.zeros(nm        , dtype=bool)
           ipm = numpy.zeros(nm        , dtype=bool)
           iem = numpy.zeros(nm        , dtype=bool)
           #
           nccut, npcut, necut, mccut, mpcut, mecut = \
                                        solvshift.solpol2.mollst(rc,rm,ic,ip,ie,icm,ipm,iem,
                                                      self.__nmol[1:],
                                                      self.__ccut,
                                                      self.__pcut,
                                                      self.__ecut)
                                                      #
           # molecule coordinate lists
           rcoordc = pos[nac:][numpy.array(nccut, dtype=bool)]
           rcoordp = pos[nac:][numpy.array(npcut, dtype=bool)]
           rcoorde = pos[nac:][numpy.array(necut, dtype=bool)]
           # molecule type lists
           mtc= self.__ind[1:][numpy.array(mccut, dtype=bool)]
           mtp= self.__ind[1:][numpy.array(mpcut, dtype=bool)]
           mte= self.__ind[1:][numpy.array(mecut, dtype=bool)]
           # moleculear atom number lists
           ntc= self.__nmol[1:][numpy.array(mccut, dtype=bool)]
           ntp= self.__nmol[1:][numpy.array(mpcut, dtype=bool)]
           nte= self.__nmol[1:][numpy.array(mecut, dtype=bool)]
           #
           self.__rc = rc
           self.__rcoordc = rcoordc
           self.__rcoordp = rcoordp
           self.__rcoorde = rcoorde
           self.__mtc = mtc
           self.__mtp = mtp
           self.__mte = mte
           self.__ntc = ntc
           self.__ntp = ntp
           self.__nte = nte
        else:
           self.__rcoordc = pos
           self.__rc = pos[:self.__nmol[0]]  # the case of no EFP molecules in the environment
        return

    # PROPERTY EVALUATORS

    def _eval(self, lwrite, num, step, theory, remove_clashes, rcl_algorithm):
        if lwrite>1: self.__debug = open('__debug_file__','w')
        if self.__pairwise_all:  self._eval_mode_global (lwrite, num, step, theory)
        else:                    self._eval_mode_central(lwrite, num, step, theory, remove_clashes, rcl_algorithm) 
        return

    def _eval_dma(self, dma, lwrite):
        """Cutoffs not supported yet for ENVIRON layer - just one object is considered"""
        if lwrite: print " WARNING! No superimposition for ENVIRON yet! Please supply DMA already in appropriate position!!!"
        shift_total = 0.0
        shift_mea   = 0.0; shift_ea = 0.0; shift_corr_mea = 0.0; shift_corr_ea = 0.0

        # --- extract the multipoles for the environment
        N_dma_environ = len(dma)
        if lwrite: print " ENVIRONMENT: %10d distributed centers" % N_dma_environ

        # DMA objects support
        if   isinstance(dma, libbbg.dma.DMA):
             assert not dma.is_traceless, " The DMA is in traceless form already! Provide primitive form instead!"
             dma_copy = dma.copy(); dma_copy.trac()
             xyz_environment = dma.get_pos()
             chg_environment = dma.get_charges()
             dip_environment = dma.get_dipoles()
             qad_environment = dma_copy.get_quadrupoles()
             oct_environment = dma_copy.get_octupoles(); del dma_copy
        # NumPy NDArray support (format: [[X,Y,Z,q],....] all in AU!)
        elif isinstance(dma, numpy.ndarray):
             assert dma.shape[1]==4, " Incorrect shape of ndarray! %s, The last dimension needs to be 4!" % dma.shape
             xyz_environment = dma[:,:3]
             chg_environment = dma[:,3]
             dip_environment = numpy.zeros((N_dma_environ, 3),numpy.float64)
             qad_environment = numpy.zeros((N_dma_environ, 6),numpy.float64)
             oct_environment = numpy.zeros((N_dma_environ,10),numpy.float64)
        # handle different data types
        else: 
           print " This data structure is not supported!\n Only libbbg.dma.DMA or numpy.ndarray objects are handled.\n"; exit()

        # central molecule
        frg = self.__bsm[0].copy(); parc= frg.copy().get()
        ndmac= parc['ndma']

        # superimposition of central molecule
        self.__rms_central = frg.sup( self._reorder_xyz( self.__rc, self.__reordlist_c),
                                      self.__suplist_c, dxyz=None)

        self.__pos_c = frg.get_pos()

        # properties of central molecule
        parc= frg.get()
        lvec = parc['lvec'].ravel()
        gijj = self.__gijk[:,self.__mode-1,self.__mode-1]
        freq   = self.__freq
        redmss = self.__redmss
        nmodes = self.__nmodes
        qadc, octc = frg.get_traceless()
        #
        freqc = freq[self.__mode-1]
        redmssc=redmss[self.__mode-1]

        # parse analytical derivatives of DMTP
        chgc1 = parc['dmac1'].ravel()
        dipc1 = parc['dmad1'].ravel()
        qadc1, octc1 = frg.get_traceless_1()
        chgc2 = parc['dmac2'][self.__mode_sder_index].ravel()
        dipc2 = parc['dmad2'][self.__mode_sder_index].ravel()
        qadc2, octc2 = frg.get_traceless_2(self.__mode)
        qadc2 = qadc2.ravel(); octc2 = octc2.ravel()

        ndma = [ parc['ndma'] ]; ndma.append(N_dma_environ)
        ndmas= sum(ndma)
        rdma = numpy.concatenate  ([ parc['rdma'], xyz_environment ]).reshape(ndmas*3)
        chg  = numpy.concatenate  ([ parc['dmac'], chg_environment ]).reshape(ndmas)
        dip  = numpy.concatenate  ([ parc['dmad'], dip_environment ]).reshape(ndmas*3)
        qad  = numpy.concatenate  ([       qadc  , qad_environment ]).reshape(ndmas*6)
        oct  = numpy.concatenate  ([       octc  , oct_environment ]).reshape(ndmas*10)


        # --- compute the frequency shifts
        # mechanical anharmonicity
        mea,a,b,c,d,e,fi = libbbg.qm.clemtp.sdmtpm(rdma,ndma,chg,dip,qad,oct,
                           chgc1,dipc1,qadc1,octc1,redmss,freq,gijj,ndmac,self.__mode,lwrite=False)
        self.__fi_dma_el = fi
        shift_mea= mea

        # electronic anharmonicity
        if self.__eval_elect_ea:
           ea ,a,b,c,d,e = libbbg.qm.clemtp.sdmtpe(rdma,ndma,chg,dip,qad,oct,
                           chgc2,dipc2,qadc2,octc2,redmss,freq,self.__mode,lwrite=False)
           shift_ea  = ea

        # potential-derivative correction terms
        if self.__eval_corr:
           corr,rf2,rf3,rf4,rk2,rk3,rk4,corr_b,corr_c,corr_d = \
                        libbbg.qm.clemtp.dmtcor(rdma,ndma,chg,dip,qad,oct,
                        chgc1,dipc1,qadc1,octc1,redmss,freq,gijj,lvec,ndmac,self.__mode,lwrite=False)
           shift_corr_mea = rf2#+rf3#+rf4
           shift_corr_ea  = rk2#+rk3#+rk4


        if self.__cunit:
           shift_mea      *= self.HartreePerHbarToCmRec
           shift_ea       *= self.HartreePerHbarToCmRec
           shift_corr_mea *= self.HartreePerHbarToCmRec
           shift_corr_ea  *= self.HartreePerHbarToCmRec

        shift_total = shift_mea + shift_ea + shift_corr_mea + shift_corr_ea

        # store the values
        self.__shift['ele_env_mea' ] = shift_mea      ; self.__shift['ele_env_ea' ] = shift_ea
        self.__shift['cor_env_mea' ] = shift_corr_mea ; self.__shift['cor_env_ea' ] = shift_corr_ea
        self.__shift['tot_env'     ] = shift_total
        return
    
    def _eval_mode_global(self, lwrite, num, step, theory):
        merror = ' Global mode is not implemented yet!'
        raise NotImplementedError, merror

        # compute exchange-repulsion energy
        for parA in PAR:
            for parB in PAR:
                e_exrep += self._pair_rep(parA,parB)

        if self.__cunit:
           e_exrep *= self.HartreeToKcalPerMole

        self.__energy_exrep = e_exrep
        return
   
    def _eval_mode_central(self, lwrite, num, step, theory, remove_clashes, rcl_algorithm, dxyz=None):
        """Includes only EFP layers (DMA environment is not accounted for!)"""
        assert step==0.006, ' Step %10.5f Angstrom is not supported yet!' % step
    
        # initialize frequency shifts
        shift_total  = 0.0 
        shift_ele_mea= 0.0; shift_pol_mea = 0.0; shift_rep_mea = 0.0; shift_dis_mea = 0.0
        shift_ele_ea = 0.0; shift_pol_ea  = 0.0; shift_rep_ea  = 0.0; shift_dis_ea  = 0.0
        shift_ele_corr_mea = 0.0 ; shift_ele_corr_ea = 0.0
        shift_dis_mea_iso = 0.0
        shift_pol_add_mea = 0.0  ;  shift_pol_add_ea = 0.0

        # central molecule
        frg = self.__bsm[0].copy(); parc= frg.copy().get()
        if num: self.__lvec = parc['lvec']
        ndmac= parc['ndma']  
    
        # superimposition of central molecule
        self.__rms_central = frg.sup( self._reorder_xyz( self.__rc, self.__reordlist_c), 
                                                          self.__suplist_c, dxyz=dxyz)
        self.__pos_c = frg.get_pos()
    
        if lwrite:
           print " %s" % frg.get()['name']
           print " Central molecule rms: %10.5f" % self.__rms_central
           xx = frg.xyz(units='angs'); print xx
        if lwrite>1: self.__debug.write(xx)
    
        # properties of central molecule
        if not num:
           parc= frg.get()
           self.__lvec = parc['lvec']
           gijj = self.__gijk[:,self.__mode-1,self.__mode-1]
        else:
           gijk   = self.__gijk
        freq   = self.__freq
        redmss = self.__redmss
        lvec   = self.__lvec.ravel()
        nmodes = self.__nmodes
        qadc, octc = frg.get_traceless()
        #
        freqc = freq[self.__mode-1]
        redmssc=redmss[self.__mode-1]
    
        N = len(self.__ntc)
        nmols = N+1
        
    
        # determine the reverse transformation tensors and read FD fragments
        if num:
           rot, transl, rms = frg.get_rotranrms()
           rot_inv = rot.transpose()
           transl_inv = -transl.copy()
    
           _dir_ = frg.get_dir(); N_FD = 1 + 9*self.__nata**2 + 3*self.__nata
           files = glob.glob(_dir_+'/num_0.006/*.frg')
           files.sort()
           
           PAR_FD = numpy.zeros( N_FD, dtype=object )
           QO_FD  = numpy.zeros( N_FD, dtype=object )
           PAR_SOL = list()
           QO_SOL  = list()
           for i in xrange(N_FD):
               frg_i     = solvshift.slvpar.Frag(files[i])
               PAR_FD[i] = frg_i.get()
               QO_FD [i] = frg_i.get_traceless()
        # parse analytical derivatives of DMTP
        else:
          PAR = list(); QO = list()
          PAR.append(parc)
          QO.append( (qadc, octc) )
    
          chgc1 = parc['dmac1'].ravel()
          dipc1 = parc['dmad1'].ravel()
          qadc1, octc1 = frg.get_traceless_1()
          chgc2 = parc['dmac2'][self.__mode_sder_index].ravel()
          dipc2 = parc['dmad2'][self.__mode_sder_index].ravel()
          qadc2, octc2 = frg.get_traceless_2(self.__mode)
          qadc2 = qadc2.ravel(); octc2 = octc2.ravel()
           
        # loop over other molecules
        rms_max = 0.0 
        rms_ave = self.__rms_central
        if lwrite: print " ELECT LAYER: %10d molecules" % N
        for i in xrange(N):
            im = self.__mtc[i]
            nm_prev = sum(self.__ntc[:i])
            nm_curr = sum(self.__ntc[:i+1])
            STR = self._reorder_xyz( self.__rcoordc[nm_prev:nm_curr]    , self.__reordlist[im] )  # reorder to BSM-FRG atom order
            if num: STR = numpy.dot(STR+transl_inv, rot_inv)
            frg = self.__bsm[im].copy()
            rms = frg.sup( STR , suplist= self.__suplist[im] ); rms_ave += rms
            if lwrite: print " %s" % frg.get()['name']
            if lwrite: print " rms C: ",rms
            if lwrite:
               xx = frg.xyz(units='angs'); print xx
            if lwrite>1: self.__debug.write(xx)
            if rms > rms_max: rms_max = rms
            par = frg.get()
            qad, oct = frg.get_traceless()
            if num:
               PAR_SOL.append( par )
               QO_SOL.append( (qad, oct) )
            else:
               PAR.append( par )
               #qad, oct = solvshift.efprot.tracls( par['dmaq'], par['dmao'] )
               QO.append( (qad, oct) )
            self.__rms_solvent_max = rms_max
        self.__rms_ave = rms_ave / numpy.float64(N+1) 

        # ELECTROSTATICS
        if self.__eval_elect:
           # COULOMB
           if num:
              func_el  = numpy.zeros(N_FD, numpy.float64)
              ndma_sol = [ x['ndma'] for x in PAR_SOL ]
              rdma_sol = numpy.concatenate  ([ x['rdma'] for x in PAR_SOL ])
              chg_sol  = numpy.concatenate  ([ x['dmac'] for x in PAR_SOL ])
              dip_sol  = numpy.concatenate  ([ x['dmad'] for x in PAR_SOL ])
              qad_sol  = numpy.concatenate  ([ QO_SOL[x][0]  for x in range(N)   ])
              oct_sol  = numpy.concatenate  ([ QO_SOL[x][1]  for x in range(N)   ])
              for i in xrange(N_FD):
                  p    = PAR_FD[i]
                  q,o  = QO_FD[i]
                  ndma = [p['ndma'],] + ndma_sol
                  rdma = numpy.concatenate ( (p['rdma'], rdma_sol) ).ravel()
                  chg  = numpy.concatenate ( (p['dmac'],  chg_sol) ).ravel()
                  dip  = numpy.concatenate ( (p['dmad'],  dip_sol) ).ravel()
                  qad  = numpy.concatenate ( (q        ,  qad_sol) ).ravel()
                  oct  = numpy.concatenate ( (o        ,  oct_sol) ).ravel()
                  func_el[i],A,B,C,D,E = libbbg.qm.clemtp.edmtpc(rdma,chg,dip,qad,oct,ndma,lwrite=False) 
              # differentiate and obtain forces f and hessian K
              fx_1 = func_el[0:self.__nata*6+1]
              fx_2 = func_el[  self.__nata*6+1: ]
              fd_calc = libbbg.utilities.diff(func=None, step=0.006*self.AngstromToBohr, 
                                                   DIM=self.__nata*3, scheme='3pt')
              fi_cart, kij_cart = fd_calc.eval( (fx_1, fx_2), symm=True)
              # transform the derivatives to normal coordinate space
              l = self.__lvec.reshape(self.__nmodes, self.__nata*3)
              fi_mode = numpy.dot(l, fi_cart)
              kij_mode= numpy.dot( numpy.dot(l, kij_cart), l.transpose() )
              # store forces and K-matrices
              self.__fi_el  = fi_mode
              self.__kij_el = kij_mode
              # diagonalize the hessian
              hess, freq, redmass, U = self._get_slv(type='e', theory=theory)
              libbbg.utilities.PRINT(freq[::-1])
              self.__hess_el = hess
              # compute approximated frequencies at SolEFP-0 theory
              shift_ele_mea = self.get_shift_from_fi (fi_mode, iter=False, K4=None, kij=None)
              shift_ele_ea  = self.get_shift_from_kii(kij_mode)
              shift_ele_corr_mea = 0.0
              shift_ele_corr_ea  = 0.0
              #
              del PAR_SOL, QO_SOL
           else:
              ndma = [ x['ndma'] for x in PAR ]
              ndmas= sum(ndma)
              rdma = numpy.concatenate  ([ x['rdma'] for x in PAR ]).reshape(ndmas*3)
              chg  = numpy.concatenate  ([ x['dmac'] for x in PAR ]).reshape(ndmas)
              dip  = numpy.concatenate  ([ x['dmad'] for x in PAR ]).reshape(ndmas*3)
              qad  = numpy.concatenate  ([ QO[x][0]  for x in range(N+1)   ]).reshape(ndmas*6)
              oct  = numpy.concatenate  ([ QO[x][1]  for x in range(N+1)   ]).reshape(ndmas*10)
              # mechanical anharmonicity
              mea,a,b,c,d,e,fi = libbbg.qm.clemtp.sdmtpm(rdma,ndma,chg,dip,qad,oct,
                                 chgc1,dipc1,qadc1,octc1,redmss,freq,gijj,ndmac,self.__mode,lwrite=False)
              self.__fi_el = fi
              shift_ele_mea= mea
              # electronic anharmonicity
              if self.__eval_elect_ea:
                 ea ,a,b,c,d,e = libbbg.qm.clemtp.sdmtpe(rdma,ndma,chg,dip,qad,oct,
                                 chgc2,dipc2,qadc2,octc2,redmss,freq,self.__mode,lwrite=False)
                 shift_ele_ea  = ea
              # potential-derivative correction terms
              if self.__eval_corr:
                 corr,rf2,rf3,rf4,rk2,rk3,rk4,corr_b,corr_c,corr_d = \
                              libbbg.qm.clemtp.dmtcor(rdma,ndma,chg,dip,qad,oct,
                              chgc1,dipc1,qadc1,octc1,redmss,freq,gijj,lvec,ndmac,self.__mode,lwrite=False)
                 shift_ele_corr_mea = rf2#+rf3#+rf4
                 shift_ele_corr_ea  = rk2#+rk3#+rk4
              #
              del PAR, QO
    
           # POLARIZATION 
           if self.__eval_pol:
              npolc = parc['npol']
              N = len(self.__ntp)
              if num:
                 PAR_SOL = list(); QO_SOL  = list()
              else:
                 PAR = list(); QO = list()
                 PAR.append( parc ); QO.append( (qadc, octc) )
              if lwrite: print " POL   LAYER: %10d molecules" % N

              # loop over solvent molecules in POL layer
              for i in xrange(N):
                  im = self.__mtp[i]
                  nm_prev = sum(self.__ntp[:i])
                  nm_curr = sum(self.__ntp[:i+1])
                  STR = self._reorder_xyz( self.__rcoordp[nm_prev:nm_curr]      , self.__reordlist[im])
                  if num: STR = numpy.dot(STR+transl_inv, rot_inv)
                  frg = self.__bsm[im].copy()
                  rms = frg.sup( STR, suplist= self.__suplist[im] )
                  if 0:#lwrite:
                     xx = frg.xyz(units='angs'); print xx
                  par = frg.get()
                  if num:
                     PAR_SOL.append( par ); QO_SOL.append( frg.get_traceless() )                     
                  else:
                     PAR.append( par )
                     QO.append( frg.get_traceless() )

              # eliminate clashing molecules
              if remove_clashes:
                 PAR, QO, PAR_ADD, QO_ADD = self._remove_clashes(PAR, QO, lwrite, rcl_algorithm)
                 N_ADD = len(PAR_ADD) - 1   # these count only solvents (minus solute)
                 N     = len(PAR)     - 1
                 if lwrite: 
                    print "            : %10d molecules after clash removal" % N    
                    print "            : %10d molecules removed            " % N_ADD    
                 if lwrite:
                    print
                    print " Coordinates of POL LAYER after removal of clashes"
                    for par in PAR:
                        print self._xyz(par['pos'], par['atno'], units='Angs'),
                    print
                

              # prepare the data structures for passing to Fortran soubroutines
              if num:
                 func_pol = numpy.zeros(N_FD, numpy.float64)
                 ndma_sol =                     [ x['ndma'] for x in PAR_SOL ]
                 npol_sol =                     [ x['npol'] for x in PAR_SOL ]
                 rdma_sol = numpy.concatenate  ([ x['rdma'] for x in PAR_SOL ])
                 rpol_sol = numpy.concatenate  ([ x['rpol'] for x in PAR_SOL ])
                 pol_sol  = numpy.concatenate  ([ x['dpol'] for x in PAR_SOL ])
                 chg_sol  = numpy.concatenate  ([ x['dmac'] for x in PAR_SOL ])
                 dip_sol  = numpy.concatenate  ([ x['dmad'] for x in PAR_SOL ])
                 qad_sol  = numpy.concatenate  ([ QO_SOL[x][0]  for x in range(N)   ])
                 oct_sol  = numpy.concatenate  ([ QO_SOL[x][1]  for x in range(N)   ])
                 npols= sum(npol_sol) + PAR_FD[0]['npol']
                 DIM  = npols*3
                 for i in xrange(N_FD):
                     dmat = numpy.zeros((DIM,DIM), numpy.float64)
                     flds = numpy.zeros( DIM, numpy.float64)
                     dipind=numpy.zeros( DIM, numpy.float64)
                     p    = PAR_FD[i]
                     q,o  = QO_FD[i]
                     ndma = [p['ndma'],] + ndma_sol
                     npol = [p['npol'],] + npol_sol
                     rdma = numpy.concatenate ( (p['rdma'], rdma_sol) ).ravel()
                     chg  = numpy.concatenate ( (p['dmac'],  chg_sol) ).ravel()
                     dip  = numpy.concatenate ( (p['dmad'],  dip_sol) ).ravel()
                     qad  = numpy.concatenate ( (q        ,  qad_sol) ).ravel() 
                     oct  = numpy.concatenate ( (o        ,  oct_sol) ).ravel()
                     rpol = numpy.concatenate ( (p['rpol'], rpol_sol) ).ravel()
                     pol  = numpy.concatenate ( (p['dpol'],  pol_sol) ).ravel()
                     func_pol[i] = solvshift.solpol.solpol(rdma,chg,dip,qad,oct,
                                   rpol,pol,dmat,flds,dipind,ndma,npol,lwrite=False)
                 # differentiate and obtain forces f and hessian K                            
                 fx_1 = func_pol[0:self.__nata*6+1]                                             
                 fx_2 = func_pol[  self.__nata*6+1: ]                                           
                 fd_calc = libbbg.utilities.diff(func=None, step=0.006*self.AngstromToBohr,     
                                                 DIM=self.__nata*3, scheme='3pt')               
                 fi_cart, kij_cart = fd_calc.eval( (fx_1, fx_2), symm=True)                     
                 # transform the derivatives to normal coordinate space                         
                 l = self.__lvec.reshape(self.__nmodes, self.__nata*3)                          
                 fi_mode = numpy.dot(l, fi_cart)                                                
                 kij_mode= numpy.dot( numpy.dot(l, kij_cart), l.transpose() )                   
                 # store forces and K-matrices                                                  
                 self.__fi_pol = fi_mode                                                        
                 self.__kij_pol = kij_mode                                                      
                 # diagonalize the hessian                                                      
                 hess, freq, redmass, U = self._get_slv(type='p', theory=theory)                
                 libbbg.utilities.PRINT(freq[::-1])                                             
                 self.__hess_pol = hess                                                         
                 # compute approximated frequencies at SolEFP-0 theory
                 shift_pol_mea = self.get_shift_from_fi (fi_mode, iter=False, K4=None, kij=None)
                 shift_pol_ea  = self.get_shift_from_kii(kij_mode)
                 #                                                                              
                 del PAR_SOL, QO_SOL
              else:
                 ndma = [ x['ndma'] for x in PAR ]                                                 
                 npol = [ x['npol'] for x in PAR ]
                 ndmas= sum(ndma)
                 npols= sum(npol)

                 rdma = numpy.concatenate  ([ x['rdma'] for x in PAR ]).reshape(ndmas*3)
                 chg  = numpy.concatenate  ([ x['dmac'] for x in PAR ]).reshape(ndmas)
                 dip  = numpy.concatenate  ([ x['dmad'] for x in PAR ]).reshape(ndmas*3)
                 qad  = numpy.concatenate  ([ QO[x][0]  for x in range(N+1)   ]).reshape(ndmas*6)
                 oct  = numpy.concatenate  ([ QO[x][1]  for x in range(N+1)   ]).reshape(ndmas*10)
                 rpol = numpy.concatenate  ([ x['rpol'] for x in PAR ]).reshape(npols*3)
                 pol  = numpy.concatenate  ([ x['dpol'] for x in PAR ])#.reshape(npols*9)
                 polinv = numpy.concatenate([ numpy.linalg.inv(x) for x in pol ]).reshape(npols*9)
                 pol  = pol.ravel()
                 rpol1= parc['lmoc1'].reshape(nmodes*npolc*3)
                 pol1 = parc['dpol1'].reshape(nmodes*npolc*9)
                 DIM  = npols*3
                 dmat = numpy.zeros((DIM,DIM), numpy.float64)
                 dimat= numpy.zeros((DIM,DIM), numpy.float64)
                 mat1 = numpy.zeros((DIM,DIM), numpy.float64)
                 flds = numpy.zeros( DIM, numpy.float64)
                 vec1  =numpy.zeros( DIM, numpy.float64)
                 vec2  =numpy.zeros( DIM, numpy.float64)
                 fivec =numpy.zeros( DIM, numpy.float64)
                 mivec =numpy.zeros( DIM, numpy.float64)
                 #dipind= flds.copy(); avec  = flds.copy()
                 epol, spol, fi, avec, dipind = solvshift.solpol2.sftpli(rdma, chg, dip, qad, oct,
                       chgc1, dipc1, qadc1, octc1,rpol, polinv, mivec, dmat, flds, dimat, fivec,
                       vec1, vec2, mat1, redmss, freq, gijj, rpol1, pol1, lvec, ndma, npol, 
                       self.__mode, ndmac, npolc, lwrite=False)
                 #spol = 0.0
                 self.__fi_pol = fi
                 shift_pol_mea = spol
                 # store field values at all polarizable centers
                 self.__flds   =(flds  .reshape(npols,3), rpol.reshape(npols,3))
                 # store induced dipole moments at all polarizable centers
                 self.__dipind =(dipind.reshape(npols,3), rpol.reshape(npols,3))
                 # store solvatochromic induced dipole moments at all polarizable centers
                 self.__avec   =(avec  .reshape(npols,3), rpol.reshape(npols,3))
                 # 
                 if remove_clashes and N_ADD:
                    npols_add=sum([ x['npol'] for x in PAR_ADD ])
                    #
                    flds_add   = numpy.zeros(0, numpy.float64)
                    dipind_add = numpy.zeros(0, numpy.float64)
                    avec_add   = numpy.zeros(0, numpy.float64)
                    #
                    #self.__fi_pol = numpy.zeros(self.__nmodes, numpy.float64)
                    #
                    for imol in range(N_ADD):
                        if lwrite: " Polarization due to molecule %d" % (imol+1)
                        PARLOC = [PAR_ADD[0], PAR_ADD[imol+1]]
                        QOLOC  = [QO_ADD [0],  QO_ADD[imol+1]]
                        #
                        ndma = [ x['ndma'] for x in PARLOC ]
                        npol = [ x['npol'] for x in PARLOC ]
                        ndmas= sum(ndma)
                        npols= sum(npol)
                        #
                        rdma = numpy.concatenate  ([ x['rdma'] for x in PARLOC ]).reshape(ndmas*3)           
                        chg  = numpy.concatenate  ([ x['dmac'] for x in PARLOC ]).reshape(ndmas)
                        dip  = numpy.concatenate  ([ x['dmad'] for x in PARLOC ]).reshape(ndmas*3)
                        qad  = numpy.concatenate  ([ x[ 0    ] for x in QOLOC  ]).reshape(ndmas*6)
                        oct  = numpy.concatenate  ([ x[ 1    ] for x in QOLOC  ]).reshape(ndmas*10)
                        rpol = numpy.concatenate  ([ x['rpol'] for x in PARLOC ]).reshape(npols*3)
                        pol  = numpy.concatenate  ([ x['dpol'] for x in PARLOC ])#.reshape(npols*9)
                        polinv = numpy.concatenate([ numpy.linalg.inv(x) for x in pol ]).reshape(npols*9)
                        pol  = pol.ravel()
                        rpol1= parc['lmoc1'].reshape(nmodes*npolc*3)
                        pol1 = parc['dpol1'].reshape(nmodes*npolc*9)
                        DIM  = npols*3
                        dmat = numpy.zeros((DIM,DIM), numpy.float64)
                        dimat= numpy.zeros((DIM,DIM), numpy.float64)
                        mat1 = numpy.zeros((DIM,DIM), numpy.float64)
                        flds = numpy.zeros( DIM, numpy.float64)
                        vec1  =numpy.zeros( DIM, numpy.float64)
                        vec2  =numpy.zeros( DIM, numpy.float64)
                        fivec =numpy.zeros( DIM, numpy.float64)
                        mivec =numpy.zeros( DIM, numpy.float64)
                        #
                        epol, spol, fi, avec, dipind = solvshift.solpol2.sftpli(rdma, chg, dip, qad, oct,
                              chgc1, dipc1, qadc1, octc1,rpol, polinv, mivec, dmat, flds, dimat, fivec,
                              vec1, vec2, mat1, redmss, freq, gijj, rpol1, pol1, lvec, ndma, npol, 
                              self.__mode, ndmac, npolc, lwrite=False)
                        #
                        shift_pol_add_mea += spol
                        self.__fi_pol += fi
                        # ammend the fields on solute due to additive POL layer
                        #print self.__flds  [0][:npolc].shape, flds  .reshape(npols,3).shape; exit()
                        self.__flds  [0][:npolc]+= flds  .reshape(npols,3)[:npolc]
                        self.__dipind[0][:npolc]+= dipind.reshape(npols,3)[:npolc]
                        self.__avec  [0][:npolc]+= avec  .reshape(npols,3)[:npolc]
                        # fields due to additive POL layer
                        if not imol:  # initiate the vector fields
                           flds_add   = numpy.concatenate( (flds_add  , flds  ) )
                           dipind_add = numpy.concatenate( (dipind_add, dipind) )
                           avec_add   = numpy.concatenate( (avec_add  , avec  ) )
                        else:         # append new solvent fields and add contributions to solute fields
                           flds_add   = numpy.concatenate( (flds_add  , flds  [npolc*3:]) )
                           dipind_add = numpy.concatenate( (dipind_add, dipind[npolc*3:]) )
                           avec_add   = numpy.concatenate( (avec_add  , avec  [npolc*3:]) )
                           #
                           flds_add  [:npolc*3] += flds[:npolc*3]
                           dipind_add[:npolc*3] += flds[:npolc*3]
                           avec_add  [:npolc*3] += avec[:npolc*3]
                    #
                    assert len(flds_add  )==npols_add*3
                    assert len(dipind_add)==npols_add*3
                    assert len(avec_add  )==npols_add*3
                    #
                    shift_pol_mea += shift_pol_add_mea
                    del PAR_ADD, QO_ADD
                 # 
                 del PAR, QO
                 #
             
           # DISPERSION
           if self.__eval_disp:
              if num: print ' **** WARNING: Dispersion will not be evaluated!!!!'
              else:
                 npolc = parc['npol']
                 N = len(self.__ntp)
                 PAR = list(); PAR.append( parc )
                 # loop over solvent molecules in DISP layer (=POL layer)
                 if lwrite: print " DISP  LAYER: %10d molecules" % N
                 for i in xrange(N):
                     im = self.__mtp[i]
                     nm_prev = sum(self.__ntp[:i])
                     nm_curr = sum(self.__ntp[:i+1])
                     STR = self._reorder_xyz( self.__rcoordp[nm_prev:nm_curr]      , self.__reordlist[im] )                 
                     frg = self.__bsm[im].copy() 
                     rms = frg.sup( STR, suplist= self.__suplist[im] )
                     par = frg.get()
                     PAR.append( par )
                 npol = [ x['npol'] for x in PAR ]
                 npols= sum(npol)
                 rpol = numpy.concatenate  ([ x['rpol'] for x in PAR ]).reshape(npols*3)
                 pol  = numpy.concatenate  ([ x['dpoli'] for x in PAR ])#.reshape(npols*9)
                 pol  = pol.ravel()
                 rpol1= parc['lmoc1'].reshape(nmodes*npolc*3)
                 pol1 = parc['dpoli1'].reshape(nmodes*npolc*12*9)
                 disp_iso = libbbg.qm.clemtp .sdisp6(rpol,rpol1,pol,npol,pol1,gijj,redmss,freq,npolc,self.__mode)
                 disp_ani = libbbg.qm.clemtp2.tdisp6(rpol,rpol1,pol,npol,pol1,gijj,redmss,freq,npolc,self.__mode)
                 self.__fi_disp = None
                 shift_dis_mea = disp_ani
                 shift_dis_mea_iso = disp_iso
                 shift_dis_ea  = 0.0
                 #
                 del PAR
    
        # EXCHANGE-REPULSION
        rms_max = 0.0
        if self.__eval_rep:
           N = len(self.__nte)
           PAR = list(); PAR_SOL = list()
           # loop over solvent molecules in EXREP layer
           if lwrite: print " EXREP LAYER: %10d molecules" % N
           for i in xrange(N):
               im = self.__mte[i]
               nm_prev = sum(self.__nte[:i])
               nm_curr = sum(self.__nte[:i+1])
               STR = self._reorder_xyz( self.__rcoorde[nm_prev:nm_curr]   , self.__reordlist[im] )
               if num: STR = numpy.dot(STR+transl_inv, rot_inv)
               frg = self.__bsm[im].copy()
               rms = frg.sup( STR , suplist= self.__suplist[im])
               if lwrite: print " rms E: ",rms
               if rms > rms_max: rms_max = rms
               par = frg.get()
               if num: PAR_SOL.append( par )
               else:   PAR.append( par )
           if self.__rms_solvent_max is None: self.__rms_solvent_max = rms_max
           if num:
              func_rep = numpy.zeros(N_FD, numpy.float64)
              for i in xrange(N_FD):
                  varA = PAR_FD[i]
                  e_rep  = 0.0
                  for j in xrange(len(PAR_SOL)):
                      varB = PAR_SOL[j] 
                      e_rep += self._pair_rep(varA,varB)
                  func_rep[i] = e_rep
              # differentiate and obtain forces f and hessian K
              fx_1 = func_rep[0:self.__nata*6+1]
              fx_2 = func_rep[  self.__nata*6+1: ]
              fd_calc = libbbg.utilities.diff(func=None, step=0.006*self.AngstromToBohr, 
                                              DIM=self.__nata*3, scheme='3pt')
              fi_cart, kij_cart = fd_calc.eval( (fx_1, fx_2), symm=True)
              # transform the derivatives to normal coordinate space
              l = self.__lvec.reshape(self.__nmodes, self.__nata*3)
              fi_mode = numpy.dot(l, fi_cart)
              kij_mode= numpy.dot( numpy.dot(l, kij_cart), l.transpose() )
              # store forces and K-matrices
              self.__fi_rep = fi_mode
              self.__kij_rep = kij_mode
              # diagonalize the Hessian
              hess, freq, redmass, U = self._get_slv(type='x', theory=theory)
              libbbg.utilities.PRINT(freq[::-1])
              self.__hess_rep = hess
              # compute approximated frequencies at SolEFP-0 theory                          
              shift_rep_mea = self.get_shift_from_fi (fi_mode, iter=False, K4=None, kij=None)
              shift_rep_ea  = self.get_shift_from_kii(kij_mode)                              
              #
              del PAR_SOL
           else:
              # basis sets
              molA = libbbg.utilities.MakeMol(parc['atno'],parc['pos'])
              bfsA = PyQuante.Ints.getbasis(molA,parc['basis'])
              nbsa = parc['nbasis']
              nmosa= parc['nmos']
              nmodes = self.__nmodes
              # parameters for central molecule
              faij = parc['fock' ]               .ravel()
              faij1= parc['fock1']               .ravel()
              cika = parc['vecl' ]               .ravel()
              cika1= parc['vecl1']               .ravel()
              za   = parc['atno' ]
              rna  = parc['pos'  ]               .ravel()
              ria  = parc['lmoc' ]               .ravel()
              ria1 = parc['lmoc1']               .ravel()
              mlist= bfsA.get_bfsl() + 1
              redmss= self.__redmss
              gijj = self.__gijk[:,self.__mode-1,self.__mode-1]
              freq = self.__freq
              lvec = self.__lvec                 .ravel()
              RNB = list(); serp = 0.0
              for par in PAR:
                  molB = libbbg.utilities.MakeMol(par['atno'],par['pos'])
                  bfsB = PyQuante.Ints.getbasis(molB,par['basis'])
                  nbsb = par['nbasis']
                  nmosb = par['nmos']
                  # instantaneous integrals
                  skm  = PyQuante.Ints.getSAB(bfsA,bfsB)       .ravel()
                  tkm  = PyQuante.Ints.getTAB(bfsA,bfsB)       .ravel()
                  sk1m = PyQuante.Ints.getSA1B(bfsA,bfsB)      .ravel()
                  tk1m = PyQuante.Ints.getTA1B(bfsA,bfsB)      .ravel()
                  # molecule B
                  fbij = par['fock']             .ravel()
                  cikb = par['vecl']             .ravel()
                  zb   = par['atno']
                  rnb  = par['pos' ]           # .ravel() 
                  RNB.append(rnb)
                  rnb  = rnb.ravel()
                  rib  = par['lmoc']             .ravel()
                  # calculate the properties!
                  sij = numpy.zeros(nmosa*nmosb       , numpy.float64)
                  tij = numpy.zeros(nmosa*nmosb       , numpy.float64)
                  smij= numpy.zeros(nmodes*nmosa*nmosb, numpy.float64)
                  tmij= numpy.zeros(nmodes*nmosa*nmosb, numpy.float64)
                  fi  = numpy.zeros(nmodes            , numpy.float64)
                  #
                  sma, shftea = solvshift.shftex.shftex(redmss, freq, gijj, lvec,
                                ria, rib, rna, rnb, ria1, cika, cikb, cika1,
                                skm, tkm, sk1m, tk1m, za, zb, nbsb, mlist,
                                faij, fbij, faij1, self.__mode, sij, tij, smij, tmij, fi)
                  #eint = solvshift.exrep.exrep(ria.reshape(nmosa,3),
                  #                             rib.reshape(nmosb,3),
                  #                             rna.reshape(len(rna)/3,3),
                  #                             rnb.reshape(len(rnb)/3,3),
                  #                             faij.reshape(nmosa,nmosa),
                  #                             fbij.reshape(nmosb,nmosb),
                  #                             cika.reshape(nmosa,nbsa),
                  #                             cikb.reshape(nmosb,nbsb),
                  #                             skm.reshape(nbsa,nbsb),
                  #                             tkm.reshape(nbsa,nbsb),za,zb)
                  #self.eexrep = eint * self.HartreeToKcalPerMole
                  serp+= sma
                  self.__fi_rep += fi
                  #
              self.__solvent_exrep = RNB
              shift_rep_mea = serp
              shift_rep_ea  = 0.0
    
    
        # accumulate the total shift at SolEFP-0 level
        shift_total += shift_ele_mea + shift_ele_ea + shift_ele_corr_mea + shift_ele_corr_ea
        shift_total += shift_pol_mea + shift_pol_ea + shift_dis_mea + shift_dis_ea
        shift_total += shift_rep_mea + shift_rep_ea
    
        # change units from A.U. to specific units
        if self.__cunit:
           shift_ele_mea *= self.HartreePerHbarToCmRec; shift_ele_ea  *= self.HartreePerHbarToCmRec
           shift_pol_mea *= self.HartreePerHbarToCmRec; shift_pol_ea  *= self.HartreePerHbarToCmRec
           shift_dis_mea *= self.HartreePerHbarToCmRec; shift_dis_ea  *= self.HartreePerHbarToCmRec
           shift_rep_mea *= self.HartreePerHbarToCmRec; shift_rep_ea  *= self.HartreePerHbarToCmRec
           shift_pol_add_mea  *= self.HartreePerHbarToCmRec; shift_pol_add_ea  *= self.HartreePerHbarToCmRec
           shift_ele_corr_mea *= self.HartreePerHbarToCmRec; shift_ele_corr_ea *= self.HartreePerHbarToCmRec
           shift_dis_mea_iso  *= self.HartreePerHbarToCmRec
           shift_total        *= self.HartreePerHbarToCmRec
    
        # diagonalize total Hessians
        if num:
           # K-elect approximation
           hess, freq, redmass, U = self._get_slv(type='exp-app', theory=theory)
           if lwrite: libbbg.utilities.PRINT(freq[::-1])
           self.__hess_tot_app = hess
           # full K
           hess, freq, redmass, U = self._get_slv(type='exp', theory=theory)
           if lwrite: libbbg.utilities.PRINT(freq[::-1])
           self.__hess_tot = hess
           r,un = numpy.linalg.eig(hess)
           self.__u = un

        # update frequency shift dictionary
        self.__shift['ele_mea'  ] = shift_ele_mea
        self.__shift['ele_ea'   ] = shift_ele_ea     
        self.__shift['cor_mea'  ] = shift_ele_corr_mea
        self.__shift['cor_ea'   ] = shift_ele_corr_ea    ; self.__shift['ele_tot'] = shift_ele_mea + shift_ele_ea + shift_ele_corr_mea + shift_ele_corr_ea
        self.__shift['pol_mea'  ] = shift_pol_mea
        self.__shift['pol_ea'   ] = shift_pol_ea         ; self.__shift['pol_tot'] = shift_pol_mea + shift_pol_ea
        self.__shift['rep_mea'  ] = shift_rep_mea
        self.__shift['rep_ea'   ] = shift_rep_ea         ; self.__shift['rep_tot'] = shift_rep_mea + shift_rep_ea
        self.__shift['dis_mea'  ] = shift_dis_mea
        self.__shift['dis_ea'   ] = shift_dis_ea         ; self.__shift['dis_tot'] = shift_dis_mea + shift_dis_ea ; self.__shift['dis_mea_iso'] = shift_dis_mea_iso
        self.__shift['total'    ] = shift_total
        self.__shift['total_ele'] = shift_total - shift_rep_mea - shift_rep_ea
        self.__shift['total_mea'] = shift_ele_mea + shift_ele_corr_mea + shift_pol_mea + shift_rep_mea + shift_dis_mea
        self.__shift['total_ea' ] = shift_ele_ea  + shift_ele_corr_ea  + shift_pol_ea  + shift_rep_ea  + shift_dis_ea
        self.__shift['solcamm'  ] = shift_ele_mea + shift_ele_ea 
        self.__shift['total_cor'] = shift_ele_corr_mea + shift_ele_corr_ea
        self.__shift['pol_add_mea'] = shift_pol_add_mea
        self.__shift['pol_add_ea' ] = shift_pol_add_ea   ; self.__shift['pol_add_tot'] = shift_pol_add_mea + shift_pol_add_ea
    
        # final printout
        if lwrite:
           print    " Electrostatic      frequency shift: %10.2f" % (shift_ele_mea + shift_ele_ea + shift_ele_corr_mea + shift_ele_corr_ea)
           print    "                MEA           : %10.2f" % shift_ele_mea
           print    "                 EA           : %10.2f" % shift_ele_ea
           print    "               CORR           : %10.2f" %(shift_ele_corr_mea + shift_ele_corr_ea)
           print    " Polarization       frequency shift: %10.2f  %10.2f" % (shift_pol_mea + shift_pol_ea, shift_pol_add_mea + shift_pol_add_ea)
           if not num: 
              print " Dispersion         frequency shift: %10.2f  %10.2f" % (shift_dis_mea_iso, shift_dis_mea)
           print    " Exchange-repulsion frequency shift: %10.2f" % (shift_rep_mea + shift_rep_ea)
           print    " ----------------------------------------------- "
           print    " TOTAL FREQUENCY SHIFT             : %10.2f" %  shift_total
         
        return

    # REORDERING 

    def _reorder_xyz(self, xyz, reord_list):
        """Reorder the coordinates of MD atoms to EFP order. 
Keep track of unequal number of atoms between MD target and 
EFP parameter BSM (missing atoms in case of fragmented EFP superposition).
The convention is to place -1 in the reord_list for atoms that have to be removed""" 
        n = len(reord_list)
        xyz_ = numpy.zeros((n,3),numpy.float64)
        i = -1
        for I in reord_list:
            i+=1
            if I>-1: xyz_[i] = xyz[I].copy() 
            else:    xyz_[i] = numpy.zeros((1,3),numpy.float64)
        return xyz_

    # STRUCTURAL DISTORTIONS

    def _calc_dq(self, fi, theory, cart=False, kij=None):
        """calculate structural distortions according to the level of SolX theory"""
        if not theory in [0,2]:
           raise Exception('Incorrect level of theory used! Possible are 0 and 2 (---> chosen %i)'%theory)
        if not theory:
           m = self.__redmss*self.__freq*self.__freq
           dq = -fi/m
        elif theory==2: # requires kij matrix of size nmodes x nmodes
           m = self.__redmss*self.__freq*self.__freq
           N = self.__nmodes
           G = self.__gijk.copy()/2.
           B = numpy.diag(m) + kij
           B1= numpy.linalg.inv(B)
           I = numpy.identity(self.__nmodes, numpy.float64)
           a = numpy.tensordot(B1,G,(1,1))
           b = numpy.tensordot(a,B1,(2,0))
           A = numpy.tensordot(b,fi,(2,0))
           dq =-numpy.dot(numpy.dot(numpy.linalg.inv(I-A),B1),fi)

        if cart:
           # transform to Cartesian coordinates                        
           l = self.__lvec.reshape(self.__nmodes, self.__nata*3)
           dq = numpy.dot( l.transpose(), dq).reshape( self.__nata, 3)
           #temp = numpy.sqrt(self.__atms)[:, numpy.newaxis]
           #dq/= temp
        return dq

    # HESSIAN MATRICES

    def _calc_hess(self, type='e', mwc=True, theory=0):
        """Calculate the Hessian matrix from numerical evaluation of derivatives"""
        m = numpy.diag( self.__redmss * self.__freq**2.0 )
        if   type=='e': 
             fi  = self.__fi_el
             kij = self.__kij_el
        elif type=='p': 
             fi  = self.__fi_pol
             kij = self.__kij_pol
        elif type=='x': 
             fi  = self.__fi_rep
             kij = self.__kij_rep
        elif type=='exp':
             fi  = self.__fi_el + self.__fi_rep + self.__fi_pol 
             kij = self.__kij_el+ self.__kij_rep+ self.__kij_pol
        elif type=='exp-app':
             fi  = self.__fi_el + self.__fi_rep + self.__fi_pol
             kij = self.__kij_el.copy()
             #kij.fill(0.)

        if not theory: dq= self._calc_dq(fi, theory, cart=False, kij=None)
        else:          dq= self._calc_dq(fi, theory, cart=False, kij=kij )
        s = numpy.tensordot(self.__gijk, dq,(0,0))
        #t = numpy.outer(numpy.sqrt(self.__redmss), numpy.sqrt(self.__redmss))
        h = m + s + kij  #/ t
        h_w = None
        #h = numpy.diag(h.diagonal())
        if mwc:
           t = numpy.outer(numpy.sqrt(self.__redmss), numpy.sqrt(self.__redmss))
           h_w= h/t
           t = numpy.sqrt(self.__redmss)[:,numpy.newaxis]
           l = self.__lvec.reshape(self.__nmodes, self.__nata * 3) / t
           l = numpy.dot(self.__M, l.transpose() )
           h_w = numpy.dot(l,numpy.dot(h_w,l.transpose()))
        return h, h_w

    def _get_slv(self, type='e', theory=0):
        """Calculate the frequencies, reduced masses and transformation matrix from numerical derivatives"""
        hess, hess_w = self._calc_hess(type, mwc=True, theory=theory)
        vib = libbbg.utilities.VIB(self.__mol, hess_w, weight=False)
        vib.eval()
        freq, redmss, U = vib.get()
        return hess, freq, redmss, U

    def _remove_clashes(self, PAR, QO, lwrite, rcl_algorithm='remove_by_name'):
        """
 Eliminate clashes by some algorithm. Default is by removing residues from the list PAR and QO.
 see .eval method for polarization frequency shifts.
 """
        PAR_NEW = [PAR[0]]; QO_NEW = [QO[0]]
        PAR_ADD = [PAR[0]]; QO_ADD = [QO[0]]
        #
        if lwrite: print " CLASH REMOVAL PROCEDURE\n Algorithm: %s\n" % rcl_algorithm
        i = 1
        if rcl_algorithm == 'remove_by_name':
           removable_pars = ['Methane', 'Ethane', 'N-Propane', 'N-methylacethamide', 'N-methylacethamide D7'] 
           removable_pars+= ['Benzene', 'Phenol', 'Formamide', 'Methane', 'Methyl amonium cation (+1)', 'Acetate anion (-1)', 'Methyl Guanidinium Cation (+1)', 'Dimethyl Sulfide'  ]
           for par in PAR[1:]:
               if lwrite: print "%30s" % par['name'], 
               if par['name'] not in removable_pars: 
                  if lwrite: print " - Included"
                  PAR_NEW.append(PAR[i])
                  QO_NEW.append(QO[i])
               else: 
                  if lwrite: print " - Rejected"
                  PAR_ADD.append(PAR[i])
                  QO_ADD.append(QO[i])
               i+=1
        return PAR_NEW, QO_NEW, PAR_ADD, QO_ADD

    # PAIR ENERGIES
    
    def _pair_elect(self,varA,varB):
        """MTP electrostatic pair energy"""
        return
    
    def _pair_rep(self,varA,varB):
        """exchange-repulsion pair energy"""
        # basis sets
        molA = libbbg.utilities.MakeMol(varA['atno'],varA['pos'])
        molB = libbbg.utilities.MakeMol(varB['atno'],varB['pos'])
        bfsA = PyQuante.Ints.getbasis(molA,varA['basis'])
        bfsB = PyQuante.Ints.getbasis(molB,varB['basis'])
        # instantaneous integrals
        skm  = PyQuante.Ints.getSAB(bfsA,bfsB)
        tkm  = PyQuante.Ints.getTAB(bfsA,bfsB)
        # parameters
        ### molecule A
        faij = varA['fock']
        cika = varA['vecl']
        za   = varA['atno']
        rna  = varA['pos']
        ria  = varA['lmoc']
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # calculate the properties!
        eint = solvshift.exrep.exrep(ria,rib,rna,rnb,faij,fbij,cika,cikb,skm,tkm,za,zb)
        return eint

    def _pair_rep_freq(self,varA,varB):
        """exchange-repulsion pair frequency shift"""
        # basis sets
        molA = libbbg.utilities.MakeMol(varA['atno'],varA['pos'])
        molB = libbbg.utilities.MakeMol(varB['atno'],varB['pos'])
        bfsA = PyQuante.Ints.getbasis(molA,varA['basis'])
        bfsB = PyQuante.Ints.getbasis(molB,varB['basis'])
        # instantaneous integrals
        skm  = PyQuante.Ints.getSAB(bfsA,bfsB)
        tkm  = PyQuante.Ints.getTAB(bfsA,bfsB)
        sk1m = PyQuante.Ints.getSA1B(bfsA,bfsB)
        tk1m = PyQuante.Ints.getTA1B(bfsA,bfsB)
        # parameters
        ### molecule A
        faij = varA['fock']
        faij1= varA['fock1']
        cika = varA['vecl']
        cika1= varA['vecl1']
        za   = varA['atno']
        rna  = varA['pos']
        ria  = varA['lmoc']
        ria1 = varA['lmoc1']
        mlist= bfsA.get_bfsl() + 1
        redmss= varA['redmass']
        gijj = varA['gijk'][:,self.__mode-1,self.__mode-1]
        freq = varA['freq']
        lvec = varA['lvec']
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # calculate the properties!
        shift ,shftea = solvshift.shftex.shftex(redmss,freq,gijj,lvec,
                                      ria,rib,rna,rnb,ria1,
                                      cika,cikb,cika1,
                                      skm,tkm,sk1m,tk1m,
                                      za,zb,mlist,
                                      faij,fbij,faij1,self.__mode)
        return shift

    # to be deprecated

    def _xyz(self, pos, atno, units='Angstrom'):
        """generate the xyz file contents (without 2 first lines). Copied from Frag class - so it should be deprecated in the future to improve code stability"""
        atomic_symbols = {1: 'H', 2: 'He', 3: 'Li', 6:'C', 7:'N', 8:'O', 9:'F', 11:'Na', 12:'Mg', 16:'S', 17:'Cl'}
        log = ''
        r = pos.copy()
        Z = atno
        if units.lower().startswith('angs'): r *= self.BohrToAngstrom
        for i,x in enumerate(r):
            log += '%3s' % atomic_symbols[Z[i]]
            log += 3*'%16.6f' % tuple(x)
            log += '\n'
        return log 

    def _eval_dq(self, fi, theory, cart=True):
        """calculate structural distortions according to the level of SolX theory"""
        m = self.__redmss*self.__freq*self.__freq
        dq = -fi/m
        if not theory in [0,2]:
           raise Exception('Incorrect level of theory used! Possible are 0 and 2 (---> chosen %i)'%theory)
        if theory==2:
           # 1 - good
           #M = numpy.diag(1./m)
           #g = self.__gijk
           #a = numpy.tensordot(M,g,(1,1))
           #b = numpy.tensordot(a,M,(2,0))
           #c = numpy.tensordot(b,fi,(2,0))
           #A = c/2.

           # 2 - not good (slightly different result than in 1)
           #A  = numpy.zeros((self.__nmodes, self.__nmodes), numpy.float64)
           #for i in range(15):
           #    for j in range(15):
           #        gk = 0.0
           #        for k in range(15):
           #            gk+= self.__gijk[i,j,k] * fi[k]/m[k]
           #        A[i,j] = gk/2.
           #    A[i,j] /= m[i]

           # 3 - incorrect, bad
           #g  = (self.__gijk / m * fi).sum(axis=2)
           #A  = g/m/2.
           #dq = numpy.dot( numpy.linalg.inv(numpy.identity(self.__nmodes)-A), dq)

           N = self.__nmodes
           G = self.__gijk.copy()/2.
           M = numpy.diag(self.__redmss*self.__freq**2)
           F = fi
           M1= numpy.linalg.inv(M)
           I = numpy.identity(self.__nmodes, numpy.float64)
           A = numpy.tensordot(M1,G,(1,1))
           B = numpy.tensordot(A,M1,(1,0))
           a = numpy.tensordot(B,F,(2,0))
           #dQ = dot(linalg.inv(M-dot(M,a)),F)
           dq = numpy.dot(numpy.dot(numpy.linalg.inv(I-a),M1),F)

        if cart:
           # transform to Cartesian coordinates                        
           l = self.__lvec.reshape(self.__nmodes, self.__nata*3)
           dq = numpy.dot( l.transpose(), dq).reshape( self.__nata, 3)
           #temp = numpy.sqrt(self.__atms)[:, numpy.newaxis]
           #dq/= temp
        return dq

        


class EFP_pair(object,libbbg.units.UNITS):
    """
=============================================================================
              EFFECTIVE FRAGMENT POTENTIAL METHOD FOR A DIMER                
=============================================================================

Usage:
A = EFP(a,b)
result = A(ct=True,nlo=False,shift=False,nmode=None,cunit=False)
rms_a, rms_b = A.sup(str_a=None,str_b=None)

Notes:
1) a and b are SLVPAR instances. If shift=True, a denotes IR-active molecule
   and should contain all necessary parameters for it
2) the a and b objects are assumed to be appropriately transformed in space 
   by rotations and translations
3) a and b in sup argument list are ndarray structures of dimension (natoms,3)
   The coordinates are assumed to be in A.U.
4) cunit - change units. If <True> then energies are returned in [kcal/mole]
   and frequency shifts in [cm-1]. Otherwise all is returned in A.U. units
"""
    def __init__(self,a,b):
        self.__molA = a
        self.__molB = b
        return
    
    def __call__(self,ct=True,nlo=False,shift=False,nmode=None,cunit=False):
        """perform all the calculation"""
        if ct: 
             ma, ea = self._ex_rep_ct(nlo,shift,nmode,cunit)
        else:
             ma, ea = self._ex_rep(nlo,shift,nmode,cunit)
        return ma, ea
    
    def __repr__(self):
        """print the two fragments"""
        log = '\n'
        log+= str(self.__molA.get_pos()*self.BohrToAngstrom)
        log+= '\n'
        log+= str(self.__molB.get_pos()*self.BohrToAngstrom)
        log+= '\n'
        return str(log)
    
    def sup(self,str_a=None,str_b=None):
        """superimpose the a and b objects to given structures"""
        rms_a, rms_b = None, None
        if str_a is not None: rms_a = self.__molA.sup(str_a)
        if str_b is not None: rms_b = self.__molB.sup(str_b)
        return rms_a, rms_b
    
    def get(self):
        """return a tuple of parameters for two fragments"""
        return self.__molA.get(), self.__molB.get()
        
    def get_mol(self):
        """returns molecular fragments"""
        return self.__molA, self.__molB
    
    # protected
    def _ex_rep_ct(self,nlo,shift,nmode,cunit):
        """calculate exchange-repulsion and charge-transfer properties"""
        # variables
        varA = self.__molA.get()
        varB = self.__molB.get()
        # basis sets
        bfsA = self.__molA.get_bfs()
        bfsB = self.__molB.get_bfs()
        # parameters
        ### molecule A
        faij   = varA['fock']
        faij1  = varA['fock1']
        cika   = varA['vecl']
        cika1  = varA['vecl1']
        za     = varA['atno']
        rna    = varA['pos']
        ria    = varA['lmoc']
        ria1   = varA['lmoc1']
        faijc  = varA['fckc'].diagonal()
        cikca  = varA['vecc']
        qa     = varA['chlpg']
        redmss = varA['redmass']
        freq   = varA['freq']
        gijj   = varA['gijk'][nmode-1,nmode-1,:]
        eiglvc = varA['lvec']
        mlist  = bfsA.get_bfsl() + 1
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        fbijc= varB['fckc'].diagonal()
        cikcb= varB['vecc']
        qb   = varB['chlpg']
        # instantaneous integrals
        #import time
        #t0 = time.time()
        skm  = PyQuante.Ints.getSAB(bfsA,bfsB)
        tkm  = PyQuante.Ints.getTAB(bfsA,bfsB)
        #t1 = time.time()
        sk1m = PyQuante.Ints.getSA1B(bfsA,bfsB)
        tk1m = PyQuante.Ints.getTA1B(bfsA,bfsB)
        tkk =  PyQuante.Ints.getT(bfsA)
        tll =  PyQuante.Ints.getT(bfsB)
        #t2 = time.time()
        vkl =  PyQuante.Ints.getVEFP(bfsA,bfsB,qb,rnb)
        vlk =  PyQuante.Ints.getVEFP(bfsB,bfsA,qa,rna)
        vkm =  PyQuante.Ints.getVEFP(bfsA,bfsA,qb,rnb)
        vln =  PyQuante.Ints.getVEFP(bfsB,bfsB,qa,rna)
        #t3 = time.time()
        #vvv = getV(bfsA,self.__molA.atoms)
        #t4 = time.time()
        #print "STAB %.2f" % (t1-t0)
        #print "XA1B %.2f" % (t2-t1)
        #print "VXLN %.2f" % (t3-t2)
        #print "TOTL %.2f" % (t3-t0)
        #print "DUPA %.2f" % (t4-t3)
        # calculate the properties!
        shftma,shftea = solvshift.shftce.shftce(redmss,freq,gijj,eiglvc,
                                      rna,rnb,ria,rib,ria1,                      
                                      cika,cikb,cikca,cikcb,cika1,
                                      skm,tkm,tkk,tll,vkm,vkl,vlk,vln,sk1m,tk1m,
                                      faij,fbij,faijc,fbijc,faij1,
                                      za,zb,mlist,nmode)
        if cunit:
           shftma *= self.HartreePerHbarToCmRec
           shftea *= self.HartreePerHbarToCmRec
           
        return shftma, shftea

    def _ex_rep(self,nlo,shift,nmode,cunit):
        """calculate exchange-repulsion property"""
        # variables
        varA = self.__molA.get()
        varB = self.__molB.get()
        # basis sets
        bfsA = self.__molA.get_bfs()
        bfsB = self.__molB.get_bfs()
        # instantaneous integrals
        skm  = PyQuante.Ints.getSAB(bfsA,bfsB)
        tkm  = PyQuante.Ints.getTAB(bfsA,bfsB)
        sk1m = PyQuante.Ints.getSA1B(bfsA,bfsB)
        tk1m = PyQuante.Ints.getTA1B(bfsA,bfsB)
        # parameters
        ### molecule A
        faij   = varA['fock']
        faij1  = varA['fock1']
        cika   = varA['vecl']
        cika1  = varA['vecl1']
        za     = varA['atno']
        rna    = varA['pos']
        ria    = varA['lmoc']
        ria1   = varA['lmoc1']
        redmss = varA['redmass']
        freq   = varA['freq']
        gijj   = varA['gijk'][nmode-1,nmode-1,:]
        lvec   = varA['lvec']
        mlist  = bfsA.get_bfsl() + 1
        ### molecule B
        fbij = varB['fock']
        cikb = varB['vecl']
        zb   = varB['atno']
        rnb  = varB['pos']
        rib  = varB['lmoc']
        # calculate the properties!
        shftma,shftea = solvshift.shftex.shftex(redmss,freq,gijj,lvec,
                                      ria,rib,rna,rnb,ria1,  
                                      cika,cikb,cika1,
                                      skm,tkm,sk1m,tk1m,
                                      za,zb,mlist,
                                      faij,fbij,faij1,nmode)
        if cunit:
           shftma *= self.HartreePerHbarToCmRec
           shftea *= self.HartreePerHbarToCmRec

        return shftma, shftea

 
    
class FragFactory(object,diff.DIFF):
    """
Solvatochromic Effective Fragment Potential Fragment
----------------------------------------------------

Usage:
"""
    def __init__(self,anh=None,basis=None,nae=None,
                      fchk=None,gmslog=None,
                      chelpg=None,esp=None,):
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
        self.__chlpg  = chelpg
        self.__esp    = esp
        self._init()
        self._create()

    # public
    
    def set(self,anh=None,basis=None,nae=None,
                 fchk=None,gmslog=None,
                 chelpg=None,esp=None,
                 dpol=None,dpol1=None,):
        """set the properties to the object"""
        if self.__anh    is not None: self.__anh    = anh
        if self.__fchk   is not None: self.__fchk   = fchk
        if self.__gmslog is not None: self.__gmslog = gmslog
        if self.__basis  is not None: self.__basis  = basis
        if self.__nae    is not None: self.__nae    = nae
        if self.__chlpg  is not None: self.__chlpg  = chelpg
        if self.__esp    is not None: self.__esp    = esp
        if self.__dpol   is not None: self.__dpol   = dpol
        if self.__dpol1  is not None: self.__dpol1  = dpol1
        return

    def reset(self,anh=None,basis=None,nae=None,
                   fchk=None,gmslog=None,
                   chelpg=None,esp=None,
                   dpol=None,dpol1=None,):
        """reset the properties"""
        self.__anh    = anh
        self.__fchk   = fchk
        self.__gmslog = gmslog
        self.__basis  = basis
        self.__nae    = nae
        self.__chlpg  = chelpg
        self.__esp    = esp
        self.__dpol   = dpol
        self.__dpol1  = dpol1
        return

    def get(self):
        """returns dictionary with parameters"""
        par = {}
        if self.__lmoc  is not None: par['lmoc' ] = self.__lmoc
        if self.__lmoc1 is not None: par['lmoc1'] = self.__lmoc1
        if self.__fock  is not None: par['fock' ] = self.__fock
        if self.__fock1 is not None: par['fock1'] = self.__fock1
        if self.__vecl  is not None: par['vecl' ] = self.__vecl
        if self.__vecl1 is not None: par['vecl1'] = self.__vecl1
        if self.__vecc  is not None: par['vecc' ] = self.__vecc
        if self.__vecc1 is not None: par['vecc1'] = self.__vecc1
        if self.__fckc  is not None: par['fckc' ] = self.__fckc
        if self.__fckc1 is not None: par['fckc1'] = self.__fckc1
        if self.__esp   is not None: par['esp'  ] = self.__esp
        if self.__chlpg is not None: par['chlpg'] = self.__chlpg
        if self.__dpol  is not None: par['dpol' ] = self.__dpol
        if self.__dpol1 is not None: par['dpol1'] = self.__dpol1
        if self.__rpol  is not None: par['rpol' ] = self.__rpol
        return par
    
    def eval(self, ct=False, cvgloc=1.0E-14):
        """Parses AO-LMO transformation matrix and Fock matrix.
Transforms the latter from AO to LMO space. Computes also 
overlap integrals and parses density matrix. If ct is True
the canonical Fock matrix and vectors will be saved."""
        assert self.__mol is not None, 'molecule not specified! (no fchk file)'
        # evaluate transformation matrices and LMO centroids
        SAO   = PyQuante.Ints.getS(self.__bfs)
        dmat  = libbbg.utilities.ParseDmatFromFchk(self.__fchk,self.__basis_size)
        vecc  = libbbg.utilities.ParseVecFromFchk(self.__fchk)
        veccocc= vecc.copy()[:self.__nae,:]
        tran, veclmo = libbbg.utilities.get_pmloca(self.__natoms,mapi=self.__bfs.LIST1,sao=SAO,
                                            vecin=veccocc,nae=self.__nae,
                                            maxit=100000,conv=cvgloc, 
                                            lprint=True,
                                            freeze=None)
        # calculate LMTPs
        camm = coulomb.multip.MULTIP(molecule=self.__mol,
                                     basis=self.__basis,
                                     method='b3lyp',
                                     matrix=dmat,
                                     transition=False,
                                     bonds=None,vec=veclmo)
        camm.camms()
        dma = camm.get()[0]
        # parse Fock matrix
        if self.__gmslog is not None:
           Fock = libbbg.utilities.ParseFockFromGamessLog(self.__gmslog,interpol=False)
        else:
           epsi = numpy.diag(libbbg.utilities.ParseAlphaOrbitalEnergiesFromFchk(self.__fchk))
           Fock = numpy.dot(vecc.T, numpy.dot(epsi, vecc))
           Fock = numpy.dot(numpy.dot(SAO, Fock), SAO)
        fock = numpy.tensordot(veclmo,numpy.tensordot(veclmo,Fock,(1,0)),(1,1))
        if ct: fckc = numpy.tensordot(vecc  ,numpy.tensordot(vecc  ,Fock,(1,0)),(1,1))
        # save
        self.__lmoc = dma.get_origin()[self.__natoms:]
        self.__rpol = self.__lmoc.copy()
        self.__tran = tran
        self.__vecl = veclmo
        if ct: self.__vecc = vecc
        self.__fock = fock
        if ct: self.__fckc = fckc
        self.__sao  = SAO
        self.__dmat = dmat
        self.__dma  = dma
        return
    
    # protected
    
    def _init(self):
        """initialize the other memorials"""
        self.__mol    = None; self.__basis_size = None
        self.__lmoc   = None; self.__lmoc1  = None; self.__natoms=None
        self.__fock   = None; self.__fock1  = None; self.__bfs  = None
        self.__vecl   = None; self.__vecl1  = None; self.__vecc = None
        self.__vecc1  = None; self.__fckc   = None; self.__fckc1= None
        self.__chlpg  = None; self.__esp    = None; self.__dpol = None
        self.__dpol1  = None; self.__rpol   = None
        return
    
    def _create(self):
        """creates the molecule"""
        if self.__fchk is not None:
           mol  = libbbg.utilities.Read_xyz_file(self.__fchk,mol=True,
                                          mult=1,charge=0,
                                          name='happy dummy molecule',
                                          basis=self.__basis)
           bfs        = PyQuante.Ints.getbasis(mol,self.__basis)
           basis_size = len(bfs)
           natoms= len(mol.atoms)
           # save
           self.__mol = mol
           self.__bfs = bfs
           self.__basis_size = basis_size
           self.__natoms = natoms
        return   
