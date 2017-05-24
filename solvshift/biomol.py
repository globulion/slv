#!*-* coding=utf-8 *-*
# -------------------------------------------------------------------------- #
#                      SOLEFP FOR BIOMOLECULES MODULE                        #
# -------------------------------------------------------------------------- #
import solvshift.slvpar
import solvshift.solefp
import solvshift.md
import libbbg.utilities
import libbbg.dma
import MDAnalysis
import numpy
import math
import sys
import re
sys.stdout.flush()

__all__     = [ 'BiomoleculeFragmentation', ]
__version__ = '1.0.1'
__author__  = 'Bartosz Błasiak'

class BiomoleculeFragmentation(object):
    """
 =============================================================================================
                     The SolEFP/EFP2 Biomolecular Fragmentation Scheme.

            Bartosz Błasiak                          Created: Summer 2015 (class form: 2017)
 =============================================================================================

 Citation: Blasiak, Ritchie, Webb, Cho; Phys. Chem. Chem. Phys. Vol.18 (2016), pp.18094-18111.

 Description: 

 Performs the vibrational solvatochromism SolEFP calculations on protein 
 in solution based on the Biomolecule EFP-Fragmentation Scheme and Rigid Molecule
 Algorithm. The protein and solvent are described by EFP2 fragments. Solvent is treated
 as in conventional SolEFP method. Protein is divided into two zones: 

  1) Near zone - atoms directly adjacent to the IR probe residue (the closest amide I subunits)

  2) Far zone - all other protein atoms (side chains and other amide I subunits). 

 Aminoacid side chains are modelled by sets of EFP2 fragments. Peptide groups (backbone) are approximated by 
 N-methylacetamide (NMA) EFP2's except for the nearest peptide subunits which are directly bonded
 to IR probe residue (near zone amide subunits). Frequency shifts due to this near zone are modelled
 by interaction of SolCAMM with 4-point charge model of CONH groups based on ESP fitting 
 derived from NMA molecule in vacuo. Coulombic, exchange-repulsion and dispersion frequency shifts
 are computed for all fragments as in conventional SolEFP-cutoff scheme. Induction is treated
 in a special way: many-body effects on induction are accounted for only for the set of EFP2's
 that are not in spacial clash with each other. The rest of EFP2's are modelled in a pairwise
 manner. The total frequency shift is written as:

          dw = dw(SolCAMM-CONH) + dw(SolEFP) + dw(through-bond) + dw(error)                     (1)

 where 
      - dw(SolCAMM-CONH)         frequency shift due to near zone

      - dw(SolEFP)               SolEFP frequency shift due to far zone

      - dw(through-bond)         is a through-bond correction due to the covalent bonding 
                                 of IR probe framework to the protein body. For SCN CN stretch
                                 mode this correction is -3.1 cm-1.

      - dw(error)                accounts for the systematic error of SolEFP model. This error
                                 is predicted to be approximately -8.0 and -15 cm-1 for 
                                 CN stretch and amide I modes, respectively.       

 The SolEFP frequency shift can be further decomposed as

         dw(SolEFP) = sum_(all EFP2's) [ dw(Coul) + dw(Ex-Rep) + dw(Disp) ] + dw(Ind)           (2)

 where 

         dw(Ind) = dw(Ind)[Unclashing EFP2's] + sum_(remaining EFP2's) dw(Ind)                  (3)

 To determine which EFP2's can be used in many-body calculations of induction effect the 
 clash removal algorithm is performed in each MD step. This algorithm has been implemented
 in the solvshift.solefp.EFP class.

 ---------------------------------------------------------------------------------------------

 Usage: 

 method = BiomoleculeFragmentation(res, probe, mode,
                                   elect=True, polar=True, repul=True, disp=True, correc=True,
                                   ccut=20.0, pcut=13.0, ecut=11.0,
                                   solcamm=None, lprint=False, is_probe_terminal=False,
                                   report='report.dat', write_solefp_input=False)

 method.run(top, traj, nframes=10)

 ---------------------------------------------------------------------------------------------

 The global settings for a computation:
 ccut         - Coulomb cutoff
 pcut         - Polarization cutoff
 ecut         - Exchange-repulsion cutoff
 elect        - evaluate electrostatics
 pol          - evaluate polarization
 rep          - evaluate exchange-repulsion
 disp         - evaluate dispersion

 =============================================================================================
                                                                 Last Revision: 9 May 2017
"""
    def __init__(self, res, probe, mode, 
                       elect=True, polar=True, repul=True, disp=True, correc=True,
                       ccut=20.0, pcut=13.0, ecut=11.0,
                       solcamm=None, lprint=False, is_probe_terminal=False,
                       report='report.dat', write_solefp_input=False,write_debug_file=False):
        # printout options
        self.report            = report
        self.lprint            = lprint
        self.write_solefp_input= write_solefp_input
        self.write_debug_file  = write_debug_file
        self.ccut = ccut; self.pcut = pcut; self.ecut = ecut
        # MD topology - EFP fragmentation relationship
        self.__res_biomolecule = self._parse_res(res)
        self.__res_probe       = self._parse_res(probe)
        self.__probe_label     = self.__res_probe.keys()[0]
        # construct solefp fragment
        self.__solefp          = solvshift.slvpar.Frag(self.__res_probe[self.__probe_label][0][0])
        # read the SolCAMM data (probe)
        self.__solcamm         = self._parse_solcamm(solcamm, mode)
        # intrinsic data structures
        self._init_data()
        # initialize the SolEFP calculator
        self._init_solefp(elect, polar, repul, disp, correc, ccut, pcut, ecut, mode)
        # normal mode (helico)
        self.__mode            = mode
        # whether IR probe is located at terminus or is inside peptide chain
        self.__is_probe_terminal = is_probe_terminal
        # utility method
        self._rr               = lambda a, b: math.sqrt(sum((a-b)**2))
        pass

    def run(self, top, traj, dframes=1, nframes=10, conh='nma', out_inp='biomolecule.sol'):
        """Analyse the MD trajectory"""
        # [0] Read the topology and trajectory
        universe = MDAnalysis.Universe(top, traj)
        system   = universe.selectAtoms('all')

        if self.lprint:
           print "\n Starting SolEFP calculations from Protein Fragmentation Scheme"
           print   " IR probe: %s (model: %s)" % (self.__probe_label, self.__solefp.get_name())
           if     self.__is_probe_terminal: print " Note: Assuming IR probe is located at one of protein termini."
           if not self.__is_probe_terminal: print " Note: Assuming IR probe is located inside protein chain."
           print "    SolEFP calculations setup:"
           print "                      Coul   Ind   Disp  ExRep"; c = libbbg.units.UNITS.BohrToAngstrom
           print "       Cutoff [Bohr]  %4.1f  %4.1f  %4.1f  %4.1f" % (self.ccut  , self.pcut  , self.pcut  , self.ecut  )
           print "       Cutoff [Angs]  %4.1f  %4.1f  %4.1f  %4.1f" % (self.ccut*c, self.pcut*c, self.pcut*c, self.ecut*c)

        # [1] Iterate over the MD trajectory
        for ts in universe.trajectory:
            if self.lprint: " Reading %d frame..." % ts.frame

            # [2] Find the amide subunits in the first step and create SolEFP input 
            if ts.frame == 1: 
               amide_atoms_far, amide_atoms_close = self._find_amide_atoms(system)
               residues = system.residues
               if self.lprint:
                  print " Note: There are %d amide units detected"    % (len(amide_atoms_far) + len(amide_atoms_close)) 
                  print "       Including %d Near Zone amide units\n" %  len(amide_atoms_close)

               # [3] Create temporary SolEFP input file
               if self.lprint: print " Creating SolEFP input..."
               self._create_solefp_inp(residues, amide_atoms_far, amide_atoms_close, conh, out_inp)

               # [4] Creating MD Input objects
               self._create_mdinput(self.__mdinput_log)

            # [5] Update SolCAMM and CONH2 DMA objects
            self._update_dma(system, amide_atoms_close)

            # [6] Evaluate frequency shifts
            log = self._eval_solefp(universe.trajectory, system)

            # [7] Save results on disk
            self.__out.write(log)
            
            # [8] Terminate
            if ts.frame == nframes: 
               if self.lprint: print " WARNING: Maximum number of %d frames reached. Quitting..." % nframes
               sys.exit()
        #
        self.__out.close()
        return

    # --- Intrinsic methods

    def _create_mdinput(self, mdinput):
        self.__args, self.__idx = solvshift.md.MDInput(mdinput).get()
        return

    def _eval_solefp(self, trajectory, selection):
        """Perform single SolEFP calculation step"""
        xyz       = selection.get_positions() * libbbg.units.UNITS.AngstromToBohr
        # 
        self.__SolEFP_Calculator.set(xyz[self.__idx], *self.__args)
        self.__SolEFP_Calculator.eval_dma(self.__conhs_dma, lwrite=self.lprint)
        self.__SolEFP_Calculator.eval(self.lprint+self.write_debug_file, remove_clashes=True)
        rms_c = self.__SolEFP_Calculator.get_rms()
        rms_s = self.__SolEFP_Calculator.get_rms_sol()
        rms_a = self.__SolEFP_Calculator.get_rms_ave()
        dw    = self.__SolEFP_Calculator.get_shifts()
        DW    = [ dw[x] for x in self.__dw_list ]
        DW_ENV= [ dw[x] for x in self.__dw_env_list ] 
        # line of output frequencies and misc data
        log = "%8.4f"                 % trajectory.time
        log+= len(DW)     * " %10.2f" % tuple(DW)
        log+= len(DW_ENV) * " %10.2f" % tuple(DW_ENV)
        log+= " %10.5f %10.5f %10.5f" % (rms_c, rms_s, rms_a)
        log+= '\n'
        return log

    def _update_dma(self, sel, amides_close):
        # update SolCAMM object
        if True:
           dummy, ir_indices, suplist = self.__res_probe[self.__probe_label][0]
           ir_indices = ir_indices[suplist-1] - 1
           xyz_probe  = sel.selectAtoms('resname %s' % self.__probe_label).get_positions() *libbbg.units.UNITS.AngstromToBohr
           xyz_sup    = xyz_probe[ir_indices]
           rms        = self.__solcamm.sup(xyz_sup, suplist=suplist)
           if self.lprint: print " Note: SolCAMM superimposition RMS: %10.5f" % rms

        # update CONH2 DMA object
        xyz_close = sel.get_positions()[amides_close[0]]
        xyz_close = numpy.concatenate((xyz_close, sel.get_positions()[amides_close[1]])) *libbbg.units.UNITS.AngstromToBohr
        self.__conhs_dma.set_structure(pos=xyz_close, equal=True, atoms=','.join(2*['C','O','N','H']))
        return

    def _create_solefp_inp(self, residues, amide_atoms_far, amide_atoms_close, conh, out_inp):
        """find (Sol)EFP fragments for IR probe, far zone amide subunits and all side chains"""
        if self.write_solefp_input:
           inp = open(out_inp, 'w')
        log = ''
        # SolEFP probe
        for residue in residues:
            if residue.name==self.__probe_label:
               probe_key = residue.name+str(residue.id)
               if self.lprint: print " Note: IR Probe is located in %s" % probe_key
               log += self._ap(self.__res_probe[self.__probe_label], self.__probe_label, residue.indices()+1)
        
        # Write the input for non-IR probe fragments
        res_logs = dict() 

        # Far amide units
        for residue in amide_atoms_far:
            if conh == 'nma':
               res_data = (('nma', [0,3,1,0,2,0,0,4,0,0,0,0], [3,5,2,8]),)           
               atom_numbers = numpy.arange(residue[0],residue[0]+12)+1
               name = 'CONH'
               res_logs[name+str(residue)] = self._app(res_data, name, atom_numbers)
            else: 
               raise NotImplementedError, " Other models of CONH group than NMA are not implemented yet!"
       
        # Aminoacid sidechains
        for residue in residues:
            try: 
                res_data = self.__res_biomolecule[residue.name]
                atom_numbers = residue.indices() + 1
            except KeyError:
                if residue.name != self.__probe_label:
                   if self.lprint: print " WARNING: Residue %5s not found in topology conversion file" % residue.name
            #
            key = residue.name+str(residue.id)
            if key not in res_logs.keys():                                                     
               res_logs[key] = self._ap(res_data, residue.name, atom_numbers) 
            else:                                                                              
               res_logs[key]+= self._qp(res_data, residue.name, atom_numbers)                        
        
        for residue, res_log in sorted(res_logs.items()):
            if residue!=probe_key:
               log += '! === %s RESIDUE === \n' % residue 
               log += res_log
        #
        if self.write_solefp_input:
           inp.write(log) 
           inp.close()
        # save
        self.__mdinput_log = log
        return 

    # --- Determine amide subunits closest to IR probe (Near Zone), and further from it (Far Zone)
    def _find_amide_atoms(self, sel):
        """Returns the atomic indices of CONH groups in far and close regions from IR_PROBE residue"""
        amides= list()
        atoms = sel.names()
        xyz   = sel.get_positions()
        vic   = sel.selectAtoms('around 2.4 (resname %s and name CA)' % self.__probe_label)
        close_atoms_indices = vic.indices()
    
        for i in xrange(len(atoms)):
            if atoms[i]=='C': 
               try:
                  if atoms[i+1]=='O':
                     if atoms[i+2]=='N':
                        if atoms[i+3]=='H':
                           amides.append([i,i+1,i+2,i+3])
               except IndexError:
                  print " Omiting %d-th 'C' atom" % (i+1)
                  continue
        # remove improper choices of CONH 
        i = 0; chosen_indices = list()
        for conh in amides:
            c,o,n,h = xyz[conh]
            if self._rr(c,n) > 2.0: 
               chosen_indices.append(i)
            i+=1
        amides = numpy.array(amides)
        chosen_indices = numpy.array(chosen_indices,int)
        try:
           amides_detected = libbbg.utilities.choose(amides,chosen_indices)
        except IndexError:
           amides_detected = amides.copy()
           if self.lprint: print " Note: No improper CONH pairs found."
        amides_rejected = amides[chosen_indices]
        amides = list(amides_detected)
        if self.lprint: 
           if len(amides_rejected): print " Note: Removed %d improper CONH pairs." % len(chosen_indices)
        
        # From among the rejected CONH groups compare CO and NH
        # and determine if some amide groups were missing
        UNDETECTED_AMIDES=0
        for i in xrange(len(chosen_indices)):
            conh_i = xyz[amides_rejected[i]]
            for j in xrange(len(chosen_indices)):
                conh_j = xyz[amides_rejected[j]]
                if i!=j:
                   ci,oi,ni,hi = conh_i
                   cj,oj,nj,hj = conh_j
                   if self._rr(ci,nj) < 2.0: UNDETECTED_AMIDES+=1
                   if self._rr(cj,ni) < 2.0: UNDETECTED_AMIDES+=1
    
        if UNDETECTED_AMIDES: print " WARNING!!! %d CONH groups were not included due to deviations from atom numbering!" % UNDETECTED_AMIDES
    
        # remove CONH units very close to the IR probe
        i = 0; close_amide_indices = list()
        for conh in amides:
            for x in conh:
                if x in close_atoms_indices: 
                   close_amide_indices.append(i)
                   break
            i+=1
        amides = numpy.array(amides, int)
        if self.__is_probe_terminal:
           assert len(close_amide_indices)==1, \
 " It must be just 1 CONH group close to %s because IR probe is assumed to be located in one of protein termini!" % self.__probe_label
        else:
           assert len(close_amide_indices)==2, \
 " It must be just 2 CONH groups close to %s because IR probe is assumed to be inside a peptide chain!" % self.__probe_label
        amides_far = list(libbbg.utilities.choose(amides, close_amide_indices))
        amides_close= list(amides[close_amide_indices])
    
        # create dma file containing 2 close CONH groups
        xyz_close = xyz[amides_close[0]]; xyz_close = numpy.concatenate((xyz_close, xyz[amides_close[1]])) *libbbg.units.UNITS.AngstromToBohr
        a = libbbg.dma.DMA(nfrag=8)
        a.set_charges(2*self.__CONH_CHARGES)
        a.set_structure(pos=xyz_close, equal=True, atoms=','.join(2*['C','O','N','H']))

        if self.lprint: print " Saving DMA's for the 2 closest CONH subunits"
        self.__conhs_dma = a
    
        return amides_far, amides_close

    def _init_solefp(self, elect, polar, repul, disp, correc, 
                           ccut, pcut, ecut, mode):
        # initialize SolEFP calculator
        self.__SolEFP_Calculator = solvshift.solefp.EFP(\
                        elect = elect, pol   = polar, rep  = repul, disp = disp, corr = correc,
                        ccut  = ccut , pcut  = pcut , ecut = ecut , ea   = True,
                        freq  = True , cunit = True , mode = mode                             )
        return

    def _init_data(self):
        """Set up the intrinsic data structures"""
        # CONH fitting charges
        self.__CONH_CHARGES = [ 4.846770E-01,-5.688321E-01,-3.846247E-02, 1.192387E-01 ] # from nma-fit/4-cites-camm-dipole-1.62292.par
        # frequency shift terms to calculate
        #                  1         2         3         4         5 
        self.__dw_list = ('solcamm','ele_mea','ele_ea' ,'cor_mea','cor_ea',
        #                  6         7         8         9         10
                          'pol_mea','pol_ea' ,'rep_mea','rep_ea' ,'dis_mea',
        #                  11        12            13 
                          'dis_ea' ,'dis_mea_iso','total')
        #                      1              2             3              4             5          6
        self.__dw_env_list = ('ele_env_mea', 'ele_env_ea', 'cor_env_mea', 'cor_env_ea', 'tot_env', 'total+env')
        # initiate Report file
        self.__out = open(self.report, 'w')
        line = '%8s'                              % '#  Time'.ljust(8 )
        line+= ' %10s' * len(self.__dw_list)      % self.__dw_list
        line+= ' %10s' * len(self.__dw_env_list)  % self.__dw_env_list
        line+= ' %10s %10s %10s\n'                % ('RMS-C'.rjust(10), 'RMS-S'.rjust(10), 'RMS-AVE'.rjust(10))
        self.__out.write(line)
        self.__out.flush()

        return

    def _parse_solcamm(self, solcamm, mode):
        """Create SolCAMM object from SolEFP library or external DMA file"""
        if solcamm is None: return self.__solefp.get_property('solcamm', ea=True, mode=mode)
        else:               return libbbg.dma.DMA(solcamm)

    def _parse_res(self, res_file):
        """Reads the SolEFP-MD topology conversion file"""                                                                              
        res = dict()
        f = open(res_file); data = self._chop_comment(f.read(),'#'); f.close()
        t = re.compile('\[', re.DOTALL)
        sections = t.split(data)
        for sec in sections[1:]:
            frags = list()
            t = re.compile('@', re.DOTALL)
            sub_sections= t.split(sec)
            res_name = sub_sections[0].split()[0]
            for sub_sec in sub_sections[1:]:
                lines = sub_sec.split('\n')
                frag_name = lines[0].strip()
                is_reorder = is_supl = False
                for line in lines:
                    if 'reorder' in line: 
                        reorder = numpy.array(line.split()[1:], int)
                        is_reorder = True
                    if 'supl'    in line:    
                        supl = numpy.array(line.split()[1:], int)
                        is_supl = True
                if not is_reorder: 
                   print '\n ERROR: No reordering list (reorder) provided for fragment < %s > in residue < %s > !\n' % (frag_name, res_name)
                   sys.exit()
                if not is_supl   : 
                   print '\n ERROR: No superimposition list (supl) provided for fragment < %s > in residue < %s > !\n' % (frag_name, res_name)
                   sys.exit()
                frags.append( (frag_name, reorder, supl) )
            res[res_name] = frags
        return res
       
    # --- Input parsing/creation helper methods ---
    def _ap(self, res_data, res_name, atom_numbers):
        log = ''
        for frag in res_data:
            log += '$Frag\n'      
            log += '%10s %10s\n' % (res_name.ljust(10),frag[0].ljust(10))
            log += '%10s %50s\n' % ('reorder'.ljust(10), ' '.join(['%i'%x for x in frag[1]]).ljust(50))
            log += '%10s %50s\n' % ('supl'.ljust(10), ','.join(['%i'%x for x in frag[2]]).ljust(50))
            log += '%10s %50s    1\n' % ('atoms'.ljust(10), ','.join(['%i'%x for x in atom_numbers]).ljust(50))
            log += self._cp()
        return log

    def _app(self, res_data, res_name, atom_numbers):
        log = ''
        for frag in res_data:
            log += '$Frag\n'      
            log += '%10s %10s\n' % (res_name.ljust(10),frag[0].ljust(10))
            log += '%10s %50s\n' % ('reorder'.ljust(10), ' '.join(['%i'%x for x in frag[1]]).ljust(50))
            log += '%10s %50s\n' % ('supl'.ljust(10), ','.join(['%i'%x for x in frag[2]]).ljust(50))
            log += '%10s %50s    1\n' % ('atoms'.ljust(10), ','.join(['%i'%x for x in atom_numbers]).ljust(50))
            log += self._cp()
        return log
    
    def _qp(self, res_data, res_name, atom_numbers):
        log = ''
        for frag in res_data:
            log += '%10s %50s    1\n' % ('atoms'.ljust(10), ','.join(['%i'%x for x in atom_numbers]).ljust(50))
        return log
    
    def _cp(self): return '$endFrag\n\n'
    
    def _chop_comment(self, text, comment_mark='#'):
        new_text = []
        for line in text.split('\n'):
            if not line.strip().startswith(comment_mark):
               new_text.append(line)
        return '\n'.join(new_text)
