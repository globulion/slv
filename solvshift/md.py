# ---------------------------------------------------- #
#           SOLVATOCHROMIC FREQUENCY SHIFT             #
#           FOR MOLECULAR DYNAMICS MODULE              #
# ---------------------------------------------------- #

#from numpy     import *
#from units     import *
#from dma       import *
#from libbbg.utilities import *
#from solvshift.slvpar import Frag
#import sys, copy, pp, libbbg.utilities, re, \
#       units, numpy, time, dma, \
#  	      MDAnalysis.coordinates.xdrfile.libxdrfile2
#from MDAnalysis.coordinates.xdrfile.libxdrfile2 import xdrfile_open, xdrfile_close,\
#                                                       read_xtc_natoms, read_xtc  ,\
#                                                       read_xtc, DIM, exdrOK
#from MDAnalysis.coordinates.TRJ import TRJReader, NCDFReader

import numpy, libbbg, solvshift.slvpar, \
       sys, copy, pp, time, re, MDAnalysis.coordinates.xdrfile.libxdrfile2, \
       MDAnalysis.coordinates.TRJ
import os

DIM    = MDAnalysis.coordinates.xdrfile.libxdrfile2.DIM
exdrOK =  MDAnalysis.coordinates.xdrfile.libxdrfile2.exdrOK

__all__ = ['SLV_MD',]
__version__ = '3.0.3'

sys.stdout.flush()

class MDInput:
   """
 -------------------------------------------------------------------
                   Represents SLV MD input file
 -------------------------------------------------------------------
                                               April 23, 2014
                                               
 Import:
 
   from solvshift.md import MDInput
   
 Usage:
 
   inp  = MDInput(file)
   args, idx = inp.get()
 
 The <args> list and slicing indices <idx> can be now set to 
 solvshift.slvpar.Frag.set method:

   Frag_instance.set(frame[idx],*args)

 <args> is a list and contains the following information:

   args = [ ind, nmol, bsm, supl, reord ]

 For the details about these memorials read Frag class documentation
 by typing 

   from solvshift.slvpar import Frag
   help(Frag)

 Note: argument of MDInput class can be also a string containing 
       the contents of input file (instead of providing path to
       a file)

 -------------------------------------------------------------------
                              Last Revision:  9 May 2017 @Globulion
"""
   def __init__(self,file_name):
       self.__file_name = file_name if isinstance(file_name,str) else '<<from file>>'
       self.__bsm  = list()
       self.__nmol = list()
       self.__ind  = list()
       self.__supl = list()
       self.__reord= list()
       self.__frag_idx = 0
       self.__frame_slice = list()

       self.__call__(file_name)

   def get(self):
       """
Return the options for Frag.set method and slicing list

Usage:

    args, idx = MDInputInstance.get()
"""
       args = [ self.__ind, self.__nmol, self.__bsm, self.__supl, self.__reord ]
       idx  = self.__frame_slice
       return args, idx

   # auxiliary methods
   def get_idx(self):
       """get slicing list"""
       return self.__frame_slice
   
   def get_args(self):
       """get arguments for Frag.set method"""
       args = [ self.__ind, self.__nmol, self.__bsm, self.__supl, self.__reord ]
       return args

   def get_efp(self,frame):
       """parse EFP coordinates from MD frame"""
       return frame[self.__frame_slice]

   def __call__(self,file_name):
       if os.path.exists(file_name):
          fileo = open(file_name)
          text  = fileo.read()
          fileo.close()
       elif '\n' in file_name:
          text = file_name
       else:
          raise IOError, " The input file does not exists."
       self.__input_text = text
       templ = re.compile(r'\$Frag',re.DOTALL)
       sections = re.split(templ,text)[1:]
       for section in sections:
           self._read(section)
         
       return

   def _read(self,section):
       """read one unique (Sol)EFP fragment"""
       b_reorder = False
       b_supl    = False
       if 'reorder' in section: b_reorder = True
       if 'supl'    in section: b_supl    = True
       lines = section.split('\n')[1:]
       if lines[0].startswith('!'): return

       # read the fragment parameters
       frag = solvshift.slvpar.Frag(lines[0].split()[1])

       # reorder the fragment
       reord = None
       for line in lines:
           if line.startswith('!'):continue
           if 'reorder' in line:
               reord = numpy.array( map(int,line.split()[1:]), int ) - 1
       #if b_reorder:   # deprecated! Now reordering applies to MD structure
       #   frag.reorder(reord)

       # add superimpositon list
       supl = None
       for line in lines:
           if line.startswith('!'):continue
           if 'supl' in line:
               ll = line.split('supl')[-1]
               if ',' in ll: supl = libbbg.utilities.text_to_list(ll, delimiter=',') - 1
               else:         supl = libbbg.utilities.text_to_list(ll) - 1

       # parse atoms for (Sol)EFP fragment
       for line in lines:
           if line.startswith('!'):continue
           if 'atoms' in line:
               atoms, n_frags = line.split()[1:]
               if ',' in atoms: atoms   = libbbg.utilities.text_to_list(atoms, delimiter=',')
               else:            atoms   = libbbg.utilities.text_to_list(atoms)
               n_frags = int(n_frags)
               
               n_atoms = atoms[-1]-atoms[0]+1
               merror  = 'MDInputError: Invalid atomic indices or fragment numbers in fragment %i' % (self.__frag_idx + 1)
               merror += '\n -----> %s' % line
               assert n_atoms%n_frags==0, merror
               n_atoms_per_mol = (atoms[-1]-atoms[0]+1)/n_frags
               for i in range(n_frags):
                   self.__ind.append(self.__frag_idx)
                   self.__nmol.append(n_atoms_per_mol)
               atoms = list(numpy.array(atoms)-1)
               self.__frame_slice += atoms

       # --- in the case of no reordering, set the appropriate list
       if reord is None:
          reord = numpy.array( [ i for i in range(n_atoms_per_mol) ], int)

       # --- append to bsm and supl memorial lists
       self.__bsm.append(frag)
       self.__supl.append(supl)
       self.__reord.append(reord)
      
       # go to the next fragment
       self.__frag_idx += 1 
       return

   def __repr__(self):
       """print the input file"""
       N = 50
       log = '-'*N+'\n'
       log+= 'SLV-MD input file: '
       log+= self.__file_name+'\n\n'+self.__input_text
       log+= '-'*N+'\n'
       return str(log)


class layer(libbbg.units.UNITS):
    """represents single layer in SLV-MD environment"""
    def __init__(self,charges,nmol=1,):
        self.nmol = nmol
        self.charges = numpy.array(charges,dtype=numpy.float64)
        self.nat = len(charges)
        
    def set_frame(self,frame):
        """sets the coordinates for layer"""
        self.frame = frame
        return
    
    def get_shift(self,thr=None):
        """calculate frequency shifts"""
        pass

def read_inp(input):
    """parse charge input file for SLV-MD routines"""
    file = open(input)
    eprotein = True
    solvent = True
    ions = True
    # search for eprotein
    querry = '[ eprotein'
    line = file.readline()
    while True:
        if (not line.startswith('!') and querry in line): break
        elif ('[ ' in line and 'eprotein' not in line):
           eprotein = False; break
        line = file.readline()
    line = file.readline()

    solat = []; solch=[]; solchind=[]
    if eprotein:
       querry1 = 'solute atoms'
       querry2 = 'charges'
       while '[ ' not in line:
          if not line.startswith('!'):
             if querry1 in line:
                solat = numpy.array(line.split('=')[-1].split(','),dtype=int)
             elif querry2 in line:
                line = file.readline()
                while 'end charges' not in line:
                    if line.startswith('!'): 
                       line = file.readline()
                    else:
                       l = line.split()
                       solchind.append(l[0])
                       solch   .append(eval(l[1]))
                       line = file.readline()
          line = file.readline()
       solchind = numpy.array(solchind,dtype=int)
       solch    = numpy.array(solch   ,dtype=numpy.float64)
       
    # search for ions
    querry = '[ ions'
    while True:
        if (not line.startswith('!') and querry in line): break
        elif ('[ ' in line and 'ions' not in line):
           ions = False; break
        line = file.readline()
    line = file.readline()
    querry1 = 'name'
    querry2 = 'charges'
    querry3 = 'nmol'
    querry4 = 'threshold'
    iname = 'default ion'; ich=[]; inmol=0; ithr=60
    if ions:
       while '[ ' not in line:
          if not line.startswith('!'):
             if   querry1 in line:
                  iname = line.split('=')[1]
             elif querry2 in line:
                if '=' in line:
                  ich = numpy.array(line.split('=')[-1].split(','),dtype=numpy.float64)
                else:
                  line = file.readline()
                  while 'end charges' not in line:
                     if line.startswith('!'): 
                        line = file.readline()
                     else:
                        l = line.split()
                        ich   .append(eval(l[0]))
                        line = file.readline()
                  ich = numpy.array(ich,dtype=numpy.float64)
             elif querry3 in line:
                inmol = int(line.split('=')[1])
             elif querry4 in line:
                ithr = numpy.float64(line.split('=')[1])
          line = file.readline()

    # search for solvent
    querry = '[ solvent'
    while True:
        if (not line.startswith('!') and querry in line): break
        elif ('[ ' in line and 'solvent' not in line):
           solvent = False; break
        line = file.readline()
    line = file.readline()
    querry1 = 'name'
    querry2 = 'charges'
    querry3 = 'nmol'
    querry4 = 'threshold'
    sname = 'default solvent'; sch=[]; snmol=1; sthr=30
    if solvent:
       while '[ ' not in line:
          if not line.startswith('!'):
             if   querry1 in line:
                  sname = line.split('=')[1]
             elif querry2 in line:
                if '=' in line:
                  sch = numpy.array(line.split('=')[-1].split(','),dtype=numpy.float64)
                else:
                  line = file.readline()
                  while 'end charges' not in line:
                     if line.startswith('!'): 
                        line = file.readline()
                     else:
                        l = line.split()
                        sch   .append(eval(l[0]))
                        line = file.readline()
                  sch = numpy.array(sch,dtype=numpy.float64)
             elif querry3 in line:
                snmol = int(line.split('=')[1])
             elif querry4 in line:
                sthr = numpy.float64(line.split('=')[1])
          line = file.readline()
    #
    file.close()
    #print "solat, solchind, solch, ich, inmol, iname, ithr, sch, snmol, sname, sthr",solat, solchind, solch, ich, inmol, iname, ithr, sch, snmol, sname, sthr
    return solat, solchind, solch, ich, inmol, iname, ithr, sch, snmol, sname, sthr



def md_shifts_pp(frame,frame_idx,
                 prot_no,slt_no,ion_no,water_no,solute_atoms,solvent_nat,
                 epdma,idma,wdma_l,suplist,solute_parameters,
                 camm,reference_structure,
                 threshold,non_atomic,rings):
    """calculate shifts for parallel SLV-MD run"""
    ### [1] __updateDMA
    N = prot_no-slt_no
    M = prot_no
    K = M+ion_no
    # eprotein
    epframe = libbbg.utilities.choose(frame,solute_atoms)[:N]
    # ions
    iframe  = frame[M:M+ion_no]
    # waters
    wframe  = frame[M+ion_no:]
    #
    if prot_no-slt_no:
       epdma.set_structure(pos=epframe,equal=True)
       epdma.set_name('eprotein')
       #epdma.write('eprotein.xyz','xyz')
       
    if ion_no:
       idma.set_structure(pos= iframe,equal=True)
       idma.set_name('ions')
       #idma.write('ions.xyz','xyz')

    for mol in xrange(water_no):
        wdma_l[mol].set_structure(pos=frame[K+mol*solvent_nat:K+mol*solvent_nat+solvent_nat],equal=True)
        wdma_l[mol].set_name('solvent-%i'%(mol+1))
        #wdma_l[mol].write('solvent-%i.xyz'%(mol+1),'xyz')
        
    ### [2] __eval
    solute_pos = []
    for i in solute_atoms:
        solute_pos.append(frame[i])
    solute_pos = numpy.array(solute_pos)
    ### find solute's geometric center
    sltcnt = numpy.sum(solute_pos,axis=0)/slt_no
        
    ### add sites to solute
    if non_atomic:
       X = reference_structure[suplist]
       Y = solute_pos[suplist]
       sup = libbbg.utilities.SVDSuperimposer()
       sup.set(Y,X)
       sup.run()
       rms = sup.get_rms()
       rot, transl = sup.get_rotran()
       transformed = numpy.dot(reference_structure,rot) + transl
       
    ### add rings
    if rings:
       pass
    
    ### rotate the solute!
    if camm:
        ### control the rms of solute wrt reference (gas phase)
        rot, rms_ref  = libbbg.utilities.RotationMatrix(final=solute_pos[suplist],
                                                 initial=reference_structure[suplist])
        ### superimpose solute and parameters
        rot, rms_inst = libbbg.utilities.RotationMatrix(final=solute_pos[suplist],
                                                 initial=solute_parameters.pos[suplist])
            
        ### rotate parameters
        solute_parameters.MAKE_FULL()
        solute_parameters.Rotate(rot)

        ### update positions and origins of parameters
        if non_atomic: solute_parameters.set_structure( pos=transformed, equal=True )
        else: solute_parameters.set_structure( pos=solute_pos, equal=True )
    #solute_parameters.write('param.xyz','xyz')
        
    ### calculate frequency shifts!!!
    shift   = numpy.zeros(5,dtype=numpy.float64)
    shift_e = numpy.zeros(5,dtype=numpy.float64)
    shift_i = numpy.zeros(5,dtype=numpy.float64)
    shift_w = numpy.zeros(5,dtype=numpy.float64)

    # eprotein
    if prot_no-slt_no:
       shift_e+= libbbg.utilities.FrequencyShift(solute=solute_parameters,
                                          solvent=epdma,
                                          solute_structure=solute_parameters.get_origin())
    # ions
    if ion_no:
       shift_i+= libbbg.utilities.FrequencyShift(solute=solute_parameters,
                                          solvent=idma,
                                          solute_structure=solute_parameters.get_origin())
    # solvent molecules
    for mol in xrange(water_no):
        solcnt = numpy.sum(wdma_l[mol].get_pos(),axis=0)/numpy.float64(solvent_nat)
        R  = numpy.sqrt(numpy.sum((solcnt-sltcnt)**2))
        if R < threshold:
            s = libbbg.utilities.FrequencyShift(solute=solute_parameters,
                                         solvent=wdma_l[mol],
                                         solute_structure=solute_parameters.get_origin())
            #if R>50.: print " solvent-%4d  %36.4f %10.3f"%((mol+1),s[3],R)
            #if numpy.abs(s[3])>0.: wdma_l[mol].write('solvent-%d.xyz'%(mol+1),'xyz' )
            #if numpy.abs(s[3])>3.: print wdma_l[mol]
            shift_w+=s
    #print "TEST",water_no,ion_no,prot_no,slt_no
    shift = shift_e + shift_i + shift_w
    ### update the report
    if not camm: 
        rms_inst = 0
        rms_ref  = 0
    report_line   = "  %5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.4f %10.4f\n"\
                        % (frame_idx-1,shift[0],shift[1],shift[2],shift[3],shift[4],
                            rms_inst, rms_ref)
    report_line_e = "  %5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.4f %10.4f\n"\
                        % (frame_idx-1,shift_e[0],shift_e[1],shift_e[2],shift_e[3],shift_e[4],
                            rms_inst, rms_ref)
    report_line_i = "  %5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.4f %10.4f\n"\
                        % (frame_idx-1,shift_i[0],shift_i[1],shift_i[2],shift_i[3],shift_i[4],
                            rms_inst, rms_ref)
    report_line_w = "  %5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.4f %10.4f\n"\
                        % (frame_idx-1,shift_w[0],shift_w[1],shift_w[2],shift_w[3],shift_w[4],
                            rms_inst, rms_ref)

    return shift  , report_line  , \
           shift_e, report_line_e, \
           shift_i, report_line_i, \
           shift_w, report_line_w

class SLV_MD(libbbg.units.UNITS):
    """\
Represents MD-derived frequency shift distribution.
Calculates frequency shifts withdrawed from
MD simulation data providing:
- file with charges
- file with trajectories
- solute MCHO's parameters
- solute atomic indices

Usage: will be added soon!"""
    
    ### frame index, frame status and box dimensions
    frame_idx = 1
    status = exdrOK
    DIM = DIM

    
    def __init__(self,pkg="amber",charges="",trajectory="",
                 solute_parameters=None,nframes=16,
                 camm=False,suplist=[],ncpus=None,
                 non_atomic=False,inp=None,
                 natoms=None,nprotein=None,
                 report_name='report',):
        solat, solchind, solch, ich, inmol,\
        iname, ithr, sch, snmol, sname, sthr = read_inp(inp)
        self.solchind = solchind
        self.solch = solch
        self.wcharges = sch
        self.pkg = pkg
        self.charges = charges
        self.trajectory = trajectory
        self.solute_atoms = solat
        self.solvent_atno = len(sch)
        self.solute_parameters = solute_parameters
        self.threshold = sthr
        self.log = '\n'
        self.camm = camm
        self.reference_structure = solute_parameters.pos.copy()
        self.suplist = suplist
        self.slt_no  = len(solat)
        self.ion_no = inmol
        #
        self.__rings = None
        self.nprotein = nprotein
        self.natoms = natoms
        self.report_name = report_name
        #
        try:
            self.ion_charge = ich[0]
        except IndexError:
            self.ion_charge = None
            
        ### read total number of atoms
        if pkg=='gromacs': self.natoms = read_xtc_natoms(self.trajectory)
        else: self.natoms = natoms
        
        self.__non_atomic = non_atomic
        ### parallel computing processes
        if ncpus is not None:
           self.job_server = pp.Server(ncpus)
           self.jobs = []
        ### write action report on the disk
        self.__init_report()
        ### withdraw charges
        epc, ic, wc = self._Charges()
        ### initialize eprot, ions and solvent DMA distributions
        self.epdma, self.idma, self.wdma_l = self._CreateDMA(epc, ic, wc)
        ### proceed the frames
        to = time.time()
        if ncpus is None: 
           self._ProceedTheFrames_no_pp(nframes)
        else:
           self._ProceedTheFrames(nframes,ncpus=ncpus)
        print time.time()-to, "  :  TIME"
        ### report on average shift and std
        self.report  .write(self.__repr__())
        self.report_e.write(self._print(self.averages_e,self.stds_e))
        self.report_i.write(self._print(self.averages_i,self.stds_i))
        self.report_w.write(self._print(self.averages_w,self.stds_w))
        ### close the files
        self.report  .close()
        self.report_e.close()
        self.report_i.close()
        self.report_w.close()
        ### print cpu workers report on screen
        if ncpus is not None:
           self.job_server.print_stats()
    
    # Private methods

    def _CreateDMA(self, epcharges, icharges, wcharges):
        """initialize DMA objects appropriately"""
        # protein environment
        if (self.prot_no-self.slt_no):
           epdma = libbbg.dma.DMA(nfrag=self.prot_no-self.slt_no)
           epdma.set_moments(charges=epcharges)
        else:
           epdma=None
        # ions environment
        if self.ion_no:
           idma = libbbg.dma.DMA(nfrag=self.ion_no)
           idma.set_moments(charges=icharges)
        else:
           idma=None
        # solvent environment
        wdma_l = []
        for mol in xrange(self.water_no):
            dmai = libbbg.dma.DMA(nfrag=self.solvent_atno)
            dmai.set_moments(charges=wcharges)
            wdma_l.append(dmai)
 
        return epdma, idma, wdma_l
        
    def _Charges(self):
        """withdraw charges from MD files"""
        
        epcharges = [] # protein environment
        icharges  = [] # ions
        wcharges  = [] # solvent molecules
        charges_file = open(self.charges)
        
        if self.pkg.lower()=="amber":
           querry = "CHARGE"
           line = charges_file.readline()
           while querry not in line:
                 line = charges_file.readline()

           line = charges_file.readline()
           line = charges_file.readline()

           while not line.startswith("%FLAG MASS"):
                 epcharges+= [ numpy.float64(x) for x in line.split() ]
                 line = charges_file.readline()
                 
           #self.prot_no = len(epcharges)
           self.prot_no = self.nprotein
           # add the ion and solvent charges
           for i in range(self.ion_no): 
               icharges.append(self.ion_charge)
           self.water_no = (self.natoms - self.prot_no - self.ion_no) / self.solvent_atno
           wcharges = self.wcharges
           
           ### convert charges to au units!
           epcharges = numpy.array(epcharges)[:self.nprotein]/18.2223
           epcharges = list(epcharges)
           
        elif self.pkg.lower()=="gromacs":     
           querry = "[ atoms ]"
           line=charges_file.readline()

           while 1:
             if querry in line: break
             else: line=charges_file.readline()

           line=charges_file.readline()
           line=charges_file.readline()
           while line!='\n':
             if not line.startswith(';'):
                epcharges.append(line.split()[6])
             line=charges_file.readline()
           
           self.prot_no = len(epcharges)
           # add the ion and solvent charges
           for i in range(self.ion_no): 
               icharges.append(self.ion_charge)
               
           self.water_no = (self.natoms - self.prot_no - self.ion_no) / self.solvent_atno
           wcharges = self.wcharges
        
        ### substitute the charges by new ones provided in the input file
        for i,index in enumerate(self.solchind):
            epcharges[index] = self.solch[i]
            
        ### remove the charges for a probe from a list!
        k = 0
        for atom in self.solute_atoms:
            epcharges.pop(atom-k)
            k+= 1
            
        epcharges = numpy.array(epcharges,dtype=numpy.float64)
        icharges  = numpy.array( icharges,dtype=numpy.float64)
        wcharges  = numpy.array( wcharges,dtype=numpy.float64)
        
        ### warn if the sum of eprotein charges is not integer
        sum_chg = numpy.sum(epcharges)
        a_chg   = numpy.abs(sum_chg)
        d_chg   = numpy.abs(a_chg-int(round(a_chg)))
        if d_chg>0.0001:
            print "\n WARNING! The sum of charges for eprotein is not integer! sum=%5.6f\n" % sum_chg
            print " Quitting... check carefully the charges and restart your task\n" % sum_chg
            raise ValueError
        else:
            print "\n The sum of charges for eprotein is %5.6f\n" % sum_chg
        return epcharges, icharges, wcharges

    def _ProceedTheFrames(self,nframes,ncpus):
        """proceeds frame by frame to collect frequency shifts"""
        
        print "\n SLV FREQUENCY SHIFT DISTRIBUTION CALCULATION MODE\n"
        self.frequency_shifts   = []
        self.frequency_shifts_e = []
        self.frequency_shifts_i = []
        self.frequency_shifts_w = []
        self.shift_corrections  = []
                
        if self.pkg == 'amber':
           # read frames
           md = MDAnalysis.coordinates.TRJ.TRJReader(self.trajectory,self.natoms)
           md.open_trajectory()

           for i in range(nframes):
               print " * Reading frame %10i from total %i frames"%(i+1,nframes)
               frame = numpy.array(md.next()) * self.AngstromToBohr

               group = "group-%i"%(i/ncpus)
               
               self.jobs.append( self.job_server.submit(func=md_shifts_pp,
                    args=(frame,self.frame_idx,self.prot_no,self.slt_no,
                            self.ion_no,self.water_no,self.solute_atoms,self.solvent_atno,
                            self.epdma,self.idma,self.wdma_l,self.suplist,
                            self.solute_parameters,self.camm,
                            self.reference_structure,self.threshold,self.__non_atomic,
                            self.__rings),
                    depfuncs=(),
                    group=group,
                    modules=("libbbg.utilities","numpy","libbbg.units","libbbg.dma",
                             "MDAnalysis.coordinates.xdrfile.libxdrfile2","clemtp"),
                           ) )
               if (i%ncpus==0 and i!=0): self.job_server.wait()
               i+=1
               #
               self.frame_idx += 1
               
        elif self.pkg == 'gromacs':
             # initialize the data
             natoms = MDAnalysis.coordinates.xdrfile.libxdrfile2.read_xtc_natoms(self.trajectory)
             frame = numpy.zeros((natoms,DIM),dtype=numpy.float32)
             box = numpy.zeros((DIM,DIM),dtype=numpy.float32)
             XTC = MDAnalysis.coordinates.xdrfile.libxdrfile2.xdrfile_open(self.trajectory,'r')
             # read frames
             i=0
             #while self.status == exdrOK:
             for i in range(nframes):
                 group = "group-%i"%(i/ncpus)
                 self.__read_xtc(XTC,natoms,frame,box,DIM)
                 self.jobs.append( self.job_server.submit(func=md_shifts_pp,
                      args=(frame,self.frame_idx,self.prot_no,self.slt_no,
                            self.ion_no,self.water_no,self.solute_atoms,self.solvent_atno,
                            self.epdma,self.idma,self.wdma_l,self.suplist,
                            self.solute_parameters,self.camm,
                            self.reference_structure,self.threshold,self.__non_atomic,
                            self.__rings),
                      depfuncs=(),
                      group=group,
                      modules=("libbbg.utilities","numpy","libbbg.units","libbbg.dma",
                               "MDAnalysis.coordinates.xdrfile.libxdrfile2","clemtp"),
                           ) )
                 if (i%ncpus==0 and i!=0): self.job_server.wait()
                 i+=1
             
             # close the trajectory file
             MDAnalysis.coordinates.xdrfile.libxdrfile2.xdrfile_close(XTC)
                    
        # gather the results
        for job in self.jobs:
            shift, report_line, shift_e, report_line_e, shift_i, report_line_i, shift_w, report_line_w = job()
            #
            for II,JJ,KK,LL in [(self.report  ,  report_line  , self.frequency_shifts   , shift  ),
                                (self.report_e,  report_line_e, self.frequency_shifts_e , shift_e), 
                                (self.report_i,  report_line_i, self.frequency_shifts_i , shift_i),
                                (self.report_w,  report_line_w, self.frequency_shifts_w , shift_w)]:
               II.write(JJ)
               II.flush()
               KK.append(LL)

        ### finalize the results
        self.frequency_shifts   = numpy.array( self.frequency_shifts   )
        self.frequency_shifts_e = numpy.array( self.frequency_shifts_e )
        self.frequency_shifts_i = numpy.array( self.frequency_shifts_i )
        self.frequency_shifts_w = numpy.array( self.frequency_shifts_w )

        self.averages  = numpy.array([numpy.average(self.frequency_shifts[:,0]),
                                      numpy.average(self.frequency_shifts[:,1]),  
                                      numpy.average(self.frequency_shifts[:,2]),
                                      numpy.average(self.frequency_shifts[:,3]),
                                      numpy.average(self.frequency_shifts[:,4])])

        self.averages_e= numpy.array([numpy.average(self.frequency_shifts_e[:,0]),
                                      numpy.average(self.frequency_shifts_e[:,1]),  
                                      numpy.average(self.frequency_shifts_e[:,2]),
                                      numpy.average(self.frequency_shifts_e[:,3]),
                                      numpy.average(self.frequency_shifts_e[:,4])])

        self.averages_i= numpy.array([numpy.average(self.frequency_shifts_i[:,0]),
                                      numpy.average(self.frequency_shifts_i[:,1]),  
                                      numpy.average(self.frequency_shifts_i[:,2]),
                                      numpy.average(self.frequency_shifts_i[:,3]),
                                      numpy.average(self.frequency_shifts_i[:,4])])

        self.averages_w= numpy.array([numpy.average(self.frequency_shifts_w[:,0]),
                                      numpy.average(self.frequency_shifts_w[:,1]),
                                      numpy.average(self.frequency_shifts_w[:,2]),
                                      numpy.average(self.frequency_shifts_w[:,3]),
                                      numpy.average(self.frequency_shifts_w[:,4])])

        self.stds  = numpy.array([numpy.std(self.frequency_shifts[:,0]),
                                  numpy.std(self.frequency_shifts[:,1]),  
                                  numpy.std(self.frequency_shifts[:,2]),
                                  numpy.std(self.frequency_shifts[:,3]),
                                  numpy.std(self.frequency_shifts[:,4])])

        self.stds_e= numpy.array([numpy.std(self.frequency_shifts_e[:,0]),
                                  numpy.std(self.frequency_shifts_e[:,1]),  
                                  numpy.std(self.frequency_shifts_e[:,2]),
                                  numpy.std(self.frequency_shifts_e[:,3]),
                                  numpy.std(self.frequency_shifts_e[:,4])])

        self.stds_i= numpy.array([numpy.std(self.frequency_shifts_i[:,0]),
                                  numpy.std(self.frequency_shifts_i[:,1]),  
                                  numpy.std(self.frequency_shifts_i[:,2]),
                                  numpy.std(self.frequency_shifts_i[:,3]),
                                  numpy.std(self.frequency_shifts_i[:,4])])

        self.stds_w= numpy.array([numpy.std(self.frequency_shifts_w[:,0]),
                                  numpy.std(self.frequency_shifts_w[:,1]),  
                                  numpy.std(self.frequency_shifts_w[:,2]),
                                  numpy.std(self.frequency_shifts_w[:,3]),
                                  numpy.std(self.frequency_shifts_w[:,4])])

    def _ProceedTheFrames_no_pp(self,nframes):
        """proceeds frame by frame to collect frequency shifts"""
        
        print "\n SLV FREQUENCY SHIFT DISTRIBUTION CALCULATION MODE\n"
        self.frequency_shifts = []
        self.shift_corrections= []
        
        if self.pkg == 'amber':           
           # read frames
           traj = open(self.trajectory)
           line = traj.readline()
           nframes = int(line.split("=")[-1])
           line = traj.readline()

           for j in range(nframes):
               print " * Reading frame %10i from total %i frames"%(j+1,nframes)
               frame = []
               while line:
                     if len(line)==25:
                        line = traj.readline()
                        break
                     else:
                        try:frame+= [ line[8*n:8*n+8] for n in range(10) ]
                        except IndexError: continue
                        except ValueError: pass
                        line = traj.readline()
               frame = frame[:-2]
               #print len(frame),frame[-3:]
               ### reshape and convert to au units!
               frame = numpy.array(frame,dtype=numpy.float64).reshape( len(frame)/3, 3 ) * self.AngstromToBohr
            
               # evaluate frequency shift for the current frame and save
               self.__eval(frame)
               
        elif self.pkg == 'gromacs':
             # initialize the data
             natoms = MDAnalysis.coordinates.xdrfile.libxdrfile2.read_xtc_natoms(self.trajectory)
             frame = numpy.zeros((natoms,DIM),dtype=numpy.float32)
             box = numpy.zeros((DIM,DIM),dtype=numpy.float32)
             XTC = MDAnalysis.coordinates.xdrfile.libxdrfile2.xdrfile_open(self.trajectory,'r')
             # read frames
             for i in range(nframes):
             #while self.status == exdrOK:
                   self.__read_xtc(XTC,natoms,frame,box,DIM)
                   
                   # evaluate frequency shift for the current frame and save
                   self.__eval(frame)
                   
             # close the trajectory file      
             MDAnalysis.coordinates.xdrfile.libxdrfile2.xdrfile_close(XTC)
             
        ### finalize the results
        self.frequency_shifts = numpy.array( self.frequency_shifts )
        self.averages  = numpy.array([numpy.average(self.frequency_shifts[:,0]),
                                      numpy.average(self.frequency_shifts[:,1]),  
                                      numpy.average(self.frequency_shifts[:,2]),
                                      numpy.average(self.frequency_shifts[:,3]),
                                      numpy.average(self.frequency_shifts[:,4])])
        self.stds = numpy.array([numpy.std(self.frequency_shifts[:,0]),
                                 numpy.std(self.frequency_shifts[:,1]),  
                                 numpy.std(self.frequency_shifts[:,2]),
                                 numpy.std(self.frequency_shifts[:,3]),
                                 numpy.std(self.frequency_shifts[:,4])])
                           
    def __read_dcd(self):
        """reads the DCD trajectory"""
        pass
    
    def __read_xtc(self,XTC,natoms,frame,box,DIM=DIM):
        """reads the XTC trajectory from XTC file object"""
        self.status, step, time, prec = MDAnalysis.coordinates.xdrfile.libxdrfile2.read_xtc(XTC,box,frame)
        frame *= self.NanometerToBohr
        centre = frame.mean(axis=0)
        print " * Reading frame %10i"%(self.frame_idx)
        self.frame_idx+=1
        return
    
    def __updateDMA(self,frame):
        """divide frame for various layers and puts them into DMAs"""
        N = self.prot_no-self.slt_no
        M = self.prot_no
        K = M+self.ion_no
        # eprotein
        epframe = libbbg.utilities.choose(frame,self.solute_atoms)[:N]
        # ions
        iframe  = frame[M:M+self.ion_no]
        # waters
        wframe  = frame[M+self.ion_no:]
        #
        self.epdma.set_structure(pos=epframe,equal=True)
        self. idma.set_structure(pos= iframe,equal=True)
        self.epdma.set_name('eprotein')
        self. idma.set_name('ions')
        for mol in xrange(self.water_no):
            o = mol*3 + 0
            h1= mol*3 + 1
            h2= mol*3 + 2
            self.wdma_l[mol].set_structure(pos=frame[K+mol*3:K+mol*3+3],equal=True)
            self.wdma_l[mol].set_name('solvent-%i'%(mol+1))
            #
        return
        
        
    def __eval(self,frame):
        """evaluates frequency shifts for one MD frame"""
        ### update DMA environment layers position
        self.__updateDMA(frame)

        ### update solute position
        solute_pos = []
        for i in self.solute_atoms:
            solute_pos.append(frame[i])
        solute_pos = numpy.array(solute_pos)
        #for i in range(self.prot_no-self.slt_no):
        #    if self.epdma.get_pos()[i] in solute_pos: print "YYY", self.epdma.get_pos()[i]
        #print 
        #print solute_pos
        ### find solute's geometric center
        sltcnt = numpy.sum(solute_pos,axis=0)/self.slt_no
        
        ### rotate the solute!
        if self.camm:
            ### control the rms of solute wrt reference (gas phase)
            rot, rms_ref  = libbbg.utilities.RotationMatrix(final=solute_pos[self.suplist],
                                                     initial=self.reference_structure[self.suplist])
            ### superimpose solute and parameters
            rot, rms_inst = libbbg.utilities.RotationMatrix(final=solute_pos[self.suplist],
                                                     initial=self.solute_parameters.pos[self.suplist])
            
            ### rotate parameters
            self.solute_parameters.MAKE_FULL()
            self.solute_parameters.Rotate(rot)

            ### update positions and origins of parameters
            self.solute_parameters.set_structure( pos=solute_pos, equal=True )
        
        #sys.exit()       
        ### calculate frequency shifts!!!
        # eprotein
        shift = libbbg.utilities.FrequencyShift(solute=self.solute_parameters,
                                         solvent=self.epdma,
                                         solute_structure=solute_pos,
                                         threshold=1e10)
        # ions
        shift+= libbbg.utilities.FrequencyShift(solute=self.solute_parameters,
                                         solvent=self.idma,
                                         solute_structure=solute_pos,
                                         threshold=1e10)
        # solvent molecules
        for mol in xrange(self.water_no):
            solcnt = numpy.sum(self.wdma_l[mol].get_pos(),axis=0)/3.
            R  = numpy.sqrt(numpy.sum((solcnt-sltcnt)**2))
            if R < self.threshold:
               shift+= libbbg.utilities.FrequencyShift(solute=self.solute_parameters,
                                                solvent=self.wdma_l[mol],
                                                solute_structure=solute_pos,
                                                threshold=1e10)
        ### save the results
        self.frequency_shifts.append( shift )
           
        ### update the report
        if not self.camm: 
            rms_inst = 0
            rms_ref  = 0
        self.report.write( "%5i %10.3f %10.3f %10.3f %10.3f %10.3f %10.4f %10.4f\n"\
                            % (self.frame_idx-1,shift[0],shift[1],shift[2],shift[3],shift[4],
                                rms_inst, rms_ref) )
        self.report.flush()
        return
        
    # Public methods
    
    def get_shifts(self):
        """returns results as arrays of: 
           - frequency shifts for each particular frame
           - averages of frequency shifts
           - standard deviations of frequency shifts
        """
        return (self.frequency_shifts, self.averages, self.stds )

    def __init_report(self):
        """initialize the report"""
        self.report = open(self.report_name+'_T.md','w')
        self.report.write("# SLV Molecular Dynamics Shifts - TOTAL\n")
        self.report.write("# Units: shifts in [cm-1], rms in [Bohr]\n")
        self.report.write("# %5s %10s %10s %10s %10s %10s %10s %10s\n"\
                                       %("Frame"    .rjust(5),
                                         "1"        .rjust(10),
                                         "1+2"      .rjust(10),
                                         "1+2+3"    .rjust(10),
                                         "1+2+3+4"  .rjust(10),
                                         "1+2+3+4+5".rjust(10),
                                         "rms-inst" .rjust(10),
                                         "rms-ref"  .rjust(10)))
        self.report.flush()

        # eprotein
        self.report_e = open(self.report_name+'_e.md','w')
        self.report_e.write("# SLV Molecular Dynamics Shifts due to EPROTEIN\n")
        self.report_e.write("# Units: shifts in [cm-1], rms in [Bohr]\n")
        self.report_e.write("# %5s %10s %10s %10s %10s %10s %10s %10s\n"\
                                       %("Frame"    .rjust(5),
                                         "1"        .rjust(10),
                                         "1+2"      .rjust(10),
                                         "1+2+3"    .rjust(10),
                                         "1+2+3+4"  .rjust(10),
                                         "1+2+3+4+5".rjust(10),
                                         "rms-inst" .rjust(10),
                                         "rms-ref"  .rjust(10)))
        self.report_e.flush()

        # ions
        self.report_i = open(self.report_name+'_i.md','w')
        self.report_i.write("# SLV Molecular Dynamics Shifts due to IONS\n")
        self.report_i.write("# Units: shifts in [cm-1], rms in [Bohr]\n")
        self.report_i.write("# %5s %10s %10s %10s %10s %10s %10s %10s\n"\
                                       %("Frame"    .rjust(5),
                                         "1"        .rjust(10),
                                         "1+2"      .rjust(10),
                                         "1+2+3"    .rjust(10),
                                         "1+2+3+4"  .rjust(10),
                                         "1+2+3+4+5".rjust(10),
                                         "rms-inst" .rjust(10),
                                         "rms-ref"  .rjust(10)))
        self.report_i.flush()

        # eprotein
        self.report_w = open(self.report_name+'_w.md','w')
        self.report_w.write("# SLV Molecular Dynamics Shifts due to SOLVENT\n")
        self.report_w.write("# Units: shifts in [cm-1], rms in [Bohr]\n")
        self.report_w.write("# %5s %10s %10s %10s %10s %10s %10s %10s\n"\
                                       %("Frame"    .rjust(5),
                                         "1"        .rjust(10),
                                         "1+2"      .rjust(10),
                                         "1+2+3"    .rjust(10),
                                         "1+2+3+4"  .rjust(10),
                                         "1+2+3+4+5".rjust(10),
                                         "rms-inst" .rjust(10),
                                         "rms-ref"  .rjust(10)))
        self.report_w.flush()
        return

    def _print(self, av, st):
        """print happy me!"""
        
        log = "\n      Your dreamed frequency shift distribution! [cm-1]\n"
        log+=   "      -------------------------------------------------\n"
        log+=   "         1         :  %16.1f ± %.1f\n" % ( av[0],
                                                            st[0] )
        log+=   "         1+2       :  %16.1f ± %.1f\n" % ( av[1],
                                                            st[1] )
        log+=   "         1+2+3     :  %16.1f ± %.1f\n" % ( av[2],
                                                            st[2] )
        log+=   "         1+2+3+4   :  %16.1f ± %.1f\n" % ( av[3],
                                                            st[3] )
        log+=   "         1+2+3+4+5 :  %16.1f ± %.1f\n" % ( av[4],
                                                            st[4] )
        log+=   "      -------------------------------------------------\n"
        log+= "\n\n"
          
        return str(log)                       
                       
    def __repr__(self):
        """print me!"""
        
        self.log+= "\n  Your dreamed frequency shift distribution! [cm-1]\n"
        self.log+=   "  -------------------------------------------------\n"
        self.log+=   "     1         :  %16.1f ± %.1f\n" % ( self.averages[0],
                                                             self.stds[0] )
        self.log+=   "     1+2       :  %16.1f ± %.1f\n" % ( self.averages[1],
                                                             self.stds[1] )
        self.log+=   "     1+2+3     :  %16.1f ± %.1f\n" % ( self.averages[2],
                                                             self.stds[2] )
        self.log+=   "     1+2+3+4   :  %16.1f ± %.1f\n" % ( self.averages[3],
                                                             self.stds[3] )
        self.log+=   "     1+2+3+4+5 :  %16.1f ± %.1f\n" % ( self.averages[4],
                                                             self.stds[4] )
        self.log+=   "  -------------------------------------------------\n"
        self.log+= "\n\n"
          
        return str(self.log)
