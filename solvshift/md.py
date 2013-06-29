# ---------------------------------------------------- #
#           SOLVATOCHROMIC FREQUENCY SHIFT             #
#           FOR MOLECULAR DYNAMICS MODULE              #
# ---------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
from MDAnalysis.coordinates.xdrfile.libxdrfile import xdrfile_open, xdrfile_close,\
                                                      read_xtc_natoms, read_xtc  ,\
                                                      read_xtc, DIM, exdrOK
                                                      
sys.stdout.flush()

class SLV_MD(UNITS):
    """represents MD-derived frequency shift distribution.
    Calculates frequency shifts withdrawed from
    MD simulation data providing:
    - file with charges
    - file with trajectories
    - solute MCHO's parameters
    - solute atomic indices
    """
    
    ### frame index, frame status and box dimensions
    frame_idx = 1
    status = exdrOK
    DIM = DIM
    
    def __init__(self,pkg="amber",charges="",trajectory="",
                 solute_atoms=(),solute_parameters=0,threshold=20,
                 camm=False,suplist=[]):

        self.pkg = pkg
        self.charges = charges
        self.trajectory = trajectory
        self.solute_atoms = solute_atoms
        self.solute_parameters = solute_parameters
        self.threshold = threshold
        self.log = '\n'
        self.camm = camm
        self.reference_structure = solute_parameters.pos.copy()
        self.suplist = suplist
        
        # write action report on the disk
        self.report = open("report.md",'w')
        self.report.write("SLV Molecular Dynamics Shifts\n")
        self.report.write("Units: shifts in [cm-1], rms in [Bohr]\n")
        self.report.write("%5s %10s %10s %10s %10s %10s %10s %10s\n"\
                                       %("Frame"    .rjust(5),
                                         "1"        .rjust(10),
                                         "1+2"      .rjust(10),
                                         "1+2+3"    .rjust(10),
                                         "1+2+3+4"  .rjust(10),
                                         "1+2+3+4+5".rjust(10),
                                         "rms-inst" .rjust(10),
                                         "rms-ref"  .rjust(10)))
        self.report.flush()
        # withdraw charges
        self._Charges()
        # proceed the frames
        from time import time
        to = time()
        self._ProceedTheFrames()
        print time()-to, "  :  TIME"
        # report on average shift and std
        self.report.write(self.__repr__())
        # close the files
        self.report.close()
    
    # Private methods
    
    def _Charges(self,ion_no=4,ion_charge=1):
        """withdraw charges from MD files"""
        
        charges = []   
        charges_file = open(self.charges)
        
        if self.pkg.lower()=="amber":
           querry = "CHARGE"
           line = charges_file.readline()
           while querry not in line:
                 line = charges_file.readline()

           line = charges_file.readline()
           line = charges_file.readline()

           while not line.startswith("%FLAG MASS"):
                 charges+= [ float64(x) for x in line.split() ]
                 line = charges_file.readline()

           ### convert charges to au units!
           charges = array(charges)/18.2223
           
        elif self.pkg.lower()=="gromacs":     
           querry = "[ atoms ]"
           line=charges_file.readline()

           while 1:
             if querry in line: break
             else: line=charges_file.readline()

           line=charges_file.readline()
           line=charges_file.readline()
           while line!='\n':
             charges.append(line.split()[6])
             line=charges_file.readline()
           
           prot_no = len(charges)
           natoms = read_xtc_natoms(self.trajectory)           
           # add the ion and solvent charges
           for i in range(ion_no): charges.append(ion_charge)
           water_no = (natoms - prot_no - ion_no) / 3
           for i in range(water_no):
               charges.append(-0.784994)
               charges.append( 0.392497)
               charges.append( 0.392497)
        
        ### zero-out the charges for a probe!       
        for atom in self.solute_atoms:
            charges[atom] = 0
        #for atom in [196,197]:#,194,195,204,205]:
        #    charges[atom] = 0 
        charges = array(charges,dtype=float64)
        
        ### construct DMA charges for solvent environment
        solvent = DMA(nfrag=len(charges))
        solvent.DMA[0] = array(charges)
       
        ### remember the solvent DMA object
        self.solvent = solvent
        print charges
        return

    def _ProceedTheFrames(self):
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
               frame = array(frame,dtype=float64).reshape( len(frame)/3, 3 ) * self.AngstromToBohr
            
               # evaluate frequency shift for the current frame and save
               self.__eval(frame)
               
        elif self.pkg == 'gromacs':
             # initialize the data
             natoms = read_xtc_natoms(self.trajectory)
             frame = zeros((natoms,DIM),dtype=float32)
             box = zeros((DIM,DIM),dtype=float32)
             XTC = xdrfile_open(self.trajectory,'r')
             # read frames
             while self.status == exdrOK:
                   self.__read_xtc(XTC,natoms,frame,box,DIM)
                   
                   # evaluate frequency shift for the current frame and save
                   self.__eval(frame)
                   
                   xdrfile_close(XTC)
                   break
                   
             # close the trajectory file      
             #xdrfile_close(XTC)
             
        ### finalize the results
        self.frequency_shifts = array( self.frequency_shifts )
        self.averages  = array([average(self.frequency_shifts[:,0]),
                                average(self.frequency_shifts[:,1]),
                                average(self.frequency_shifts[:,2]),
                                average(self.frequency_shifts[:,3]),
                                average(self.frequency_shifts[:,4])])
        self.stds = array([std(self.frequency_shifts[:,0]),
                           std(self.frequency_shifts[:,1]),
                           std(self.frequency_shifts[:,2]),
                           std(self.frequency_shifts[:,3]),
                           std(self.frequency_shifts[:,4])])

    def __read_dcd(self):
        """reads the DCD trajectory"""
        pass
    
    def __read_xtc(self,XTC,natoms,frame,box,DIM=DIM):
        """reads the XTC trajectory from XTC file object"""
        self.status, step, time, prec = read_xtc(XTC,box,frame)
        frame *= self.NanometerToBohr
        centre = frame.mean(axis=0)
        print " * Reading frame %10i"%(self.frame_idx)
        self.frame_idx+=1
        return
        
    def __eval(self,frame):
        """evaluates frequency shifts for one MD frame"""
        ### update DMA solvent position
        self.solvent.pos = frame

        ### update solute position
        solute_pos = []
        for i in self.solute_atoms:
            solute_pos.append(frame[i])
        solute_pos = array(solute_pos)
        
        ### rotate the solute!
        if self.camm:
            ### control the rms of solute wrt reference (gas phase)
            rot, rms_ref  = RotationMatrix(final=solute_pos[self.suplist],
                                        initial=self.reference_structure[self.suplist])
            ### superimpose solute and parameters
            rot, rms_inst = RotationMatrix(final=solute_pos[self.suplist],
                                        initial=self.solute_parameters.pos[self.suplist])
            ### update positions and origins of parameters
            self.solute_parameters.set_structure( pos=solute_pos, equal=True )
               
            ### rotate parameters
            self.solute_parameters.MAKE_FULL()
            self.solute_parameters.Rotate(rot)
               
        ### calculate frequency shifts!!!
        shift = FrequencyShift(solute=self.solute_parameters,
                               solvent=self.solvent,
                               solute_structure=solute_pos,
                               threshold=self.threshold)
                                  
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
                
    def __repr__(self):
        """print me!"""
        
        self.log+= "\n Your dreamed frequency shift distribution! [cm-1]\n"
        self.log+=   " -------------------------------------------------\n"
        self.log+=   "    1         :  %16.1f ± %.1f\n" % ( self.averages[0],
                                                            self.stds[0] )
        self.log+=   "    1+2       :  %16.1f ± %.1f\n" % ( self.averages[1],
                                                            self.stds[1] )
        self.log+=   "    1+2+3     :  %16.1f ± %.1f\n" % ( self.averages[2],
                                                            self.stds[2] )
        self.log+=   "    1+2+3+4   :  %16.1f ± %.1f\n" % ( self.averages[3],
                                                            self.stds[3] )
        self.log+=   "    1+2+3+4+5 :  %16.1f ± %.1f\n" % ( self.averages[4],
                                                            self.stds[4] )
        self.log+=   " -------------------------------------------------\n"
        self.log+= "\n\n"
          
        return str(self.log)