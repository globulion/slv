# ---------------------------------------------------- #
#           SOLVATOCHROMIC FREQUENCY SHIFT             #
#           FOR MOLECULAR DYNAMICS MODULE              #
# ---------------------------------------------------- #

from numpy     import *
from units     import *
from dma       import *
from utilities import *
import sys, copy
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
    
    def __init__(self,pkg="amber",charges="",trajectory="",
                 solute_atoms=(),solute_parameters=0,threshold=20,
                 camm=False):

        self.pkg = pkg
        self.charges = charges
        self.trajectory = trajectory
        self.solute_atoms = solute_atoms
        self.solute_parameters = solute_parameters
        self.threshold = threshold
        self.log = '\n'
        self.camm = camm
        self.reference_structure = solute_parameters.pos.copy()
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
    
    def _Charges(self):
        """withdraw charges from MD files"""
    
        if self.pkg.lower()=="amber":
           charges_file = open(self.charges)
           querry = "CHARGE"
           line = charges_file.readline()
           while querry not in line:
                 line = charges_file.readline()

           line = charges_file.readline()
           line = charges_file.readline()

           charges = []
           while not line.startswith("%FLAG MASS"):
                 charges+= [ float64(x) for x in line.split() ]
                 line = charges_file.readline()

           ### convert charges to au units!
           charges = array(charges)/18.2223
           for atom in self.solute_atoms:
               charges[atom] = 0
           #zerowalista = arange(2415-1,2486-1+1)
           #zerowalista = arange(1-1,12-1+1)
           #for atom in zerowalista:
           #    data[atom] = 0

           ### construct DMA charges for solvent environment
           solvent = DMA(nfrag=len(charges))
           solvent.DMA[0] = array(charges)
       
           ### remember the solvent DMA object
           self.solvent = solvent
           charges_file.close()
       
       
    def _ProceedTheFrames(self):
        """proceeds frame by frame to collect frequency shifts"""
        
        print "\n SLV FREQUENCY SHIFT DISTRIBUTION CALCULATION MODE\n"
        self.frequency_shifts = []
        self.shift_corrections= []
           
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
               rot, rms_ref  = RotationMatrix(final=solute_pos[:5],
                                            initial=self.reference_structure[:5])
               ### superimpose solute and parameters
               rot, rms_inst = RotationMatrix(final=solute_pos[:5],
                                            initial=self.solute_parameters.pos[:5])
               ### update positions and origins of parameters
               self.solute_parameters.pos = array(solute_pos)
               self.solute_parameters.origin = array(solute_pos)
               
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
                               % ((j+1),shift[0],shift[1],shift[2],shift[3],shift[4],
                                  rms_inst, rms_ref) )
            self.report.flush()

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