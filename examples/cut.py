#!/usr/bin/python
"""
  --- SLV-MD Cut-off estimator ---

reads the water simulation trajectory
and cuts spheres of water calculating
frequency shifts for that placing there
MeSCN molecule in the center of sphere
for various spheres of increasing radii

Usage:

cut.py [xtc] [xyz solute] [anh] [solpar]
"""
# ---------------------------------------------------------------------------------
import sys, gaussfreq, utilities, dma
import numpy as np
from MDAnalysis.coordinates.xdrfile.libxdrfile import xdrfile_open, xdrfile_close,\
                                                      read_xtc_natoms, read_xtc  ,\
                                                      read_xtc, DIM, exdrOK
from solvshift.slv import SLV
from units import UNITS
# ---------------------------------------------------------------------------------

print __doc__
if len(sys.argv) < 2: sys.exit()


# --- input parameters
xtcfile = sys.argv[1]
n       = 199 # no of frame form MD to take to build spheres
solute  = sys.argv[2]
solute_atoms = ['C','N','S','C','H','H','H'] # MeSCN
anh_file= sys.argv[3]
parameters= sys.argv[4]
mode_id = 3 # SCN vibration

# --- utilities
def frame_n(file,n):
    """read nth frame"""
    natoms = read_xtc_natoms(file)
    frame = np.zeros((natoms,DIM),dtype=np.float32)
    box = np.zeros((DIM,DIM),dtype=np.float32)
    XTC = xdrfile_open(file,'r')
    for i in range(n):
          status, step, time, prec = read_xtc(XTC,box,frame)

    xdrfile_close(XTC)
    return natoms, frame

def calc_center(xyz):
    """return center of xyz structure"""
    n=len(xyz)
    return np.sum(xyz,axis=0)/np.float64(n)

def write_waters(sphere,center,name='waters'):
    """writes waters """
    out = open(name,'w')
    n = len(sphere) + 1
    print >> out, "%i\n"%n
    print >> out,"%s"%'X',
    print >> out, 3*"%10.6f "%tuple(center)
    for mol in range((n-1)/3):
        o = mol*3+0
        h1= mol*3+1
        h2= mol*3+2
        print >> out, "%s"%'O',
        print >> out, 3*"%10.6f "%tuple(sphere[o ])
        print >> out, "%s"%'H',
        print >> out, 3*"%10.6f "%tuple(sphere[h1])
        print >> out, "%s"%'H',
        print >> out, 3*"%10.6f "%tuple(sphere[h2])


    out.close()
    return

def get_dma(waters):
    """prepare Charge DMA object for waters"""
    n = len(waters)
    soldma= dma.DMA(nfrag=n)
    soldma.set_structure(pos=waters,equal=True)
    #
    charges = np.zeros(n,dtype=np.float64)
    for mol in xrange(n/3):
        o = mol*3+0
        h1= mol*3+1
        h2= mol*3+2
        charges[o ] = -0.784994
        charges[h1] =  0.392497
        charges[h2] =  0.392497
    #charges[-1]=0;charges[7]=0 # --- check what happens if we reset some charges...
    
    soldma.set_moments(charges=charges)
    return soldma

# [1] SOLUTE
# --- extract coordinates of solute
solute_xyz = utilities.Read_xyz_file(solute,ar=1)[1]
sltcnt = calc_center(solute_xyz)

# --- SolDMA for solute
slt_dma = utilities.ParseDMA(parameters,'c')

# --- rotate the parameters
slt_init_pos = slt_dma.get_pos()
slt_dma.MAKE_FULL()
rot,rms = utilities.RotationMatrix(initial=slt_init_pos,final=solute_xyz)
slt_dma.Rotate(rot)
print "    RMS of solute rotation: ", rms
slt_dma.set_structure(pos=solute_xyz,equal=True)

# --- anharmonic data for frequency shift
anh = gaussfreq.FREQ(anh_file)

# [2] SPHERES
# --- read the n-th frame
natoms, frame = frame_n(xtcfile,n) 
frame *= UNITS.NanometerToBohr * UNITS.BohrToAngstrom
center = np.sum(frame,axis=0) / natoms
write_waters(frame,center)

# --- translate solute to center of geometry of the spheres
solute_xyz += center-sltcnt

# --- iterate over all water molecules
#     frame, natoms and center objects are static 
def cut_shell(r_min,r_max):
    """cut shell of water molecues from the center of frame"""
    shell = []
    for mol in xrange(natoms/3):
        o = mol*3+0
        h1= mol*3+1
        h2= mol*3+2
        R1= np.sqrt( np.sum((frame[o ]-center)**2) )
        R2= np.sqrt( np.sum((frame[h1]-center)**2) )
        R3= np.sqrt( np.sum((frame[h2]-center)**2) )
        if (r_min<R1<r_max):
         if (r_min<R2<r_max):
          if (r_min<R3<r_max):
           shell.append( frame[o ] )
           shell.append( frame[h1] )
           shell.append( frame[h2] )
    #
    return np.array(shell,dtype=np.float64)

# --- iterate over sphere radii
radii = np.array([[ 5, 8],# in Angstroem !
                  [ 6, 9],
                  [ 7,10], 
                  [ 8,11],
                  [ 9,13],
                  [10,14],
                  [11,16],
                  [12,17], 
                  [13,18],
                  [14,19],
                  [15,20],
                  [16,22],
                  [17,23],
                  [18,24],
                  [ 4,30]])
print "\n FREQUENCY SHIFTS [cm-1] \n"
for radius in radii:
    sphere = cut_shell(radius[0],radius[1])

    # --- print the structures of water spheres for check
    print " * Writing for sphere %02i-%02i"%tuple(radius),
    write_waters(sphere,center,name='sphere-%02i-%02i'%tuple(radius))
    
    # --- obtain DMA distributions for water shells
    #slt_in_sphere = concatenate([solute_xyz,sphere])
    sol_dma = get_dma(sphere)

    # --- calculate frequency shift!
    shift = SLV(nsatoms=['O','H','H'],                    
                nstatoms=solute_atoms,
                L=anh.L,
                mode_id=mode_id,
                solute=slt_dma,
                solute_structure=slt_dma.get_pos(),
                solvent=sol_dma,
                mixed=True,
                camm=slt_dma,
                redmass=anh.redmass,
                freq=anh.freq,
                #bsm=sol_dma,
                lprint=False,
                ref_structure=slt_dma.get_pos() )
    print 4*"%10.4f"%tuple(shift.shift[0][:4])
