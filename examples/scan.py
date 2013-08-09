#!/usr/bin/python

"""
make scans of 1 water molecule 
placed in various distences from
solute molecule apart

Usage:

scans.py [xyz_solute] [xyz_solvent] 
         [parameters] [bsm file]
         [anharmonic file]
"""

from random import random
from sys import argv
import utilities, math, gaussfreq, sys
from solvshift.slv import SLV
import numpy as np
from units import *

print __doc__
if len(argv) < 5: sys.exit()

# --- input parameters
solute = argv[1]
solvent= argv[2]
parameters= argv[3]
bsm= argv[4]
anh_file= argv[5]

n_directions = 40
n_pass = 65 # integer!; in bohr
solsol = ['C','N','S','C','H','H','H','O','H','H']
mode_id= 3

# --- utilities
def calc_center(xyz):
    """return center of xyz structure"""
    n=len(xyz)
    return np.sum(xyz,axis=0)/np.float64(n)

def norm(vec):
    """normalize vector"""
    n=np.sqrt(sum(vec**2))
    return vec/n

def write_xyz(coord,atoms,name):
    """write xyz file"""
    n=len(coord)
    out=open(name,'w')
    out.write(str(n)+"\n\n")
    for i in range(n):
        print >> out, "%s %12.6f %12.6f %12.6f"%(atoms[i],coord[i,0],coord[i,1],coord[i,2])
    out.close()
    return

# --- extract coordinates of solute and solvent
solute_xyz = utilities.Read_xyz_file(solute,ar=1)[1]
solvent_xyz = utilities.Read_xyz_file(solvent,ar=1)[1]

# --- compute center of mass of solute and solvent
sltcnt = calc_center(solute_xyz)
solcnt = calc_center(solvent_xyz)

# --- translate solvent to center of geometry of solute
solvent_xyz -= solcnt-sltcnt 

#print solvent_xyz*UNITS.BohrToAngstrom
#print solute_xyz*UNITS.BohrToAngstrom

# --- for each direction calculate the position of
#     water molecule and calculate shifts; write 
#     shifts to a nice report
slt_dma = utilities.ParseDMA(parameters,'c')
sol_dma = utilities.ParseDMA(bsm,'c')
# --- rotate the parameters
slt_init_pos = slt_dma.get_pos()
slt_dma.MAKE_FULL()
rot,rms = utilities.RotationMatrix(initial=slt_init_pos,final=solute_xyz)
slt_dma.Rotate(rot)
print rms
slt_dma.set_structure(pos=solute_xyz,equal=True)

sol_init_pos = sol_dma.get_pos()
sol_dma.MAKE_FULL()
rot,rms = utilities.RotationMatrix(initial=sol_init_pos,final=solvent_xyz)
sol_dma.Rotate(rot)
print rms
sol_dma.set_structure(pos=solvent_xyz,equal=True)

#
anh = gaussfreq.FREQ(anh_file)

for ndir in range(n_directions):
    dirctn = norm( np.array([random()*2.-1. for i in range(3)],dtype=np.float64) )
    out_shift1 = open('shiftscan-'+str(ndir+1)+'.out','w')
    out_shift2 = open('shiftscan+'+str(ndir+1)+'.out','w') 
    for nps in range(n_pass):
        solnew = solvent_xyz - dirctn * (nps+1+7)
        sol_dma.set_structure(pos=solnew,equal=True)
        shift = SLV(nsatoms=solsol[7:],
                    nstatoms=solsol[:7],
                    L=anh.L,
                    mode_id=mode_id,
                    solute=slt_dma,
                    solute_structure=slt_dma.get_pos(),
                    solvent=sol_dma,
                    mixed=False,
                    camm=slt_dma,
                    redmass=anh.redmass,
                    freq=anh.freq,
                    bsm=sol_dma,
                    lprint=False,
                    ref_structure=slt_dma.get_pos() )
        out_shift1.write("%10.2f"%((nps+8)*UNITS.BohrToAngstrom))
        out_shift1.write("%12.4f %12.4f %12.4f %12.4f\n"% tuple(shift.shift[0][:4]))

        solnew = np.concatenate([solute_xyz,solnew]) * UNITS.BohrToAngstrom
        write_xyz(solnew,solsol,name='./str/scan-%i-%i' % (ndir,nps))
        #
        solnew = solvent_xyz + dirctn * (nps+1+7)
        sol_dma.set_structure(pos=solnew,equal=True)
        shift = SLV(nsatoms=solsol[7:],
                    nstatoms=solsol[:7],
                    L=anh.L,
                    mode_id=mode_id,
                    solute=slt_dma,
                    solute_structure=slt_dma.get_pos(),
                    solvent=sol_dma,
                    mixed=False,
                    camm=slt_dma,
                    redmass=anh.redmass,
                    freq=anh.freq,
                    bsm=sol_dma,
                    lprint=False,
                    ref_structure=slt_dma.get_pos() )
        out_shift2.write("%10.2f"%((nps+8)*UNITS.BohrToAngstrom))
        out_shift2.write("%12.4f %12.4f %12.4f %12.4f\n"% tuple(shift.shift[0][:4]))


        solnew = np.concatenate([solute_xyz,solnew]) * UNITS.BohrToAngstrom
        write_xyz(solnew,solsol,name='./str/scan-%i+%i' % (ndir,nps))


    out_shift1.close()
    out_shift2.close()
