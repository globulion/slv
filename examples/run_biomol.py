#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
 SolEFP/EFP2 fragmentation model of biomolecule. 
 Development script.

 Usage: [pdb (MD snapshot)]

"""
from sys import argv
if len(argv)==1: 
   print __doc__
   exit()

import solvshift.biomol
import sys

if __name__=='__main__':
     res    = 'gmx.res'
     probe  = 'probe.res'
     mode   =  4
     traj   = sys.argv[1]
     solefp = solvshift.biomol.BiomoleculeFragmentation(\
                       res, probe, mode, solcamm=None, lprint=True, ccut=45.0, pcut=16.0, ecut=13.0,
                       repul=True,
                       write_solefp_input=True, write_debug_file=True,
                       rcl_algorithm='remove_by_name', include_ions=True, include_polar=True)

     solefp.run(traj, traj, dframes=1, nframes=1, conh='nma', out_inp='biomolecule.sol')
