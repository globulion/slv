#!/usr/bin/python
"""
Calculates EDS frequency shifts

usage:

./eds.py [anh] [eds dirs] [mode] [method]

note:
eds dirs type for example:
w20:sder/20 - parse w20 as fdir
              and w20/sder/20 as sdir
mode        - number of mode (Helico, normal numbers)
method      - either HF or MP2
"""
# -----------------------------------
print __doc__
from sys import argv, exit
if len(argv)==1: exit()
import os
os.environ['__IMPORT__COULOMB__']='1'
from numpy import *
from libbbg.utilities import *
from solvshift.mcho import MCHO
from solvshift.diff import DIFF
from libbbg.gaussfreq import FREQ
from libbbg.units import UNITS
# -----------------------------------
anh    = argv[1]
eds    = argv[2]
mode   = int(argv[3])
method = argv[4]
ANH = FREQ(anh)
nmodes = len(ANH.freq)
CALC=DIFF(freq=ANH.freq,
          dir='.',
          L=ANH.L,
          eds=eds,
          eds_method=method)

os.system('export SLV_MODE_ID=%d' % (mode-1))
#PRINT(CALC.fEDS)
t = sqrt(ANH.redmass*UNITS.AmuToElectronMass)
print CALC.fEDS.shape
f_el = CALC.fEDS[:,1]
f_mtp= CALC.fEDS[:,2]
f_rep= CALC.fEDS[:,4]
f_pol= CALC.fEDS[:,5]
print f_el  * t
print f_mtp * t

PARAM = MCHO(gijj=ANH.K3,freq=ANH.freq,
             redmass=ANH.redmass,
             MM_mode_id=mode-1,
             eds=True,
             fEDS=CALC.fEDS,
             sEDS=CALC.sEDS,
             nmodes=nmodes)
shifts = PARAM.eds_shifts

log = "\n"
log+= " E(EL,10)  %20.2f %10.2f %10.2f\n"%tuple(shifts[:,1])
log+= "    E(EL,MTP) %17.2f %10.2f %10.2f\n"%tuple(shifts[:,2])
log+= "    E(EL,PEN) %17.2f %10.2f %10.2f\n"%tuple(shifts[:,3])
log+= " E(EX,HL)  %20.2f %10.2f %10.2f\n"%tuple(shifts[:,4])
log+= " E(DEL)    %20.2f %10.2f %10.2f\n"  %tuple(shifts[:,5])
log+= " DE(HF)    %20.2f %10.2f %10.2f\n\n" %tuple(shifts[:,6])

if method.lower()=='mp2':
   log+= 'happy mp2\n'
   log+= " E(MP,2)   %20.2f %10.2f %10.2f\n"%tuple(shifts[:,7])
   log+= " E(EL,R,12)%20.2f %10.2f %10.2f\n"%tuple(shifts[:,8])
   log+= "    E(EL,M,2)%17.2f %10.2f %10.2f\n"%tuple(shifts[:,9])
   log+= "    E(EL,P,2)%17.2f %10.2f %10.2f\n"%tuple(shifts[:,10])
   log+= " E(DS,20)  %20.2f %10.2f %10.2f\n"%tuple(shifts[:,11])
   log+= " DE(EX-DEL,2)%18.2f %10.2f %10.2f\n"  %tuple(shifts[:,12])
   log+= " DE(MP2)   %20.2f %10.2f %10.2f\n\n" %tuple(shifts[:,13])


print log
