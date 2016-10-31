#!/usr/bin/python
"""
 Creates a scan ...
"""
from solvshift.slvpar import Frag
from solvshift.solefp import EFP
from libbbg.utilities import QMFile
import numpy, libbbg
from numpy import cos, sin

def cog(xyz):
    c = xyz.sum(axis=0)/len(xyz)
    return c

# 
xyz_nma = QMFile('nma.xyz').get_pos()
xyz_wat = QMFile('water.xyz').get_pos()

cog_nma = cog(xyz_nma)
cog_wat = cog(xyz_wat)


xyz_wat-= cog_wat
xyz_nma-= cog_nma

r = xyz_nma[7] - xyz_nma[1]
vec_uni = r/numpy.linalg.norm(r)

# translate water molecule
xyz_wat += 6.0 * vec_uni

xyz = numpy.concatenate( [xyz_nma, xyz_wat] )

print xyz * libbbg.units.UNITS.BohrToAngstrom

# construct SolEFP solver
mode = 8
ccut, pcut, ecut = 10e10, 10e10, 10e10
supl_1 = numpy.array([2,8,3,5]) - 1
supl_2 = numpy.array([1,2,3]) - 1
reord_1= numpy.array([1,2,3,4,5,6,7,8,9,10,11,12]) - 1
reord_2= numpy.array([1,2,3]) - 1

frg_solute = Frag('nma')
frg_solvent= Frag('water')


nmol = [  frg_solute.get_natoms() , frg_solvent.get_natoms() ]
ind  = [  0                       , 1 ]
bsm  = ( frg_solute, frg_solvent )
supl = [ supl_1, supl_2 ]
reord= [ reord_1, reord_2 ]

efp = EFP(elect=1,pol=1,disp=1,rep=1,corr=1,freq=True,mode=mode,
          ccut=ccut,pcut=pcut,ecut=ecut,
          cunit=True)

efp.set(xyz,ind,nmol,bsm,supl,reord)
efp.eval(0)
shifts = efp.get_shifts()

results = [ [ shifts['total'], shifts['ele_tot'] ] ]

# make scan
t = 0.1
rot_mat = numpy.array([cos(t), -sin(t), 0, sin(t), cos(t), 0, 0, 0, 1]).reshape(3,3)

for i in range(20):
    print " Making scan %d" % (i+1)
    xyz_wat = numpy.dot(xyz_wat, rot_mat)
    frg_solvent.sup(xyz_wat)
    xyz = numpy.concatenate( [xyz_nma, xyz_wat] )
    bsm  = ( frg_solute, frg_solvent )

    efp.set(xyz,ind,nmol,bsm,supl,reord)
    efp.eval(0)
    shifts = efp.get_shifts()
    results.append([ shifts['total'], shifts['ele_tot'] ])

results = numpy.array(results)

# write result to the file
out = open('dat','w')
for i in range(len(results)):
    out.write('%14.2f %14.2f\n' % (results[i,0], results[i,1]) )
out.close()
print results
print efp.get_rms_max()
