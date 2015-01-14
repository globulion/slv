#!/usr/bin/python
"""
 Allign the molecule 
"""
from libbbg.utilities import QMFile, ParseDmatFromFchk, ParseVecFromFchk, PotentialContourMap
from libbbg.dma import DMA
from numpy import array
import time 

dma = DMA('mescn.dma')
dma_pot = DMA('mescn-camm.dma')
atoms = [1,2,3,4]
bonds = [[1,2],  [1,3] , [3,4]]
atid  = (3,1,4)
radii = array([1.5, 1.4, 1.8, 1.5, 2, 2, 2])
levels = [-0.03, -0.018,-0.013,-0.01,-0.009,-0.008,-0.007,-0.006,-0.005,-0.004,-0.003,
          -0.002,-0.001,-0.0005, -0.0002 -0.0001, -0.00005,0.0,0.00005, 0.0001, 0.0002, 
          0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 
                   0.01, 0.013, 0.018, 0.03]
#levels = levels[3:-3]
#levels=None
mol = QMFile('mescn.xyz', mol=True, basis='6-311++G**')
#dmat = mol.get_dmat()
bfs = mol.get_bfs()
mol = mol.get_mol()
print mol
dmat = ParseDmatFromFchk('mescn_alligned.fchk', 121)
vec = ParseVecFromFchk('mescn.fchk')[:19]
#vec = None
dmat= None
radii=None
basis='6-311++G**'
atnos = [6.0, 7.0, 16.0, 6.0, 1.0 ,1.0, 1.0]
atom_colors=['sienna','darkseagreen','orange','sienna']

print " Generating the map..."
t0 = time.time()
map = PotentialContourMap(dma, atoms, bonds, atid, allign_axes=(1,2,0),
          pad_x=10.0, pad_y=10.0, dx=0.3, dy=0.3, levels=levels, 
          #pad_left=100.0, pad_right=100.0, 
          #pad_up=100.0, pad_down=100.0,
          dma_pot=dma_pot, qm_mask_thr=0.15,  linthresh=0.0001, 
          radii=radii, vec=vec, dmat=dmat, basis=basis, atnos=atnos, block=True, name=None,
          atom_colors=atom_colors, bond_width=4, bond_fmt='k--')
map.make()
t1 = time.time()
print "\n It took %10.4f seconds\n" % (t1-t0)
