#!/usr/bin/env python
"""
Computes the CAMM parameters from all FCHK files in the working directory.

Usage:

  slv_gen-camm [method] [basis] 

Methods for reading one-particle density matrix: 
   SCF - SCF wavefunction
   MP2 - Moller-Plesset
   CC  - Coupled cluster
basis set has to be fully compatible with Gaussian (turn of shells
and all coefficients!)
"""
from sys  import argv
if len(argv)==1: 
   print __doc__
   exit()

import libbbg.utilities
import PyQuante.Ints
import coulomb.multip
import glob

def bua(file_fchk,basis,bonds,vec,vec_ref,method):
    molecule = libbbg.utilities.Read_xyz_file(file_fchk,mol=True,
                                       mult=1,charge=0,
                                       name='happy dummy molecule')
    
    bfs        = PyQuante.Ints.getbasis(molecule,basis)
    basis_size = len(bfs)
    #print " - basis size= ", basis_size
    print " - parsing %s density matrix" % method
    dmat = libbbg.utilities.ParseDmatFromFchk(file_fchk,basis_size,method)
    def check_sim(l):
        """check the sim list"""
        for x,y in l:
            i=0;j=0
            for a,b in l:
                if a==x: i+=1
                if b==y: j+=1
            if (i>1 or j>1): 
                print " --- !ERROR! --- "
                break
                                                                              
    ### parse vectors and make Pipek-Mezey transformation
    if vec is not None:
       natoms= len(molecule.atoms)
       SAO   = PyQuante.Ints.getS(bfs)
       print " - ovelrap AO matrix evaluation..."
       nae = vec
       vec = libbbg.utilities.ParseVecFromFchk(file_fchk)[:nae,:]
       
       print " - Pipek-Mezey localization of %i orbitals..." %nae
       tran, vec = libbbg.utilities.get_pmloca(natoms,mapi=bfs.LIST1,sao=SAO,
                                        vecin=vec,nae=nae,
                                        maxit=100000,conv=1.0E-10,
                                        lprint=False,
                                        freeze=None)
       vec, sim = libbbg.utilities.order(vec_ref,vec,start=0)
       print sim
       check_sim(sim)
    ### calculate CAMMs
    print " - multipole integrals in AO basis evaluation..."
    camm = coulomb.multip.MULTIP(molecule=molecule,
                                 basis=basis,
                                 method='b3lyp',
                                 matrix=dmat,
                                 transition=False,
                                 bonds=bonds,vec=vec,hexadecapoles=False)
    print " - calculation of %s"%camm.operation
    camm.camms()
    #camm.mmms()
    #camm.__printMMMs__()
    #CAMM.__printCAMMs__()
    
    dma = camm.get()[0]
    dma.write(file_fchk[:-5]+'.camm')
    print " --- Writing file:  :", file_fchk[:-5]+'.camm'    
    return

def CalculateCAMM_not_pp(basis='6-311++G**',bonds=[],ncpus=4,
                   vec=None,vec_ref=None,method='SCF'): 
    """calculates CAMMs from density matrix from GAUSSIAN09 using COULOMB.py routines"""
    pliki_fchk  = glob.glob('./*_.fchk')
    pliki_fchk.sort()
    print "\n FCHK Files!\n"  
    for i in range(len(pliki_fchk)):
        print pliki_fchk[i]
    print
    
    # compute reference vectors
    if vec is not None:
       ref_mol = libbbg.utilities.Read_xyz_file(pliki_fchk[0],mol=True,
                                         mult=1,charge=0,
                                         name='happy dummy molecule')
       natoms     = len(ref_mol.get_pos())
       bfs_ref    = PyQuante.Ints.getbasis(ref_mol,basis)
       basis_size = len(bfs_ref)
       sao_ref    = PyQuante.Ints.getS(bfs_ref)
       print " - basis size= ", basis_size
       nae = vec
       print " - nae       = ", len(vec)
       vec_ref = libbbg.utilities.ParseVecFromFchk(pliki_fchk[0])[:nae,:]
       t, vec_ref = libbbg.utilities.get_pmloca(natoms,mapi=bfs_ref.LIST1,
                                         sao=sao_ref,
                                         vecin=vec_ref,nae=nae,
                                         maxit=100000,conv=1.0E-19,
                                         lprint=False,
                                         freeze=None)
    # submit the jobs!
    for file_fchk in pliki_fchk:
        bua (file_fchk,basis,bonds,vec,vec_ref,method)
    print
    return

if __name__=='__main__':
   #bonds = [map(int,(x.split(','))) for x in argv[-1].split('-')]
   #bonds=[(1,0),(2,1),(3,2),(4,2)] # NMA
   bonds = None
   vec = None
   method = argv[1]
   basis  = argv[2]
   CalculateCAMM_not_pp(basis=basis,bonds=bonds,vec=vec,method=method)

