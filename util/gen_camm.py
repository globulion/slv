#!/usr/bin/python
#*-* coding: UTF-8 *-*
from libbbg.utilities import Read_xyz_file, reorder, ParseVecFromFchk, ParseDmatFromFchk, get_pmloca, order
import PyQuante.Ints, coulomb.multip, glob, os

def bua(file_fchk,basis,bonds,vec,vec_ref,method,mult,charge):
    """helper temporary function from solvshift.diff"""
    molecule = Read_xyz_file(file_fchk,mol=True,
                             mult=mult,charge=charge,
                             name='happy dummy molecule')
    
    bfs        = PyQuante.Ints.getbasis(molecule,basis)
    basis_size = len(bfs)
    #print " - basis size= ", basis_size
    print " - parsing %s density matrix" % method
    dmat = ParseDmatFromFchk(file_fchk,basis_size,method)
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
       vec = ParseVecFromFchk(file_fchk)[:nae,:]
       
       print " - Pipek-Mezey localization of %i orbitals..." %nae
       tran, vec = get_pmloca(natoms,mapi=bfs.LIST1,sao=SAO,
                                        vecin=vec,nae=nae,
                                        maxit=100000,conv=1.0E-10,
                                        lprint=False,
                                        freeze=None)
       vec, sim = order(vec_ref,vec,start=0)
       print sim
       check_sim(sim)
    ### calculate CAMMs
    print " - multipole integrals in AO basis evaluation..."
    camm = coulomb.multip.MULTIP(molecule=molecule,
                                 basis=basis,
                                 method='b3lyp',
                                 matrix=dmat,
                                 transition=False,
                                 bonds=bonds,vec=vec)
    print " - calculation of %s"%camm.operation
    camm.camms()
    #camm.mmms()
    #camm.__printMMMs__()
    #CAMM.__printCAMMs__()
    
    dma = camm.get()[0]
    dma.write(file_fchk[:-5]+'.camm')
    print " --- Writing file:  :", file_fchk[:-5]+'.camm'    
    return

def gen_camm(file=None, basis='6-311++G**',bonds=None,ncpus=4,
             vec=None,vec_ref=None,natoms=7,method='SCF',mult=1,charge=0,
             fchk_search='*_.fchk', fchk_dir=os.environ['PWD']): 
    """
Temporary utility:
Calculates CAMMs/CBAMMs/LMTPs from density matrix from GAUSSIAN09.
The DMA file with the appropriate distribution is saved.


Usage: 
 1) calculating CAMM: gen_camm(file=..., basis=..., method=..., fchk_search='*_.fchk', fchk_dir='./')

    where 
       o file is a FCHK file. If it is set to None (default) then it is assumed
         you have several FCHK files in the directory specified in fchk_dir. 
         In this case FCHK files will be searched according to fchk_search wildcard.
         The wildcard can be set to any bash-type wildcard.
       o basis is a string defining basis set name compatible with PyQuante
       o method refers to density level from FCHK (SCF, MP2, CC etc)

 2) calculating CBAMM: the same as above but add bonds (see in the code of this function)
 3) calculating LMTP: the same as above but specify vec (LCAO-MO coefficients for canonical set of orbitals)
    and natoms (number of atoms). 

Wanring: Parallel run (for many fchk files) is not working at present 
because of problem with administration of some shared libraries.

Notes:
Note that for the purpose of differentiation the option with file=None is handy
since the FCHK files are first always sorted. Then the first FCHK is assumed to be the reference
(not perturbed) file. The option vec_ref is not implemented yet but it is trivial.
"""
    if vec_ref is not None: raise NotImplementedError, 'vec_ref is not supported yet! The first FCHK file (after sorting) is assumed to be the reference!'
    if file is None:
       pliki_fchk  = glob.glob(fchk_dir+'/'+fchk_search)
       pliki_fchk.sort()
       print "\n Kolejność plików. Sprawdź czy się zgadzają!\n"  
       for i in range(len(pliki_fchk)):
           print pliki_fchk[i]
       print
    else:
       pliki_fchk = [file]
    
    # compute reference vectors
    if vec is not None:
       ref_mol = Read_xyz_file(pliki_fchk[0],mol=True,
                                         mult=mult,charge=charge,
                                         name='happy dummy molecule')
       bfs_ref    = PyQuante.Ints.getbasis(ref_mol,basis)
       basis_size = len(bfs_ref)
       sao_ref    = PyQuante.Ints.getS(bfs_ref)
       print " - basis size= ", basis_size
       nae = vec
       print " - nae       = ", len(vec)
       vec_ref = ParseVecFromFchk(pliki_fchk[0])[:nae,:]
       t, vec_ref = get_pmloca(natoms,mapi=bfs_ref.LIST1,
                               sao=sao_ref,
                               vecin=vec_ref,nae=nae,
                               maxit=100000,conv=1.0E-19,
                               lprint=False,
                               freeze=None)
    # submit the jobs!
    for file_fchk in pliki_fchk:
        bua (file_fchk,basis,bonds,vec,vec_ref,method,mult,charge)
    print
    return

if __name__ == '__main__':
   print gen_camm.__doc__
   from sys import argv
   fchk = argv[1]
   basis= argv[2]
   method=argv[3]
   # tests:
   #gen_camm('t.fchk', basis='6-311++G**', method='SCF')
   #gen_camm(None    , basis='6-31G*'    , method='SCF',fchk_search='*9.fchk')
   gen_camm(fchk, basis, method=method, fchk_search='*.fchk')
   exit()
