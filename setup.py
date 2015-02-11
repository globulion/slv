#!/usr/bin/python
# -*- coding: utf-8 -*-

# -----------------------------------------------
from numpy.distutils.core import setup, Extension
from glob import glob
import sys
# -----------------------------------------------

### extension module specifications

# exchange-repulsion EFP
SHFTEX = Extension(name='solvshift.shftex',
                   sources=['solvshift/src/shftex.f'],)
# exchange-repulsion + CT EFP
SHFTCE = Extension(name='solvshift.shftce',
                   sources=['solvshift/src/shftce.f'],)
# exchange-repulsion interaction energy only
EXREP  = Extension(name='solvshift.exrep',
                   sources=['solvshift/src/exrep.f'],)

# rotation of wave function in EFP
EFPROT = Extension(name='solvshift.efprot',
                   sources=['solvshift/src/efprot.f'],)


# data files: BSM fragment files
molecules = ['water', 'water2', 'dmso'   , 'meoh' , 'chcl3',
             'dcm'  , 'na+'   , 'me-so3-', 'so3--', 'meoac',
             'nma'  , 'nma-d7', 'mescn'  , 'li+'  , 'et2coo', ] 

numerical = { 'mescn': ('num_0.006', 'num_0.025') , }

frg_files = list()

for mol in molecules:
    files = ('solvshift-dat/frg/%s' % mol , [ i for i in glob('frg/%s/*.frg' % mol) ] +
                                            [ i for i in glob('frg/%s/*.xyz' % mol) ] +
                                            [ i for i in glob('frg/%s/*.efp' % mol) ] )

    frg_files.append(files)

    if mol in numerical.keys():
       dirs = numerical[mol]
       for dir in dirs:
           files = ('solvshift-dat/frg/%s/%s' % (mol,dir), [ i for i in glob('frg/%s/%s/*.frg' % (mol, dir)) ] )
           frg_files.append(files)
       

### install

setup(name='SOLVSHIFT',
      version='1.0.1',
      description='Solvatochromic shifts from coarse-grained SolEFP theory',
      author='Bartosz Błasiak',
      author_email='blasiak.bartosz@gmail.com',
      url='no-page-yet',
      packages=['solvshift',
                #'solvshift.util',
                #'solvshift.diff',
               ],
      ext_modules=[SHFTEX,EFPROT,SHFTCE,EXREP,],
      data_files=frg_files,
     )
