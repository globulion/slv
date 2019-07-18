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
molecules = ['water'            , 'water2'           , 'dmso'                , 'meoh'      ,     'chcl3'        ,     
             'dcm'              , 'na+'              , 'me-so3-'             , 'so3--'     ,     'meoac'        ,
             'nma'              , 'nma-d7'           , 'mescn'               , 'li+'       ,     'et2coo'       , 
             'mecn'             , '4-me-imidazole'   , '4-me-phenol'         , 'ccl4'      ,     'etoh'         ,
             'chonh2'           , 'chonhme'          , 'comenh2'             , 'ethane'    ,     'mecoo-'       ,
             'menh3+'           , 'methane'          , 'n-propane'           , 'me-guanidinium+'                ,
             'thf'              , 'dmf'              , 'cyclohexane'         , 'benzene'   ,     'mesh'         ,
             'menh2'            , 'phenol'           , 'cho-ch-nh2-ch3'      , 'dms'       ,     'imidazole'    ,
             'me-s-cho'         , 'me-s-co-ch-ch2'   , 'mecooh'              , 'phenolate-',     'imidazolium+' , 
             'cl-'              , 'me-s-13c-15n'     ,                                                           ] 

numerical = { 'mescn': ('num_0.006', 'num_0.025') , 
              'nma'  : ('num_0.006',            ) , }

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

# data files: MD-EFP2 topology conversion files
tc_files  = list()
tc_files.append( ('solvshift-dat/dat/md'       , [ 'dat/gmx.tc', ] ) )

### install

setup(name='SOLVSHIFT',
      version='2.0.0',
      description='Solvatochromic shifts from coarse-grained SolEFP theory',
      author='Bartosz Błasiak',
      author_email='blasiak.bartosz@gmail.com',
      url='no-page-yet',
      packages=['solvshift',
                #'solvshift.util',
                #'solvshift.diff',
               ],
      ext_modules=[SHFTEX,EFPROT,SHFTCE,EXREP,],
      data_files=frg_files+tc_files,
     )
