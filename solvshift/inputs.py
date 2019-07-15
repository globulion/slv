# ---------------------------------------- #
#          FF INPUT FILE MODULE            #
# ---------------------------------------- #

"""
autor wzorował się skryptem 
            ff.py 
Roberta W. Góry

https://gitlab.com/rgora/ff.py
"""

#from numpy  import numpy.array, float64, numpy.zeros
#from string import Template
#import re, os
#reflags = re.DOTALL
#from ff import *
#from math import ceil
#from units import *
#from utilities import Periodic
import re, os, ff, math, numpy, libbbg, string
from ff import *
reflags = re.DOTALL

__all__ = ['GAUSSIAN_INPUTS','GAMESS_INPUTS',]
__version__ = '3.3.5'

class INPUT_TEMPLATE(string.Template):
      """for easy handling input files from 
         template file stored on disk. 
         From ff.py of R. W. Góra"""
      delimiter = '@'

class INPUTS_FACTORY(FF,object):
      """general virtual class for inputs"""
      def __init__(self, xyz_file=0, frag_file=0, k=FF.k, L=0, redmass=0,
                   n_point=5, cartesian=False, mode_id=-1,
                   translation='',#all_cartesian=False,
                   sder_work_dir="/home/globula/programs/solvshift/sderiv/0.006/"):
          # step of differentiation and displacements
          INPUTS_FACTORY.setStep(k,n_point=n_point)
          INPUTS_FACTORY.setDispl(L=L,redmass=redmass,mode_id=mode_id)
          self.mode_id = mode_id
          # XYZ file 
          self.xyz_file = xyz_file
          # FRAGMENT FILE
          self.frag_file = frag_file
          # TRANSLATION 
          self.translation=translation
          # Mass-weighted eigenvector matrix L
          self.L = L
          # working directory for derivatives wrt normal modes
          self.sder_work_dir = sder_work_dir
          # --- 
          self.ReadTemplate()
          self.WriteStepFile()
          if cartesian:       self.WriteInputsCartesian()
          #elif all_cartesian: self.WriteInputsAllCartesian()
          else:               self.WriteInputs()
          
      def WriteStepFile(self):
          """documents differentiation step taken in the input files creation
          as well as step scaling factors"""
          
          if os.path.isfile('slv.step'):
             os.system('rm -f slv.step')
          
          step_file = open('slv.step','w')
          step_file.write('DIFFERENTIATION STEP (Å): %30.10f\n' % self.k)
          step_file.write('NUMERICAL METHOD POINTS: %36d\n' % self.n_point)
          step_file.write('PACKAGE: %57s\n' % self.package)
          step_file.write('SECOND DERIVATIVE MODE: %i\n' % self.mode_id)
          step_file.write('SECOND DERIVATIVE WORKING DIRECTORY: %s\n' % self.sder_work_dir)
          step_file.write('SECOND DERIVATIVE DIFFERENIATION STEP (DIMENSIONLESS): %.10f\n' % 10.0 )
          step_file.write('TRANSLATION: %s' % self.translation)
          step_file.close()
          ### do not allow to change the step file
          os.system('chmod a-w slv.step')
          
      def ReadFragmentsFile(self):
          """reads coordinates of fragment sites
          reading it from *.str file"""
          
      def ReadTemplate(self):
          """tries to read template from disk.
             if none, punches one on the disk. 
             From ff.py of R.W.Gora  """
          try: 
             self.templ = open( '%s.templ' % self.package ).read()
             self.templ = INPUT_TEMPLATE(self.templ)
          except IOError:
             print "\n              There is NO template file!\n I am creating it on the actual directory. Please check it!\n"
             print   "                 < %s.templ >\n\n" % self.package
             open('%s.templ' % self.package ,'w').write(self.templ)
             exit()
          except AttributeError:
             pass

      def WriteInputs(self):
          """writes the input files on the disk"""
          raise NotImplementedError 
          
      def WriteInputsCartesian(self):
          """writes the input files on the disk"""
          raise NotImplementedError 
 
      def Translate(self):
          """translates the input structure to the center 
          specified by keyword translation=...
          Note that translation is made in Angstroms!
          """
          
          translation_vector = numpy.zeros(3,dtype=numpy.float64)
          solute_structure = numpy.array(self.atoms[:,1:],dtype=numpy.float64)
          atoms = numpy.array(self.atoms[:,0])
          
          if self.translation.lower() == 'com':
             ### calculate center of mass in Angstrom
             mass_sum = 0
             for atom in range(len(atoms)):
                 #translation_vector+= libbbg.units.UNITS.mass[ libbbg.units.UNITS.atomic_numbers [atoms[atom]]  ] *\
                 translation_vector+= Atom([atoms[atom]]).mass *\
                                       solute_structure[atom]
                 #mass_sum+= libbbg.units.UNITS.mass[ libbbg.units.UNITS.atomic_numbers [atoms[atom]]  ]
                 mass_sum+= Atom([atoms[atom]]).mass
             translation_vector/=mass_sum
          elif self.translation.lower() == 'coe':
             ### calculate coe in Angstrom
             l_sum = 0
             vec = self.L[:,self.mode_id].reshape(len(atoms),3) * libbbg.units.UNITS.BohrToAngstrom
             for atom in range(len(atoms)):
                 translation_vector+= vec[atom]**2 * solute_structure[atom]
                 l_sum+= sum(vec[atom]**2)
             translation_vector/=l_sum
             
          return -translation_vector
 
 
class GAUSSIAN_INPUTS(INPUTS_FACTORY):
   """gaussian input file maker for FF calculations"""
   def __init__(self,xyz_file=0,frag_file=0,k=FF.k,L=0,redmass=0,
                n_point=5,cartesian=False,mode_id=-1,translation='',
                sder_work_dir=''):
       self.package = 'gaussian'
       self.templ = """\
%chk=@CHK
%mem=1800mb
%nproc=4
#p MP2/6-31G scf(conver=10,xqc) nosymm
# iop(7/33=1) GFInput freq(HPModes)
# integral(grid=199974) density=current
              
Finite Field computations

0 1
@DATA
""" 
       super(GAUSSIAN_INPUTS, self).__init__(xyz_file,frag_file,k,L,redmass,
                                             n_point,cartesian,mode_id,translation,
                                             sder_work_dir)

   def WriteInputs(self):
       """writes gaussian inputs on disk"""
       p=libbbg.utilities.Periodic(0)
       try:
          xyz=open(self.xyz_file,'r').readlines()
          nAtoms = int(xyz[0])
          xyz=xyz[2:int(xyz[0])+2]
       except ValueError:
          print "\n   The < %s > geometry file has incorrect syntax!\n" % self.xyz_file
          exit(1)
          
       # FRAGMENT FILE HANDLING
       try:
          frag=open(self.frag_file,'r').readlines()
          nFrags = int(frag[0])
          frag=frag[2:int(frag[0])+2]
       except ValueError:
          print "\n   The < %s > fragments' file has incorrect syntax!\n"\
                                % self.frag_file
          exit(1)
          
       # FRAGMENT FILE (COORDS IN AU!!!!!)
       data_frag=[]
       for i in range(len(frag)): data_frag.append(frag[i].split())
       inp_frag = numpy.array(data_frag[:])
       f_inp_frag = inp_frag[:,1:]
       f_inp_frag = numpy.array([ map( float, x ) for x in f_inp_frag ])
       
       ### translate the structure (if wanted)
       #f_inp_frag += self.Translate()
       
       FFFF  = f_inp_frag.copy()
       f_inp_frag*= libbbg.units.UNITS.AngstromToBohr
       #f_inp_frag+= (self.displacements[i] )
       s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
       inp_frag = inp_frag.tolist()
       inp_frag_dat = inp_frag[:]
       for j in range(nFrags):
           inp_frag[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
              % ( f_inp_frag[j][0], f_inp_frag[j][1], f_inp_frag[j][2], 
                   data_frag[j][0] )
       FRAG=''.join(inp_frag)
       FRAG= FRAG[:-1] # without blank line!
       # WRITE REFERENCE FRAGMENT FILE FOR SLV PACKAGE ROUTINES
       slv_frag_file = open('slv.frags','w')
       slv_frag_file.write(FRAG)
          
       # ------
       # write input files
       data=[]
       for i in range(len(xyz)): data.append(xyz[i].split())
       ### atom's list
       self.atoms = numpy.array(data)
       
       data_frag=[]
       for i in range(len(frag)): data_frag.append(frag[i].split())
       for i in range(len(self.displacements)):
           atom_id = math.ceil( i/4. )
           if i==0:
              displ_id = 0
           else:
              displ_id = i%4 + 4*int(not(i%4))
         #if i==self.mode_id:
           inp = numpy.array(data[:])
           f_inp = inp[:,1:]
           f_inp = numpy.array([ map( float, x ) for x in f_inp ])
           f_inp+=(self.displacements[i] * libbbg.units.UNITS.BohrToAngstrom)
           
           ### translate the structure (if wanted)
           f_inp += self.Translate()
           
           s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
           inp = inp.tolist()
           for j in range(nAtoms):
               inp[j] = "%-3s %16.10f %16.10f %16.10f\n"\
                  % (data[j][0],
                    f_inp[j][0], f_inp[j][1], f_inp[j][2] )
           XYZ=''.join(inp)
           
           inp_frag = numpy.array(data_frag[:])
           f_inp_frag = inp_frag[:,1:]
           f_inp_frag = numpy.array([ map( float, x ) for x in f_inp_frag ])

           #KKKK  = f_inp.copy()###################
           #KKKK*= libbbg.units.UNITS.AngstromToBohr################
           #for j in range(self.nAtoms):#################
           #    inp_string[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
           #       % ( KKKK[j][0], KKKK[j][1], KKKK[j][2], 
           #            data[j][0] )################
           #KKKKKKKK=''.join(inp_string)################
           #KKKKKKKK= KKKKKKKK[:-1] # without blank line!############
           # WRITE REFERENCE FRAGMENT FILE FOR SLV PACKAGE ROUTINES##############
           #slv_str_file = open('slv.str','w')#############
           #slv_str_file.write(KKKKKKKK)################           
           ### translate the structure (if wanted)
           #f_inp += self.Translate()

           #################################zdbgiSGfeosahgeasgaweg
           s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
           for j in range(nFrags):
               inp_frag_dat[j] = "%-3s %16.10f %16.10f %16.10f\n"\
                  % ( data_frag[j][0] , FFFF[j][0], FFFF[j][1], FFFF[j][2], )
           BULA=''.join(inp_frag_dat)
           BULA= BULA[:-1] # without blank line!  
                  
           f_inp_frag*= libbbg.units.UNITS.AngstromToBohr
           #f_inp_frag+= (self.displacements[i] )
           s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
           inp_frag = inp_frag.tolist()
           for j in range(nFrags):
               inp_frag[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
                  % ( f_inp_frag[j][0], f_inp_frag[j][1], f_inp_frag[j][2], 
                       data_frag[j][0] )
           FRAG=''.join(inp_frag)
           FRAG= FRAG[:-1] # without blank line!
           # WRITE REFERENCE FRAGMENT FILE FOR SLV
           if i==0: 
              slv_frag_file = open('slv.frags','w')
              slv_frag_file.write(FRAG)

           inpfile = self.xyz_file[:-4] +"_A%02.d" % (atom_id) + "_D%02.d" % (displ_id) + "_.inp"
           chkfile = self.xyz_file[:-4] +"_A%02.d" % (atom_id) + "_D%02.d" % (displ_id) + "_.chk"

           try:
              finput = self.templ.substitute(DATA=XYZ[:-1], CHK=chkfile, FRAG=BULA)
           except:
              finput = self.templ.substitute(DATA=XYZ[:-1], CHK=chkfile)
           open(inpfile,'w').write(finput+'\n')
             

       print "\n   The %d input files have been saved in the actual directory!\n"\
              % len(self.displacements)

   def WriteInputsCartesian(self):
       """writes gaussian cartesian inputs on disk"""
       p=libbbg.utilities.Periodic(0)
       try:
          xyz=open(self.xyz_file,'r').readlines()
          nAtoms = int(xyz[0])
          xyz=xyz[2:int(xyz[0])+2]
       except ValueError:
          print "\n   The < %s > geometry file has incorrect syntax!\n" % self.xyz_file
          exit(1)
          
       # FRAGMENT FILE HANDLING
       try:
          frag=open(self.frag_file,'r').readlines()
          nFrags = int(frag[0])
          frag=frag[2:int(frag[0])+2]
       except ValueError:
          print "\n   The < %s > fragments' file has incorrect syntax!\n"\
                                % self.frag_file
          exit(1)
          
       # FRAGMENT FILE (COORDS IN AU!!!!!)
       data_frag=[]
       for i in range(len(frag)): data_frag.append(frag[i].split())
       inp_frag = numpy.array(data_frag[:])
       f_inp_frag = inp_frag[:,1:]
       f_inp_frag = numpy.array([ map( float, x ) for x in f_inp_frag ])
       
       ### translate the structure (if wanted)
       #f_inp_frag += self.Translate()
       
       FFFF  = f_inp_frag.copy()
       f_inp_frag*= libbbg.units.UNITS.AngstromToBohr
       #f_inp_frag+= (self.displacements[i] )
       s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
       inp_frag = inp_frag.tolist()
       inp_frag_dat = inp_frag[:]
       for j in range(nFrags):
           inp_frag[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
              % ( f_inp_frag[j][0], f_inp_frag[j][1], f_inp_frag[j][2], 
                   data_frag[j][0] )
       FRAG=''.join(inp_frag)
       FRAG= FRAG[:-1] # without blank line!
       # WRITE REFERENCE FRAGMENT FILE FOR SLV PACKAGE ROUTINES
       slv_frag_file = open('slv.frags','w')
       slv_frag_file.write(FRAG)
       
       
       #################################zdbgiSGfeosahgeasgaweg
       s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
       for j in range(nFrags):
           inp_frag_dat[j] = "%-3s %16.10f %16.10f %16.10f\n"\
              % ( data_frag[j][0] , FFFF[j][0], FFFF[j][1], FFFF[j][2], )
       BULA=''.join(inp_frag_dat)
       BULA= BULA[:-1] # without blank line!  
       
       
       # ------
       # write input files
       data=[]
       for i in range(len(xyz)): data.append(xyz[i].split())

       # zero displacement file
       inp = numpy.array(data[:])
       inp_string = inp[:]
       inp_string = inp_string.tolist()       
       ### atom's list 
       self.atoms = inp
       #print inp
       
       f_inp = inp[:,1:]
       f_inp = numpy.array([ map( float, x ) for x in f_inp ])
       
       KKKK  = f_inp.copy()###################
       KKKK*= libbbg.units.UNITS.AngstromToBohr################
       for j in range(self.nAtoms):#################
           inp_string[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
              % ( KKKK[j][0], KKKK[j][1], KKKK[j][2], 
                   data[j][0] )################
       KKKKKKKK=''.join(inp_string)################
       KKKKKKKK= KKKKKKKK[:-1] # without blank line!############
       # WRITE REFERENCE FRAGMENT FILE FOR SLV PACKAGE ROUTINES##############
       slv_str_file = open('slv.str','w')#############
       slv_str_file.write(KKKKKKKK)################
       
       ### translate the structure (if wanted)
       #f_inp += self.Translate()
       
       s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
       inp = inp.tolist()
       for j in range(self.nAtoms):
           inp[j] = "%-3s %16.10f %16.10f %16.10f\n"\
              % (data[j][0],
                f_inp[j][0], f_inp[j][1], f_inp[j][2] )
       XYZ=''.join(inp)
       XYZ=XYZ[:-1]
       inpfile = self.xyz_file[:-4] + "_A%03.d"% (0) +"_D%02.d" % (0) + "_.inp"
       chkfile = self.xyz_file[:-4] + "_A%03.d"% (0) +"_D%02.d" % (0) + "_.chk"
       try:
           finput = self.templ.substitute(DATA=XYZ, CHK=chkfile)
       except:
           finput = self.templ.substitute(DATA=XYZ, CHK=chkfile, FRAG=BULA)
       open(inpfile,'w').write(finput+"\n")

       # DISPLACEMENTS
       for i in range(self.nAtoms):
         for d in range(len(self.DisplCart)):
           inp = numpy.array(data[:])
           f_inp = inp[:,1:]
           f_inp = numpy.array([ map( float, x ) for x in f_inp ])
           
           ### translate the structure (if wanted)
           #f_inp += self.Translate()
           
           for j in range(nAtoms):
               if i==j:
                  f_inp[j,:] += self.DisplCart[d,:]

           s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
           inp = inp.tolist()
           for j in range(self.nAtoms):
               inp[j] = "%-3s %16.10f %16.10f %16.10f\n"\
                  % (data[j][0],
                    f_inp[j][0], f_inp[j][1], f_inp[j][2] )
           XYZ=''.join(inp)
           XYZ=XYZ[:-1]
           


           inpfile = self.xyz_file[:-4] + "_A%03.d"% (i+1) +"_D%02.d" % (d+1) + "_.inp"  
           chkfile = self.xyz_file[:-4] + "_A%03.d"% (i+1) +"_D%02.d" % (d+1) + "_.chk"  

           try:
              finput = self.templ.substitute(DATA=XYZ, FRAG=BULA, CHK=chkfile)
           except error:
              finput = self.templ.substitute(DATA=XYZ, CHK=chkfile)
           open(inpfile,'w').write(finput+"\n")
           
             

       print "\n   The %d input files have been saved in the actual directory!\n"\
              % (self.nAtoms * 3 * 4 + 1)
              
 
   def WriteInputsAllCartesian(self):
       """writes gaussian cartesian inputs on disk"""

              
class GAMESS_INPUTS(INPUTS_FACTORY):
   """gaussian input file maker for FF calculations"""
   def __init__(self,xyz_file=0,frag_file=0,
                k=FF.k,L=0,redmass=0,
                n_point=5,
                cartesian=False,mode_id=-1,translation='',
                sder_work_dir=''):
                    
       self.package = 'gamess'
       self.templ = """\
 $system mwords=4 memddi=12 parall=.t. $end
 $contrl scftyp=rhf runtyp=energy icharg=0 mult=1 units=angs
         maxit=100 exetyp=check ispher=1
         icut=20 itol=30 aimpac=.f. mplevl=0 cctyp=none $end
 $scf    dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
         conv=1d-11  $end
 $basis  gbasis=STO ngauss=3 $end
 $data
<TYPE HERE YOUR JOB TITLE>
C1 0
@DATA
 $end
""" 
       super(GAMESS_INPUTS, self).__init__(xyz_file,frag_file,
                                           k,L,redmass,n_point,
                                           cartesian,mode_id,translation,
                                           sder_work_dir)

   def WriteInputs(self):
       """writes gaussian inputs on disk"""
       p=libbbg.utilities.Periodic(0)
       try:
          xyz=open(self.xyz_file,'r').readlines()
          nAtoms = int(xyz[0])
          xyz=xyz[2:int(xyz[0])+2]
       except ValueError:
          print "\n   The < %s > geometry file has incorrect syntax!\n" % self.xyz_file
          exit(1)
          
       # FRAGMENT FILE HANDLING
       try:
          frag=open(self.frag_file,'r').readlines()
          nFrags = int(frag[0])
          frag=frag[2:int(frag[0])+2]
       except ValueError:
          print "\n   The < %s > fragments' file has incorrect syntax!\n"\
                                % self.frag_file
          exit(1)
          
       # ------
       # write input files
       data=[]
       for i in range(len(xyz)): data.append(xyz[i].split())
       data_frag=[]
       for i in range(len(frag)): data_frag.append(frag[i].split())
       for i in range(len(self.displacements)):
           atom_id = math.ceil( i/4. )
           if i==0:
              displ_id = 0
           else:
              displ_id = i%4 + 4*int(not(i%4))
         #if i==self.mode_id:
           inp = numpy.array(data[:])
           f_inp = inp[:,1:]
           f_inp = numpy.array([ map( float, x ) for x in f_inp ])
           f_inp+= (self.displacements[i] * libbbg.units.UNITS.BohrToAngstrom)
           s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
           inp = inp.tolist()
           for j in range(nAtoms):
               inp[j] = "%-3s %7.1f %16.10f %16.10f %16.10f\n"\
                  % (data[j][0], libbbg.units.Atom(data[j][0]).atno,
                    f_inp[j][0], f_inp[j][1], f_inp[j][2] )
           XYZ=''.join(inp); XYZ=XYZ[:-1]
                             #% (data[j][0], libbbg.units.UNITS.atomic_numbers[data[j][0]],
           
           inp_frag = numpy.array(data_frag[:])
           f_inp_frag = inp_frag[:,1:]
           f_inp_frag = numpy.array([ map( float, x ) for x in f_inp_frag ])
           f_inp_frag*= libbbg.units.UNITS.AngstromToBohr
           #f_inp_frag+= (self.displacements[i] )
           s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
           inp_frag = inp_frag.tolist()
           for j in range(nFrags):
               inp_frag[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
                  % ( f_inp_frag[j][0], f_inp_frag[j][1], f_inp_frag[j][2], 
                       data_frag[j][0] )
           FRAG=''.join(inp_frag)
           FRAG= FRAG[:-1] # without blank line!
           # WRITE REFERENCE FRAGMENT FILE FOR SLV
           if i==0: 
              slv_frag_file = open('slv.frags','w')
              slv_frag_file.write(FRAG)

           inpfile = self.xyz_file[:-4] +"_A%02.d" % (atom_id) + "_D%02.d" % (displ_id) + "_.inp"    

           try:
              finput = self.templ.substitute(DATA=XYZ)
           except:
              finput = self.templ.substitute(DATA=XYZ)
           open(inpfile,'w').write(finput)
             

       print "\n   The %d input files have been saved in the actual directory!\n"\
              % len(self.displacements)

   def WriteInputsCartesian(self):
       """writes gaussian cartesian inputs on disk"""
       p=libbbg.utilities.Periodic(0)
       try:
          xyz=open(self.xyz_file,'r').readlines()
          nAtoms = int(xyz[0])
          xyz=xyz[2:int(xyz[0])+2]
       except ValueError:
          print "\n   The < %s > geometry file has incorrect syntax!\n" % self.xyz_file
          exit(1)
          
       # FRAGMENT FILE HANDLING
       try:
          frag=open(self.frag_file,'r').readlines()
          nFrags = int(frag[0])
          frag=frag[2:int(frag[0])+2]
       except ValueError:
          print "\n   The < %s > fragments' file has incorrect syntax!\n"\
                                % self.frag_file
          exit(1)

       # FRAGMENT FILE (COORDS IN AU!!!!!)
       data_frag=[]
       for i in range(len(frag)): data_frag.append(frag[i].split())       
       inp_frag = numpy.array(data_frag[:])
       f_inp_frag = inp_frag[:,1:]
       f_inp_frag = numpy.array([ map( float, x ) for x in f_inp_frag ])
       f_inp_frag*= libbbg.units.UNITS.AngstromToBohr
       #f_inp_frag+= (self.displacements[i] )
       s_inp_frag = numpy.array([ map( str  , x ) for x in f_inp_frag ])
       inp_frag = inp_frag.tolist()
       for j in range(nFrags):
           inp_frag[j] = "POINT %16.10f %16.10f %16.10f %9s\n"\
              % ( f_inp_frag[j][0], f_inp_frag[j][1], f_inp_frag[j][2],
                   data_frag[j][0] )
       FRAG=''.join(inp_frag)
       FRAG= FRAG[:-1] # without blank line!
       # WRITE REFERENCE FRAGMENT FILE FOR SLV
       slv_frag_file = open('slv.frags','w')
       slv_frag_file . write(FRAG)
         
       # ------
       # write input files
       data=[]
       for i in range(len(xyz)): data.append(xyz[i].split())
       # zero displacement file
       inp = numpy.array(data[:])
       f_inp = inp[:,1:]
       f_inp = numpy.array([ map( float, x ) for x in f_inp ])
       s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
       inp = inp.tolist()
       for j in range(self.nAtoms):
           inp[j] = "%-3s %7.1f %16.10f %16.10f %16.10f\n"\
              % (data[j][0], libbbg.units.UNITS.atomic_numbers[data[j][0]],
                f_inp[j][0], f_inp[j][1], f_inp[j][2] )
       XYZ=''.join(inp); XYZ=XYZ[:-1]
       inpfile = self.xyz_file[:-4] + "_A%03.d"% (0) +"_D%02.d" % (0) + "_.inp"  
       #try:
       finput = self.templ.substitute(FRAG=FRAG,DATA=XYZ)
       #except:
       #    finput = self.templ.substitute(DATA=XYZ)
       open(inpfile,'w').write(finput)      
        
       
       # DISPLACEMENTS       
       for i in range(self.nAtoms):
         for d in range(len(self.DisplCart)):
           inp = numpy.array(data[:])
           f_inp = inp[:,1:]
           f_inp = numpy.array([ map( float, x ) for x in f_inp ])
           for j in range(nAtoms):
               if i==j:
                  f_inp[j,:] += self.DisplCart[d,:]

           s_inp = numpy.array([ map( str  , x ) for x in f_inp ])
           inp = inp.tolist()
           for j in range(self.nAtoms):
               inp[j] = "%-3s %7.1f %16.10f %16.10f %16.10f\n"\
                  % (data[j][0], libbbg.units.UNITS.atomic_numbers[data[j][0]],
                     f_inp[j][0], f_inp[j][1], f_inp[j][2] )
           XYZ=''.join(inp)
           XYZ=XYZ[:-1]

           inpfile = self.xyz_file[:-4] + "_A%03.d"% (i+1) +"_D%02.d" % (d+1) + "_.inp" 
           finput = self.templ.substitute(FRAG=FRAG,DATA=XYZ)

           open(inpfile,'w').write(finput)
             

       print "\n   The %d input files have been saved in the actual directory!\n"\
              % (self.nAtoms * 3 * 4 + 1)
