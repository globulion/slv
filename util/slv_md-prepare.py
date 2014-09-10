#!/usr/bin/python
"""
 --------------------------------------------------------------------
                       SLV-MD CHARGE ENVIRONMENT
 --------------------------------------------------------------------

 This helper script is designed to prepare the environment input file
 for charge-based electrostatics calculations of frequency shifts usi
 ng SLV-MD routines. 
 
 Usage:
 
 ./slv_md-prepare.py [charges] [solute atoms] [atoms to zero-out]
 
 The output of the script is 'environment.inp' file which is to be us
 ed by SLV_MD class. The debug information for checking is stored in
 'debug.log' file.
 
 Notes:
 
  1) [charges] - file with charges. It can be: 
       a) gromacs .itp file
       b) gromacs .tpr file

  2) [solute atoms] - the string specifying numbers of atoms of a sol
     ute molecule. The numbers are in normal (not Python-like) conven
     tion.

  3) [atoms to zero-out] - the string specifying numbers of atoms who
     se charges are to be zeroed out. These atoms belong to EFP fragm
     ents. The numbers are in normal (not Python-like) convention. 
 
 Warning:
 
  - generally the 'environment.inp' file should be manually edited !
    The possible changes are neutralization of residual charges after
    cutting off EFP fragments, changing the charges for closer residu
    es and so on.

 --------------------------------------------------------------------
                                     Last revision: 10 Sep 2014
 --------------------------------------------------------------------
"""
print __doc__
from sys import argv, exit
if len(argv)==1: exit()
from utilities import text_to_list
from numpy import array, float64
import re, re_templates

# utilities
def parse_charges(c_file,
                  n_solvent_atoms=3,n_solvent_mol=0,
                  n_ion_atoms=1,n_ion_mol=0,c_ion=1.0,
                  all=False,debug=False):
    """parse charges from TPR formatted file or ITP file"""

    Q = []
    # GROMACS TPR 
    if c_file.lower().endswith('.tpr'):
       # open the formatted (gmx-dumped) TPR file                      
       a = open(c_file)
       text = a.read()
       a.close()
       
       # split to search for sections
       templ    = re.compile(r'      name=',re.DOTALL)
       sections = re.split(templ, text)[1:]
       
       print len(sections)
       
       # parse the charges
       templ = re.compile('q= ?' + re_templates.re_real_e ,re.DOTALL) 
       
       for section in sections:
           Q += re.findall(templ,section)
       for i in range(len(Q)):
           Q[i] = Q[i][2:]
       Q = array(Q,float64)

    # GROMACS ITP
    elif c_file.lower().endswith('.itp'):
       pass 

    # AMBER PRMTOP
    elif c_file.lower().endswith(".prmtop"):
       pass
 
    # return charges 
    if all:
       # debug info
       if debug:
          deb = open('debug.log','w')
          log = ' Debug info\n\n'
          log+= ' No of ALL charges: %25d\n\n' % len(Q)
          deb.write(log)
          deb.close()
       return Q
    else:
       chg_ep = Q[:-n_solvent_atoms-n_ion_atoms*n_ion_mol]
       chg_ion= Q[-n_solvent_atoms-n_ion_atoms*n_ion_mol:-n_solvent_atoms]
       chg_sol= Q[-n_solvent_atoms:] 
       chg_sol = chg_sol.tolist()
       chg_sol*= n_solvent_mol
       chg_sol = array(chg_sol,float64)
       # debug info
       if debug:
          deb = open('debug.log','w')
          log = ' Debug info\n\n'
          log+= ' No of EPROTEIN charges: %25d\n' % len(chg_ep)
          log+= ' No of ION      charges: %25d\n' % len(chg_ion)
          log+= ' No of SOLVENT  charges: %25d\n' % len(chg_sol)
          log+= ' \n'
          log+= ' Ion chrges : ' + str(chg_ion) + '\n'
          log+= ' Sol charges: ' + str(chg_sol) + '\n\n'
          deb.write(log)
          deb.close()

       return chg_ep, chg_ion, chg_sol



def write_input(chg_ep,solat,zout,
                n_ion_mol=0,c_ion=1.0,name_ion='Sodium cation'):
    """write the 'environment.inp' file"""
    # zero-out the charges
    chg_ep[solat]= 0.00
    for z in zout:
        chg_ep[z] = 0.00
    #
    out = open('environment.inp','w')
    log = '[ eprotein ]\n\n'
    log+= 'solute atoms='
    at  = ''
    for i in solat:
        at+= '%s' % i + ','
    log+= at[:-1]+'\n\n'
    log+= 'charges\n'
    for i in range(len(chg_ep)):
        log+= '%5d %10.6f\n' % (i,chg_ep[i])
    log+= 'end charges\n'
    log+= '\n'
    #
    log+= '[ ions-1 ]\n\n'
    log+= 'name=%s\n' % name_ion 
    log+= 'nmol=%d\n' % n_ion_mol
    log+= 'charges=%.4f\n' % c_ion
    log+= '\n'
    #
    log+= '[ solvent-1 ]\n\n'
    log+= '[ end ]\n'
    #
    out.write(log)
    out.close()
    return

# create the indices of atoms (in Python convention)
solat= text_to_list(argv[2],delimiter=',') - 1
zout = argv[3].split(':')

for i in range(len(zout)):
    zout[i] = text_to_list(zout[i],delimiter=',') - 1

#print solat
#print zout

# parse the charges
chg_ep, chg_ion, chg_sol = parse_charges(argv[1],
                                        n_solvent_atoms=3,n_solvent_mol=3,
                                        n_ion_atoms=1,n_ion_mol=4,c_ion=1.0,
                                        all=False,debug=True)

print len(chg_ep)
print chg_ion
print chg_sol
# create the input file
write_input(chg_ep,solat,zout,
            n_ion_mol=4,c_ion=1.0,name_ion='Sodium cation')
