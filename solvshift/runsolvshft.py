#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                ---------------------------------------------------------                    *
 *                             --- WELCOME TO SOLVSHIFT ਉ  ---                                 *
 *                                                                                             *
 *                The script has to perform computations of frequency shift                    *
 *                of a solute imerged in a solvent in terms of discrete and                    *
 *                implicit solvatochromic models and its extended versions.                    *
 *                Available models are:                                                        *
 *                  1) Cho-Onsager continuum model                                             *
 *                  2) Discrete SolX models and its arbitrary contractions.                    *
 *                     The SolX models available are SolDMA and SolMMM                         *
 *                                                                                             *
 *                ---------------------------------------------------------                    *
 *                Usage:                                                                       *
 *                      python solvshift.py -h [--help]                                        *
 *                                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 """
def Usage():
    print \
 """
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                           SOLVSHIFT OPTIONS                                               *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [1] OPERATION MODE                                                                                       *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -h, --help               : display help                                                               *
 *     -v, --version            : display package info                                                       *
 *     -i, --inputs    [file]   : input mode - create finite field input files                               *
 *                              : [file]  file containing structure (in Å! and XYZ format)                   *
 *     -c, --calc               : calculation mode                                                           *
 *     -g, --gaussian           : sets Gaussian log files to be parsed (relevant if -c option is used)       *
 *                              : as well as gaussian input files were created (relevant if -i option used)  *
 *                              :             / default is Gamess                                            *
 *     -H, --ncpus     [int]    : number of processors (default is one)                                      *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [2] DIFFERENTIATION                                                                                      *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -s, --step      [float]  : set step for the differentiation wrt cartesian coordinates                 *
 *                              :             / the step is in Å!                                            *
 *     -n, --n-point   [int]    : specifies type of Finite Field procedure (e.g. -n 7)                       *
 *                              :             / defult poinity is 5                                          *
 *     -x, --xyz       [1/0]    : differentiation in cartesian space                                         *
                                              / default --xyz=1                                              *
 *     -W, --sder-work [dir]    : sderiv working directory to be written in slv.step file                    *
 *                              :             / default './sder'                                             *
 *     -P, --print-der          : print the derivatives of (distributed) multipole moments                   *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [3] PARAMETERS SPECIFICATION                                                                             *
 *  -------------------------------------------------------------------------------------------------------  *
 *     The default is SolCharge model                                                                        *
 *     -O, --solmmm             : calculate molecular solvatochromic multipole moments (MSMA)                *
 *     -d, --solcamm            : calculate distributed solvatochromic multipole moments (DSMA)              *
 *     -y, --mixed              : create mixed solute-solvent model (solvent DSMA is constructed             *
 *                              : from atomic charges)                                                       *
 *     -u, --make-ua   [UAs]    : switch on united atom mode                                                 *
 *                              :             / examples: 2,3,5,6 or 4,5,6,9-12,14,15,18-2,1                 *
 *                              :             / the first number is UA, the next before '-' are atoms to be  *
 *                              :             / contracted. If no '-' it means only one UA group defined     *
 *     -T, --transl    [x]      : specifies translation of the input structure to the given center x         *
 *                              :             / x can be 'com' - center of mass or                           *
 *                              :             /          'coe' - center of square mass-weighted eigenvector  *
 *                              :             /          'at,X,Y,u' - weighted center between atoms X and Y  *
 *                              :             /                       with weighting ratio 'u'               *
 *     -G, --allign             : allign the multipoles along xyz axes                                       *
 *     -k, --check              : sum up the DSMA to MSMA and print                                          *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [4] SOLUTE/SOLVENT SPECIFIC INFORMATION                                                                  *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -a, --anh       [file]   : use anharmonic file for solute                                             *
 *                              : [file] containing harmonic and anharmonic data                             *
 *     -m, --mode-id   [int]    : provide normal mode ID number in Python convention (i.e. N-1)              *
 *     -F, --frag      [file]   : provide xyz file containing non-atomic fragments to be parametrized        *
 *     -R, --read      [file]   : read the parameters from file (coulomb format)                             *
 *     -b, --bsm       [file]   : parameters for benchmark solvent molecule (BSM) in coulomb format          *
 *     -z  --camm      [file]   : read electrostatic moments of solute molecule (for shoft corrections)      *
 *     -Z  --fderiv-j  [file]   : read derivatives of electrostatic moments of solute molecule (for shift    *
 *                              : corrections)                                                               *
 *     -A, --solute    [atoms]  : provide a string of atomic symbols for solute separated by a comma and     *
 *                              : without space, eg.: -A C,C,O,H,H,H,H                                       *
 *     -B, --solvent   [atoms]  : provide a string of atomic symbols for solvent separated by a comma and    *
 *                              : without space, eg.: -A C,C,O,H,H,H,H                                       *
 *     -L, --suplist   [list]   : provide superimposition indices. '5' means superimpose first 5 atoms of    *
 *                              : solute; '1,2,4,6' means superimpose first, second, fourth and sixth atom   *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [5] FREQUENCY SHIFTS                                                                                     *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -f, --freq               : calculate frequency shift                                                  *
 *     -X, --corr               : evaluate corrections to frequency shifts                                   *
 *     -t, --typ       [type]   : provide the file of solute/solvent system names (relevant if -f option     *
 *                              : is used)                                                                   *
 *     -D, --target    [dir]    : provide the directory with solute/solvent system files                     *
 *     -M, --MD        [traj]   : calculate frequency shift distribution from Molecular Dynamics simulation  *
 *                              : providing the trajectory [traj] / default are 'amber' files /              *
 *                              : charges file as the args[0]                                                *
 *     -U, --md-package         : specify MD package software (GROMACS or AMBER)                             *
 *     -Q, --struct             : evaluate solute structural distortions in the first order due              *
 *                              : to solvation                                                               *
 *     -E, --eds                : evaluate frequency shifts from direct differentiation of interaction       *
 *                              : energies (default is EDS - Hybrid Variational-Perturbational Interaction)  *
 *                              : Energy Decomposition Scheme)                                               *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [6] CHO-ONSAGER MODEL                                                                                    *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -o, --onsager            : calculate Cho-Onsager C_MA and C_EA coefficients                           *
 *     -C, --cavity    [float]  : provide cavity radius in Å for Cho-Onsager calculations                    *
 *     -e, --epsilon   [float]  : provide dielectric donstant of continuum. If not specified                 *
 *                              : the calculations are performed for 30 epsilon values from                  *
 *                              : 1.1 to 100                                                                 *
 *     -p, --pol                : include molecular polarizability                                           *
 *     -j, --max-iter  [int]    : apply flexible model (if specified) using maximum number of                *
 *                              : iterations                                                                 *
 *     -l, --threshold [float]  : set convergence threshold for solvation energy (relevant if                *
 *                              : -l option specified)                                                       *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [7] INPUT/OUTPUT FILE HANDLING                                                                           *
 *  -------------------------------------------------------------------------------------------------------  *
 *     -R, --read      [file]   : read the parameters from file (coulomb format)                             *
 *     -S, --save               : save the parameters                                                        *
 *     -N, --name      [name]   : provide the name of output file                                            *
 *     -w, --cube      [file]   : write a solvatochromic potential/ energy density in a file                 *
 *     -V, --make-sol  [X] [Y]  : create solute_typ and solvent_typ files. Usage:                            *
 *                              :   - X has a format N,typ (e.g. 12,3A+B+C). N is the number of solute atoms *
 *                              :   - Y file with CHELPG analysis for target structure                       *
 *  -------------------------------------------------------------------------------------------------------  *
 *  [!] Machine epsilon is %E for float64 type                                                     *
 *  -------------------------------------------------------------------------------------------------------  *
 *  Globulion@HaveFun!                   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                       *
 * * * * * * * * * * * * * * * * * * * * *
""" % finfo(float64).eps

from head import *

#     Copyright © 2012, Bartosz Błasiak (globula@o2.pl)
#  
#     This program is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#  
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#  
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

# ----------------------------------------------------------------------------------------------------
__authors__         = " Bartosz Błasiak (globula@o2.pl) "
__version__         =   filter(str.isdigit, "$Revision: 26 $")
# ----------------------------------------------------------------------------------------------------

# -------------------- #
#    SOLVSHIFT INFO    #
# -------------------- #
def Intro():
    """Print usage information and exit. Zgapione z ff.py Roberta W. Góry."""
    os.system('clear')
    print __doc__
    print " Machine epsilon is: ",finfo(float64).eps,"for float64 type\n"
    print;exit()

def Version():
    """Print version and other information"""
    from __init__ import __version__
    print __doc__
    print 
    print " SOLVSHIFT.py ਉ  © , version ", __version__
    print __authors__
    print

# ----------------- #
#    MAIN ROUTINE   #
# ----------------- #
def Main(argv):
    """ solvshift.py main routine """

    ## --------------------------- ##
    ##   T H E   D E F A U L T S   ##
    ## --------------------------- ##
    
    ### operational options                                            # comments and annotations
    ncpus          = 1                                                 # number of cpu units
    Inputs         = False                                             # input mode
    Calculate      = False                                             # calculation mode
    Gamess         = True                                              # Gamess input files
    Gaussian       = False                                             # Gaussian input files
    Fragments      = False                                             # non-atomic cites
    step           = 0.05                                              # step for cartesian differentiation
    n_point        = 5                                                 # pointity of differentiation
    mode_id        =-1                                                 # normal mode index
    file_type      = 'gamess'                                          # 
    cartesian      = True                                              # make cartesian displacements
    sder_work_dir  = './sder'                                          # directory for sder
    make_ua        = False                                             # contracted models
    ua_list_5      = [ ( 1,12,11, 6), ( 4, 7, 9,10), ( 2, 8) ]         # UA; methyl groups + NH amide
    ua_list        = [ ( 4, 5, 6, 7) ]                                 # UA
    ua_list        = [ ( 1,12,11, 6), ( 4, 7, 9,10) ]                  # UA; methyl groups
    suplist        = [ 0, 1, 2, 3, 4]                                  # superimposition indices
    
    ### input/output file handling
    read_file      = ''                                                # read parameters
    make_cube      = False                                             # create cube
    make_coulomb   = False                                             # 
    out_name       = 'slv.par'                                         # output file name
    make_sol       = '12,default'                                      # solute/solvent target files

    ### parameters
    SolCHELPG      = True                                              # SolCHELPG mode
    SolMMM         = False                                             # SolMMM mode
    SolCAMM        = False                                             # SolCBAMM mode
    SolPOL         = False                                             # SolPOL mode
    mixed          = False                                             # mixed mode
    translation    = ''                                                # center of MMM
    allign         = False                                             # allign DMA in some ways
    charges        = {'H': 0.34242, 'O': -0.777834, 'Na': 1.00}        # BSM charges
    zero_camm      = ''                                                # for correction terms: moments
    fderiv_j       = ''                                                # for correction terms: fderivs of moments
    
    ### checking and troubleshooting
    print_derivatives = False                                          # prints derivatives of DMTP
    check             = False                                          # summ DMTP to MMM
    
    ### frequency shift (SLV)
    frequency_shift= False                                             # calculate shifts
    corrections    = False                                             # calculate correction terms
    eds            = False                                             # calculate shifts from EDS differentiation
    types          = 'types.txt'                                       # file with stems of target system files
    bsm_file       = ''                                                # BSM file
    solute_atno    = ['C','N','C','C','O','H','H','H','H','H','H','H'] # NMA
    solvent_atno   = ['O','H','H']                                     # H2O
    structural_change = False                                          # evaluate dQ
    target_dir     = os.environ.get('TARGET_DIR')                      # directory with target system files
    
    ### Cho-Onsager
    onsager        = False                                             # calculate shifts from continuum models
    ons_pol        = False                                             # including polarization explicitly
    cavity_radius  = 7.000                                             # in Bohr
    epsilon        = 0                                                 # dielectric constant
    iterate        = False                                             # apply iterative algorithm
    max_iter       = 0                                                 # maximum no of iterations
    threshold      = 1.E-6                                             # threshold for polarization energy convergence

    ### molecular dynamics frequency shifts (SLV_MD)
    md             = False                                             # calculate shifts from MD trajectory
    md_package     = 'amber'                                           # MD package
    md_charges     = ''                                                # file with charges
    md_trajectory  = ''                                                # file with trajectory
    md_input       = ''                                                # input file with charge information
    

    ## --------------------------- ##
    ##        O P T I O N S        ##
    ## --------------------------- ##
    
    try:
        opts, args = getopt.getopt(argv, "hi:a:s:cF:n:gx:fM:ot:m:T:OC:ydu:b:pW:A:B:kQj:l:e:D:w:PN:SR:GEV:z:Z:XU:H:L:" ,
                                        ["help"      ,
                                         "inputs="   ,
                                         "anh="      ,
                                         "step="     ,
                                         "calc="     ,
                                         "frag="     ,
                                         "n-point="  ,
                                         "gaussian"  ,
                                         "xyz="      ,
                                         "freq-shift",
                                         "MD="       ,
                                         "onsager"   ,
                                         "typ="      ,
                                         "mode-id="  ,
                                         "transl="   ,
                                         "solmmm"    ,
                                         "cavity="   ,
                                         "mixed"     ,
                                         "solcamm"   ,
                                         "make-ua="   ,
                                         "bsm="      ,
                                         "pol"       ,
                                         "sder-work=",
                                         "solute="   ,
                                         "solvent="  ,
                                         "check"     ,
                                         "struct"    ,
                                         "target-dir=",
                                         "threshold=",
                                         "epsilon="  ,
                                         "max-iter=" ,
                                         "cube="     ,
                                         "print-derivatives",
                                         "name="     ,
                                         "save"      ,
                                         "read="     ,
                                         "allign"    ,
                                         "eds"       ,
                                         "make-sol=" ,
                                         "camm="     ,
                                         "fderiv-j=" ,
                                         "correction-terms",
                                         "md-package=",
                                         "ncpus="     ,
                                         "suplist="]    )
    except getopt.GetoptError, error:  
        print "\n   SLV: Error in command line! Quitting...\n"
        exit()
    if not argv: 
       Version()
    for opt, arg in opts:
        if opt in ("-h", "--help" ):
           Usage()
        if opt in ("-v", "--version" ):
           Version()
        if opt in ("-i", "--inputs" ):
           Inputs   = True 
           xyz_file = arg
        if opt in ("-W", "--sder-work" ):
           sder_work_dir = arg
        if opt in ("-H", "--ncpus"):
           ncpus = int(arg)
        if opt in ("-T", "--transl" ):
           translation   = arg 
        if opt in ("-w", "--cube" ):
           make_cube   = arg 
        if opt in ("-g", "--gaussian"):
           Gaussian = True ; file_type="gaussian"
           Gamess   = False
        if opt in ("-F", "--frag" ):
           Fragments   = True 
           frag_file = arg
        if opt in ("-a", "--anh"):
           anh = arg
        if opt in ("-l", "--threshold"):
           threshold = float64(arg)
        if opt in ("-e", "--epsilon"):
           epsilon = float64(arg)
        if opt in ("-j", "--max-iter"):
           max_iter = int(arg)
           iterate = True
        if opt in ("-s", "--step"):
           step = float64(arg)         
        if opt in ("-n", "--n-point"):
           n_point = int(arg)        
        if opt in ("-c", "--calc"):
           Calculate = True
        if opt in ("-O", "--solmmm"):
           SolMMM = True
           SolCHELPG = False
        if opt in ("-x", "--xyz"):
           if arg=='0': cartesian=False
        if opt in ("-o", "--onsager"):
           onsager = True     
        if opt in ("-f", "--freq-shift"):
           frequency_shift = True
        if opt in ("-d", "--solcamm"):
           SolCAMM   = True
           SolCHELPG = False
        if opt in ("-b", "--bsm"):
           bsm_file = arg
        if opt in ("-D", "--target-dir"):
           target_dir = arg
        if opt in ("-k", "--check"):
           check = True
        if opt in ("-A", "--solute"):
           solute_atno = array(arg.split(','),dtype=str)
        if opt in ("-B", "--solvent"):
           solvent_atno = array(arg.split(','),dtype=str)
        if opt in ("-C", "--cavity"):
           cavity_radius = float64(arg) *UNITS.AngstromToBohr
        if opt in ("-p", "--pol"):
           ons_pol = True
        if opt in ("-M", "--MD"):
           md = True
           md_trajectory = arg
           md_charges = args[0]
           md_input = args[1]
        if opt in ("-U", "--md-package" ):
           md_package   = arg.lower()
        if opt in ("-t", "--typ"):
           types = arg
        if opt in ("-m", "--mode-id"):
           mode_id = int(arg)
        if opt in ("-y", "--mixed"):
           mixed = True
        if opt in ("-u", "--make-ua"):
           make_ua = True
           ua_list = ParseUnitedAtoms(arg)
        if opt in ("-Q", "--struct"):
           structural_change = False
        if opt in ("-P", "--print-derivatives" ):
           print_derivatives   = True
        if opt in ("-N", "--name" ):
           out_name   = arg
        if opt in ("-S", "--save" ):
           make_coulomb   = True
        if opt in ("-R", "--read" ):
           read_file   = arg
        if opt in ("-E", "--eds" ):
           eds   = True
        if opt in ("-G", "--allign" ):
           allign   = True
        if opt in ("-z", "--camm" ):
           zero_camm   = arg
        if opt in ("-Z", "--fderiv-j" ):
           fderiv_j   = arg
        if opt in ("-X", "--correction-terms" ):
           corrections   = True
        if opt in ("-V", "--make-sol" ):
           N,typ = arg.split(',')
           N = int(N)
           solute_ids   = asarray(linspace(0,N-1,N),dtype=int)
           MakeSoluteAndSolventFiles(args[-1], typ, solute_ids, charges )
        if opt in ("-L", "--suplist" ):
           if ',' in arg: 
              m1 = lambda i: i-1
              suplist = map(int,arg.split(','))
              suplist = map(m1 ,suplist)
           else: suplist = [x for x in range(int(arg))]
                    
    ### --- choose the task --------------------------------------------------
    
    # ---------------------------------------------------------------------- #
    #                              INPUT MODE                                #
    # ---------------------------------------------------------------------- #
    
    if Inputs:
       ### [1.1] get gas-phase molecular specific data
       if not Fragments: frag_file = xyz_file
       ANH = FREQ(anh)
       ### [1.2] GAMESS inputs
       if Gamess:
          GAMESS_INPUTS(xyz_file=xyz_file,
                        frag_file=frag_file,
                        k=step,
                        L=ANH.L,
                        redmass=ANH.redmass,
                        n_point=n_point,
                        cartesian=cartesian,
                        mode_id=mode_id,
                        sder_work_dir=sder_work_dir)
                        
       ### [1.3] GAUSSIAN inputs
       elif Gaussian:
          GAUSSIAN_INPUTS(xyz_file=xyz_file,
                          frag_file=frag_file,
                          k=step,
                          L=ANH.L,
                          redmass=ANH.redmass,
                          n_point=n_point,
                          cartesian=cartesian,
                          mode_id=mode_id,
                          translation=translation,
                          sder_work_dir=sder_work_dir)

    # ---------------------------------------------------------------------- #
    #                           CALCULATION MODE                             #
    # ---------------------------------------------------------------------- #
    
    if Calculate:
       
       ### [2.1] get gas-phase molecular specific data
       ANH = FREQ(anh)
       redmass=ANH.redmass
       freq=ANH.freq
       if print_derivatives:
          ANH.FDeriv(Print=1,Debye=0,divide=0)
          ANH.Intens_1(Print=1,Debye=0)
          
       ### [2.2] SolX - distributed solvatochromic multipole model mode
       if not onsager:
          
          ### [2.2.1] read external files in COULOMB format
          if read_file:
             parameters = ParseDMA(read_file,'coulomb')
             if fderiv_j: 
                fderiv  = ParseDMA(fderiv_j ,'coulomb')
             
             ### allign the parameters
             if allign:
                from utilities import Allign
                #alligned_dma = parameters.OverallMoments_old()
                #coe=ANH.COE(structure=parameters.pos)[mode_id]
                #new_origin = array([[ 0.13156235,  0.53434119,  1.8956809 ]])
                #alligned_dma.MAKE_FULL()
                #alligned_dma.ChangeOrigin(new_origin_set=new_origin)
                #Allgn=Allign(xyz=parameters.copy().pos,atid=(3,5,2),dma=alligned_dma.copy())
                #alligned_dma, xyz =Allgn.get_transformed()
                #alligned_dma.MAKE_FULL()
                #alligned_dma.MakeTraceless()
                #alligned_dma.makeDMAfromFULL()
                #alligned_dma.ChangeUnits(1,UNITS.AuToStarkTunningRate,
                #                           UNITS.AuToQuadrupolarTunningRate*1e8*2,
                #                           UNITS.AuToOctupolarTunningRate*1e16*2)
                print parameters
                #print alligned_dma

          ### [2.2.2] calculate parameters!
          else:
              ### calculate the derivatives of distributed multipole moments!
              CALC=DIFF(freq=freq,
                        dir='.',
                        cartesian=cartesian,
                        L=ANH.L,
                        camm=SolCAMM,
                        eds=eds)

              fdip=0;sdip=0;dipole=0

              ### print the derivatives of distributed multipole moments
              if print_derivatives: print CALC
          
              ### [!] --- SolCHELPG, SolCAMM or SolDMA
              if not SolMMM:
                 PARAM = MCHO(fderiv=CALC.Fder,
                              sderiv=CALC.Sder,
                              fpol=CALC.fpol,
                              spol=CALC.spol,
                              pol=ons_pol,
                              gijj=ANH.K3,
                              freq=freq,
                              redmass=redmass,
                              onsager=False,
                              overall_MM=SolMMM,
                              MM_mode_id=mode_id,
                              eds=eds,
                              fEDS=CALC.fEDS,
                              sEDS=CALC.sEDS)
                 if eds: print PARAM.eds_shifts 

                 if ons_pol: 
                    for i in PARAM.solpol: print i
                    SolPOL = PARAM.solpol[-1]
                    
                 ### make contracted united-atom-based SolX-n models
                 if make_ua:
                    if SolCAMM:
                       PARAM[mode_id].MakeUa(ua_list,change_origin=True)
                    else:
                       PARAM[mode_id].MakeUa(ua_list,change_origin=False)
                 
                 ### update the properties
                 PARAM.Update(PARAM)
                 
                 ### visualize the parameters
                 if make_cube:
                    cube = PotentialMap(xyzfile=make_cube,dma=PARAM[mode_id],nxyz=(120,120,120))
                    cube.eval()
                    cube.tofile('slv.cube')
                    
                 ### print the all parameters and Stark tunning rates
                 #print PARAM
                 #print read_file
                 ### check
                 if check:
                    atoms = PARAM.atoms
                    for i in atoms: print i
                    print PARAM[mode_id].OverallMoments()
                
              ### [!] --- SolMMM
              else:
                 dipolederivatives = DipDeriv(file=anh,L=ANH.L,step=step)
                 fdip = dipolederivatives.fdip
                 xyz,sdip = dipolederivatives.sdip
                 dipole=dipolederivatives.dipole
             
                 PARAM = MCHO(fderiv=CALC.Fder,
                              sderiv=CALC.Sder,
                              gijj=ANH.K3,
                              fdip=fdip,
                              sdip=sdip,
                              sdip_full=xyz,
                              fder_MM=CALC.FDip,
                              sder_MM=CALC.sder_MM,
                              MM_mode_id=mode_id,
                              freq=freq,
                              redmass=redmass,
                              overall_MM=SolMMM)

 
                 MA,EA,full=PARAM.parameters_MM
                 SolMMM = MA+EA
                 EA_MM = EA
                 MA_MM = MA
                 log = "\n"
                 log+= " "+"="*33+"\n"
                 log+= " SOLVATOCHROMIC MOMENTS [A.U]\n\n" 
                 log+= " "+"="*33+"\n\n\n"
                 log+= str(MA)
                 log+= str(EA)
                 log+= str(full)
                 print log
                 
                 ### allign the parameters
                 if allign:
                    print " OVERALL SOLVATOCHROMIC MOMENTS IN UNITS OF GLORY\n"
                    alligned_dma = full.copy()
                    fff = PARAM[mode_id].copy()
                    #alligned_dma = fff.OverallMoments_old()
                    coe=ANH.COE(structure=PARAM.fragments)[mode_id]
                    new_origin = array([coe])
                    alligned_dma.MAKE_FULL()
                    alligned_dma.ChangeOrigin(new_origin_set=new_origin)
                 
                    from utilities import Allign
                    Allgn=Allign(xyz=PARAM.fragments,atid=(3,5,2),dma=alligned_dma)
                    alligned_dma, xyz =Allgn.get_transformed()
                    print xyz #* UNITS.BohrToAngstrom
                    alligned_dma.MAKE_FULL()
                    alligned_dma.MakeTraceless()
                    alligned_dma.makeDMAfromFULL()
                    alligned_dma.ChangeUnits(1,UNITS.AuToStarkTunningRate,
                                               UNITS.AuToQuadrupolarTunningRate*1e8*2,
                                               UNITS.AuToOctupolarTunningRate*1e16*2)
                    print alligned_dma


              ### store the parameters in memory
              if SolMMM:
                 parameters = PARAM.parameters_MM[-1]  # list [EA,MA,full]
              else:
                 parameters = PARAM[mode_id]
               
              ### write parameters
              if make_coulomb:
                 if SolMMM:
                    parameters.write(out_name) # writes total anharmonic parameters only (-1)
                 else:
                    parameters.write(out_name)


       ### [2.3] Cho-Onsager solvatochromic model mode
       else:
          ### calculate the dipole moment derivatives!
          dipolederivatives = DipDeriv(file=anh,L=ANH.L,step=step)
          fdip = dipolederivatives.fdip
          xyz,sdip = dipolederivatives.sdip
          dipole=dipolederivatives.dipole
          ANH.Intens_1(Debye=0,Print=1)
          dipolederivatives.SDeriv(Print=1)

          ### calculate the derivatives of distributed multipole moments!
          CALC=DIFF(freq=freq,
                    dir='.',
                    cartesian=cartesian,
                    L=ANH.L,
                    camm=SolCAMM)

          PARAM = MCHO(fderiv=fdip,sderiv=0,
                       gijj=ANH.K3,freq=freq,redmass=redmass,
                       onsager=onsager,
                       cavity_radius=cavity_radius,
                       fdip=fdip,
                       sdip=sdip,
                       sdip_full=xyz,
                       dipole=dipole,
                       overall_MM=SolMMM,
                       fder_MM=CALC.FDip,
                       sder_MM=CALC.sder_MM,
                       MM_mode_id=mode_id,
                       ons_pol=ons_pol,
                       step_cart=CALC.step,
                       step_mode=CALC.sderiv_step,
                       pol_1=CALC.polarizability_set_1,
                       pol_2=CALC.polarizability_set_2,
                       L=ANH.L,
                       iterate=iterate,
                       max_iter=max_iter,
                       threshold=threshold,
                       epsilon=epsilon)#,
                       #giijk=ANH.K4)


          log = "\n CHO-ONSAGER COEFFICIENTS [A.U.] \n FOR MODE %15.2f [cm-1]\n"%ANH.freq[mode_id]
          log+= "\n C_MA %10.6f"% (PARAM.onsager[mode_id][0])
          log+= "\n C_EA %10.6f"% (PARAM.onsager[mode_id][1])
          log+= "\n C_QF %10.6f"% (PARAM.onsager[mode_id][2])
          log+="\n"
          print log
       

       ## ------------------------------------------- ##
       ##        F R E Q U E N C Y   S H I F T        ##
       ## ------------------------------------------- ##
       
       if frequency_shift:
          ## ----------------------------------------------------- ##
          ##        SHIFTS FROM:                                   ##
          ##          * DISTRIBUTED SOLVATOCHROMIC CHARGES         ##
          ##          * MOLECULAR SOLVATOCHROMIC MULTIPOLES        ##
          ##          * DISTRIBUTED SOLVATOCHROMIC MULTIPOLES      ##
          ## ----------------------------------------------------- ##
          if Gaussian:
            out = open('shifts.dat','w') 
            print  >> out, " SLV FREQUENCY SHIFT REPORT. ALL VALUES IN [cm-1]"
            print  >> out, " ------------------------------------------------"
            if corrections:
               print  >> out, " %9s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s"% ('SYSTEM','R-1','R-2','R-3','R-4','R-5',
                                                                                                  'R-1','R-2','R-3','R-4','R-5')
            else:
               print  >> out, " %9s %13s %13s %13s %13s %13s"% ('SYSTEM','R-1','R-2','R-3','R-4','R-5')
            
            for typ in open(types).read().split('\n')[:-1]:
             solute  = ParseDMA( target_dir+'/solute_%s'%typ,file_type)
             solvent = ParseDMA( target_dir+'/solvent_%s'%typ,file_type)
             try:
                if not fderiv_j: fderiv = CALC.Fder[mode_id]
             except UnboundLocalError:
                fderiv = []
                
             if SolCHELPG:
                   
                f = SLV(pkg="gaussian",
                             nsatoms=solvent_atno,
                             nstatoms=solute_atno,
                             L=ANH.L,
                             mode_id=mode_id,
                             solute=parameters,
                             solute_structure=solute.get_pos(),
                             solvent=solvent,
                             solute_origin=translation,
                             mixed=mixed,
                             overall_MM=SolMMM,
                             camm=False,
                             chelpg=parameters,
                             structural_change=structural_change,
                             fderiv=fderiv,
                             redmass=ANH.redmass,
                             freq=ANH.freq,
                             ref_structure=parameters.get_pos(),
                             suplist=suplist)

                print  >> out, "%10s"%typ,
                print  >> out,"%13.2f "*5%tuple( f.shift[0] )
                print Emtp.log
                
             if SolCAMM:
                
                f = SLV(pkg="gaussian",
                             nsatoms=solvent_atno,
                             nstatoms=solute_atno,
                             L=ANH.L,
                             mode_id=mode_id,
                             solute=parameters,
                             solute_structure=solute.get_pos(),
                             solvent=solvent,
                             solute_origin=translation,
                             mixed=mixed,
                             overall_MM=SolMMM,
                             camm=parameters,
                             bsm_file=bsm_file,
                             structural_change=structural_change,
                             fderiv=fderiv,
                             redmass=ANH.redmass,
                             freq=ANH.freq,
                             ref_structure=parameters.get_pos(),
                             solpol=SolPOL,
                             gijj=ANH.K3,
                             suplist=suplist)
                
                ### evaluate corrections to the frequency shifts
                if corrections:
                   if fderiv_j: f.eval_shiftcorr(zero_camm)
                   else:        f.eval_shiftcorr(zero_camm,ua_list)

                ### nice report
                print " ===== ",typ," ===== "
                if corrections:
                   print  >> out, "%10s"%typ,
                   print  >> out,"%13.2f "*5%tuple( f.shift[0] ),
                   print  >> out,"%13.2f "*5%tuple( f.shift_corr )
                else:
                   print  >> out, "%10s"%typ,
                   print  >> out,"%13.2f "*5%tuple( f.shift[0] )
                                      
                print f
                #f.get_ShiftCorr('CH3N3_A000_D00_.camm')
                #rrr= f.get_StructuralChange(CALC.Fder,f.SOLVENT,f.solute_structure,ANH.L,PARAM.fragments)\
                #* UNITS.BohrToAngstrom
                #print
                #print f.solute_structure*UNITS.BohrToAngstrom
                #print 
                #print PARAM.fragments * UNITS.BohrToAngstrom


             if SolMMM:
                print "   ELECTRONIC+MECHANICAL "
                f = SLV(pkg="gaussian",
                               nsatoms=solvent_atno,
                               nstatoms=solute_atno,
                               L=ANH.L_,
                               mode_id=mode_id,
                               solute=parameters,
                               solute_structure=solute.get_pos(),
                               solvent=solvent,
                               solute_origin=translation,
                               mixed=mixed,
                               overall_MM=parameters,
                               bsm_file=bsm_file,
                               ref_structure=parameters.get_pos(),
                               suplist=suplist )
                        
                print Emtp.log
                #print "\n SOLUTE SOLVATOCHROMIC MULTIPOLES [A.U.]\n" 
                #print f.SOLUTE
                print  >> out, "%10s"%typ,
                print  >> out,"%13.2f "*5%tuple( f.shift[0] )
                
             if (SolCHELPG and mixed):
                print "   MIXED  "
                f = SLV(pkg="gaussian",
                               nsatoms=solvent_atno,
                               nstatoms=solute_atno,
                               L=ANH.L,
                               mode_id=mode_id,
                               solute=parameters,
                               solute_structure=solute.get_pos(),
                               solvent=solvent,
                               solute_origin=translation,
                               mixed=mixed,
                               overall_MM=parameters,
                               bsm_file=bsm_file,
                               ref_structure=parameters.get_pos(),
                               suplist=suplist )
                        
                print Emtp.log
                print "\n SOLUTE SOLVATOCHROMIC MULTIPOLES [A.U.]\n" 
                #print f.SOLUTE
             
            out.close()
              
       ## ------------------------------------------- ##
       ##             MOLECULAR DYNAMICS              ##
       ## ------------------------------------------- ##
       elif md:
          ### === RAL MONOMERS === ###
          ral_suplist = [0,3,4,5]
          #solute_atoms = [198,199,200,201,202,203]    # 36-I18C
          #solute_atoms = [236,237,238,239,240,241]    # 37-R20C
          #solute_atoms = [346,347,348,349,350,351]    # 38-N27C
          #solute_atoms = [360,361,362,363,364,365]    # 39-G28C
          #solute_atoms = [367,368,369,370,371,372]    # 40-N29C
          #solute_atoms = [398,399,400,401,402,403]    # 41-Y31C
          #solute_atoms = [419,420,421,422,423,424]    # 42-K32C
          #solute_atoms = [441,442,443,444,445,446]    # 43-S33C
          #solute_atoms = [775,776,777,778,779,780]    # 44-N54C
          suplist      = [0,1,2,3]
          #
          md_freq_shifts = SLV_MD(pkg=md_package,
                                  charges=md_charges,
                                  trajectory=md_trajectory,
                                  solute_parameters=parameters,
                                  camm=SolCAMM,
                                  suplist=suplist,
                                  ncpus=ncpus,
                                  non_atomic=False,
                                  inp=md_input)
          
          print md_freq_shifts
          shifts, averages, stds = md_freq_shifts.get_shifts()
        
          ### make a histogram!!!
          import matplotlib.pyplot as plt
          plt.hist(shifts[:,3],bins=100)
          plt.show()
# ----------------------------------------------------------------




# -------------------------------- #
#       M A I N   A C T I O N      #
# -------------------------------- #

if __name__ == '__main__':
     Main(argv[1:])

# ************************************************************************** #
# TESTS 

def sprawdz_transformacje_calkowitych_multipoli():
    normal = DMA(nfrag=1)
    point  = DMA(nfrag=1)

    normal.DMA[1][0] = array(FREQ().Dipole( argv[1] ))
    normal.DMA[2][0] = array(FREQ().Quadrupole( argv[1] ))
    normal.DMA[3][0] = array(FREQ().Octupole( argv[1] ))

    point .DMA[1][0] = array(FREQ().Dipole( argv[2] ))
    point .DMA[2][0] = array(FREQ().Quadrupole( argv[2] ))
    point .DMA[3][0] = array(FREQ().Octupole( argv[2] ))

    normal.name = 'NORMAL'
    point .name = 'POINT'

    print normal
    print point

    normal_DMA = ParseDMA( argv[1],'gaussian' )
    point_DMA  = ParseDMA( argv[2],'gaussian' )

    sup = SVDSuperimposer()
    Y = point_DMA.get_pos()
    X = normal_DMA.get_pos()
    #       rot     
    # Y <--------- X
    #               
    sup.set(Y,X)
    sup.run()
    rot, transl = sup.get_rotran()
    print sup.get_rms()

    normal.MAKE_FULL()
    normal.Rotate(rot)
    origin = array([[1.45,4.00,-11.23]]) * UNITS.AngstromToBohr
    #origin = array([[0.0,0.00,-0.0]]) * UNITS.AngstromToBohr
    normal.MAKE_FULL()
    normal.ChangeOrigin(origin)
    normal.MAKE_FULL()
    normal.DMA_FULL[2]*= UNITS.BohrElectronToDebye
    normal.DMA_FULL[3]*= UNITS.BohrElectronToDebye * UNITS.BohrToAngstrom
    normal.DMA_FULL[4]*= UNITS.BohrElectronToEsuCm * UNITS.BohrToCm**2 *1E34

    normal.MakeTraceless()
    normal.makeDMAfromFULL()

    print normal

def spwawdz_rotacje_cammow():
    from sys import argv
    normal ,a = ParseDMA(argv[1],'coulomb')
    rotate ,a = ParseDMA(argv[2],'coulomb')
    a, X = ParseDMA(argv[1][:-4]+'log','gaussian')
    a, Y = ParseDMA(argv[2][:-4]+'log','gaussian')
    del a
    sup = SVDSuperimposer()
    sup.set(Y,X)
    sup.run()
    rot, transl = sup.get_rotran()
    rotate.pos = array(Y)
    rotate.origin = array(Y)
    normal.pos = array(X)
    normal.origin = array(X)
    print "\n NORMAL\n"
    print normal
    print "\n ROTATE\n"
    print rotate
    print "\n AFTER ROTATION\n"
    print sup.get_rms()
    
    normal.MAKE_FULL()
    normal.Rotate(rot)
    
    print normal
    
    # teraz sprawdz translacje
    #print "\n SPRAWDZIAN TRANSLACJI\n"
    #print rotate
    #rotate.MAKE_FULL()
    #rotate.ChangeOrigin(zero=1)
    #rotate.MAKE_FULL()
    #rotate.ChangeOrigin(new_origin_set=Y,zero=0)
    #print "\n I CO WYSZŁO?\n"
    #print rotate   
    
    
#spwawdz_rotacje_cammow()
#sprawdz_transformacje_calkowitych_multipoli()

#from sys import argv
#f=FREQ(file=argv[1])
#f.FDeriv(Print=1,Debye=0,divide=0)
#f.Intens_1(1,0)