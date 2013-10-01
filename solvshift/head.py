# --------------------- #
#  SOLVSHIFT HEAD FILE  #
# --------------------- #

# ---- import Python built-in modules
from numpy  import *
from sys    import argv, exit
from string import Template
import re, os, getopt, sys, pylab
reflags = re.DOTALL

# ---- import BBG library
from dma        import *
from units      import *
from utilities  import *
from utilities2 import *
from gaussfreq  import *

# ---- import Coulomb.py and PyQuante packages
if os.environ['__IMPORT__COULOMB__']:
   from coulomb  import work
   from coulomb  import *
   from coulomb.cube import *
   from PyQuante import *
   
# ---- import Solvshift modules
from ff         import *
from diff       import *
from inputs     import *
from mcho       import *
from slv        import *
from dipderiv   import *
from md         import *
from hessian    import HESSIAN
from slvpar     import SLVPAR
from solefp     import SOLEFP
