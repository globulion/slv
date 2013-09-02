#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Makes histogram from SLV MD output file (report.md)

Usage:

./slv_hist.py [report file] [n] [bins]

note:
n - accuracy (1: R-1,  2: R-2, ... , 5: R-5)
"""
# ------------------------------- #
from sys import argv, exit        #
from numpy import array,float64   #
from numpy import average, std    #
from scipy import optimize        #
import pylab as pl                #
import math, numpy                #
# ------------------------------- #
print __doc__                     #
if len(argv)==1: exit()           #
# ------------------------------- #

def Gaussian(n,sigma,n_o):
    """single Gaussian distribution"""
    return (1./(sigma*math.sqrt(2*math.pi)))*numpy.exp(-(n-n_o)**2/(2*sigma**2))

def ResidGaussian(p,y,n):
    """residual function"""
    sigma, n_o = p
    return y - Gaussian(n,sigma,n_o)

def r_square(func,args,data,**kwargs):
    """return R^2 coefficient of fitting"""
    data_av = numpy.average(data)
    sse = numpy.sum((func(args,**kwargs)-data)**2)
    sst = numpy.sum((data-data_av)**2)
    return 1 - sse/sst

report =     argv[1]
N      = int(argv[2])
bins   = int(argv[3])

shifts = []
data   = open(report)
for i in range(4):line = data.readline()
while line!='\n':
   if not line.startswith('#'):
      shifts.append(line.split()[N])
   line = data.readline()
data.close()
shifts = array(shifts,dtype=float64)

# ------------------------------------------------------------------------
log = '\n'
log+= ' '+28*'-'+'\n'
log+= ' Average: %6.1f ± %2.1f [cm-1]\n' % (average(shifts), std(shifts))
log+= ' '+28*'-'+'\n'
print log

n, bins, patches = pl.hist(shifts, bins,  
                           histtype='step',normed=True)
N = len(bins)-1
X = numpy.zeros(N)
X[:N] = (bins[1:N+1]+bins[:N])/2.
# initial guess for parameters
sigma0 = 5.0 ; n_o0 = 0.;
[sigma, n_o], flag = optimize.leastsq(ResidGaussian, 
                                        [sigma0, n_o0], 
                                        args=(n,X))

r2 = r_square(func=Gaussian,args=X,data=n,sigma=sigma,n_o=n_o)
func_name = '$f(\Delta\omega) = \\frac{1}{\sigma\sqrt{2\pi}} e^{-\\frac{(\Delta\omega-\Delta\omega_o)}{2\sigma^2}}$'
func_parm = 'Parameters:\n'
func_parm+= '$\Delta\omega_0=%3.2f$\n'% n_o
func_parm+= '$\sigma=%3.2f$\n'% sigma
func_parm+= '$r^2=%3.4f$\n'% r2
# plot the distribution
ax = pl.gca()
pl.title('Frequency shift distribution')
pl.xlabel('$\Delta\omega\;[{\\rm cm}^{-1}]$',fontsize=16)
pl.ylabel('$f(\Delta\omega)$',fontsize=16)
#pl.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
pl.setp(patches,color='gray')
pl.text(0.02,0.88,func_name,fontsize=21,transform=ax.transAxes,bbox=dict(facecolor='green', alpha=0.1))
pl.text(0.02,0.58,func_parm,fontsize=16,transform=ax.transAxes,bbox=dict(facecolor='blue', alpha=0.1))
pl.plot(X, Gaussian(X, sigma, n_o),linewidth=2.0)

pl.show()
# --- END --- #
