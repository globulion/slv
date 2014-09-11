#!/usr/bin/python
"""
Shows the FRG file

Usage:
  ./show.py frg_file
"""
print __doc__
from sys import argv,exit
if len(argv)==1:exit()
from solvshift.slvpar import Frag

a = Frag(argv[1])
par = a.get()

# print LMOC
lmoc = par['lmoc']
print " LMO Centrioids:\n"
for i in lmoc:
    print " %10.4f %10.4f %10.4f" % tuple(i)
print
