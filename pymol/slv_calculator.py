#!/usr/bin/env python
#*-* coding: utf-8 -*-

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2014 by Bartosz Blasiak 
#              <blasiak.bartosz@gmail.com>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
about="""
SLV Calculator, Bartosz Błasiak Copyright

This plugin application is designed to facilitate exploring
the nature of interaction-induced frequency shifts. It couples
PyMol visualization package with Solvshift programs carrying
out coarse-grained calculations of frequency shifts using Solvatochromic
Theory Based on Effective Fragment Potential Method (SolEFP). 
This includes quantum mechanical effets as follows:

a) pseudo-Coulomb electrostatics (not dumped)
b) first-order exchange-repulsion energy
c) second-order polarization effects

The charge trasfer is generally believed to be negligible. The 
dispersio effects are not developed and implemented yet.

The evaluation of solvatochromic properties is based on molecule-specific
parameters which are obtained from fully first-principles calculations
and, therefore, are completely free from any empirical manipulation.
Solvshift already offers the builtin set of various (Sol)EFP parameters
which are read when a proper molecule name is provided. The callable names
are as follows: nma, nma-d7, water, water2, meoh, dmso, mescn, meoac, dcm,
chcl3, na+, so3-- and me-so3-. water contains 5-centered DMA distribution 
for Coulomb electrostatics whereas water2 is based on 3-centered CAMM distribution.

Usage: 

1) load molecules 
2) set the appropriate selection names
3) set the appropirate (Sol)EFP fragment molecule file name
4) specify the pick-up lists used for superimposition. They contain
the information which atoms are to be superimposed with the benchmark
molecule's parameters.
"""
from pymol import cmd
import Tkinter, Pmw, webbrowser
from solvshift.solefp import EFP
from solvshift.slvpar import Frag
from numpy import float64, array, concatenate, arange
from libbbg.utilities import text_to_list
from libbbg.units import UNITS 

def __init__(self):
    self.menuBar.addmenuitem('Plugin','command','SLV Calculator',label='SLV Calculator',
                             command=lambda s=self:SLV_Calculator(s))

class Fonts:
    #helv36rn = tkFont.Font(family='Helvetica', size=36, weight='normal',
    #                       slant='roman',underline=0,overstrike=0)
    #helv36rb = tkFont.Font(family='Helvetica', size=36, weight='bold',
    #                       slant='roman',underline=0,overstrike=0)
    unknown = '-*-lucidatypewriter-medium-r-*-*-*-140-*-*-*-*-*-*'

# ------------------- PLUGIN CODE ------------------------------------------------
class SLV_Calculator(Fonts,UNITS):
    def __init__(self,app):
        self.parent = app.root
        self.createWidgets()

    def quit(self):
        self.dialog.destroy()

    def createWidgets(self): 
        # The happy dialog
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Run','Quit','Globulion','Are you happy?'),
                                 title = 'SLV Calculator',
                                 command = self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                                text = 'SLV Calculator\nBartosz Blasiak, Sep 2014 - http://www.github.com/globulion',
                                background = 'black',
                                foreground = 'white',
                                )
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
 
        # the happy notebookyyyy
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

        # Set up the Main page
        page1 = self.notebook.add('Calculator')
        group = Pmw.Group(page1,tag_text='Main options')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        # --- selections
        self.solute = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solute selection name             : ',
                                     value='',
                                        )
        self.solvent= Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solvent selection name            : ',
                                     value='',
                                        )
        self.name_1 = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solute SLV frg molecule name      : ',
                                     value='',
                                        )
        self.name_2 = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solvent SLV fragment molecule name: ',
                                     value='water',
                                        )
        self.supl_1 = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solute supl                       : ',
                                     value='None',
                                        )
        self.supl_2 = Pmw.EntryField(group.interior(),
                                     labelpos='w',
                                     label_text='Solvent supl                      : ',
                                     value='None',
                                        )
        self.checkbuttons = Pmw.RadioSelect(group.interior(),
                                            buttontype = 'checkbutton',
                                            orient = 'vertical',
                                            labelpos = 'w',
                                            selectmode='multiple',
                                            )
        for text in ('Electrostatics',
                     'Corrections',
                     'Repulsion',
                     'Polarization',):
            self.checkbuttons.add(text)
        self.checkbuttons.setvalue(['Electrostatics','Repulsion','Polarization'])

        for entry in [self.solute, self.name_1, self.supl_1, self.solvent, self.name_2, self.supl_2, self.checkbuttons]:
            entry.pack(fill='x',padx=4,pady=1)
 
        # happy results!
        t = '%7s %7s %7s %7s %7s %7s\n' % ('Coul'  .rjust(7),
                                           'Rep'   .rjust(7),
                                           'Pol'   .rjust(7),
                                           'Disp'  .rjust(7),
                                           'TOT'   .rjust(7),
                                           'rms'   .rjust(7) )
        self.resultsFrame=Tkinter.Frame(group.interior())
        bar=Tkinter.Scrollbar(self.resultsFrame,)
        self.resultsText=Tkinter.Text(self.resultsFrame,yscrollcommand=bar.set,background="#ddddff",font="Courier 12")
        bar.config(command=self.resultsText.yview)
        self.resultsText.insert(Tkinter.END,t)
        self.resultsText.pack(side=Tkinter.LEFT,expand="yes",fill="both")
        bar.pack(side=Tkinter.LEFT,expand="yes",fill="y")
        self.resultsFrame.pack(expand="yes",fill="both")

        # Set up the About page
        page2 = self.notebook.add('About Me')
        group = Pmw.Group(page2, tag_text='About PyMOL APBS Tools')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = about
        interior_frame=Tkinter.Frame(group.interior())
        bar=Tkinter.Scrollbar(interior_frame,)
        text_holder=Tkinter.Text(interior_frame,yscrollcommand=bar.set,background="#ddddff",font="Helvetica 12")
        bar.config(command=text_holder.yview)
        text_holder.insert(Tkinter.END,text)
        text_holder.pack(side=Tkinter.LEFT,expand="yes",fill="both")
        bar.pack(side=Tkinter.LEFT,expand="yes",fill="y")
        interior_frame.pack(expand="yes",fill="both")

        self.notebook.setnaturalsize()
        self.dialog.show()


    def run(self):
        xyz_1 = array(cmd.get_model(self.solute .getvalue(), 1).get_coord_list(),float64) * self.AngstromToBohr
        xyz_2 = array(cmd.get_model(self.solvent.getvalue(), 1).get_coord_list(),float64) * self.AngstromToBohr
        frg_solute    = Frag(self.name_1.getvalue())
        frg_solvent   = Frag(self.name_2.getvalue())
        xyz           =      concatenate((xyz_1,xyz_2),)

        n_solute_atoms  = frg_solute.get_natoms()
        n_solvent_atoms = frg_solvent.get_natoms()
        
        
        efp = EFP(elect=1,pol=1,rep=1,corr=1,disp=1,freq=True,mode=4,
                  ccut=None,pcut=None,ecut=None,
                  cunit=True)
        
        n_atoms         = len(xyz)
        n_solvent_mol   = (n_atoms - n_solute_atoms)/n_solvent_atoms
        
        supl = list(); reord = list()
        s_1 = self.supl_1.getvalue(); s_2 = self.supl_2.getvalue()
        if s_1.lower() == 'none': s_1 = None
        if s_2.lower() == 'none': s_2 = None
        supl.append(s_1)  # solute
        supl.append(s_2)  # solvent
        reord = [arange(n_solute_atoms), arange(n_solvent_atoms)]
        
        # --- build nmol, ind and bsm data
        bsm  = ( frg_solute, frg_solvent )
        nmol = [  frg_solute.get_natoms() , ]
        ind  = [  0                       , ]
        for i in range(n_solvent_mol):
            nmol.append( frg_solvent.get_natoms() )
            ind .append( 1                        )
        
        # --- evaluate the shifts
        efp.set(xyz,ind,nmol,bsm,supl,reord)
        efp.eval(0)
            
        # --- get the shifts
        shifts = efp.get_shifts()
        rms = efp.get_rms()
        SHIFTS = [ shifts[x] for x in ['ele_tot','rep_tot','pol_tot','dis_tot','total'] ]
        # line of output frequency shifts
        log = ''
        log+= 5*" %7.2f" % tuple(SHIFTS)
        log+= " %7.5f\n" % rms
        print log
        self.resultsText.insert(Tkinter.END,log) 
        return

    def go_to_git(self): webbrowser.open("http://www.github.com/globulion")

    def execute(self, result):
        if   result == 'Quit'          : self.quit()
        elif result == 'Run'           : self.run() 
        elif result == 'Are you happy?': print "Yes! ! ! ! ! ! ! ! ! ! :D:D:D::DD::D:D"
        elif result == 'Globulion'     : self.go_to_git()
        else: self.quit()
