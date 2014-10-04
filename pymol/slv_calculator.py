# -*- coding: utf-8 -*-
# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2004 by Charles Moad <cmoad@indiana.edu>
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

import tkSimpleDialog
import tkMessageBox
from pymol import cmd
import sys, urllib2, zlib, copy, units
import Tkinter, Pmw, webbrowser
from solvshift.solefp import EFP
from solvshift.slvpar import Frag
from numpy import float64, array, concatenate
from utilities import text_to_list

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
class SLV_Calculator(Fonts,units.UNITS):
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
                                text = 'SLV Calculator\nBartosz Błasiak, Sep 2014 - http://www.github.com/globulion',
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
        self.resultsFrame=Tkinter.Frame(group.interior())
        bar=Tkinter.Scrollbar(self.resultsFrame,)
        self.resultsText=Tkinter.Text(self.resultsFrame,yscrollcommand=bar.set,background="#ddddff",font="Helvetica 12")
        bar.config(command=self.resultsText.yview)
        self.resultsText.insert(Tkinter.END,'')
        self.resultsText.pack(side=Tkinter.LEFT,expand="yes",fill="both")
        bar.pack(side=Tkinter.LEFT,expand="yes",fill="y")
        self.resultsFrame.pack(expand="yes",fill="both")

        # Set up the About page
        page2 = self.notebook.add('About Me')
        group = Pmw.Group(page2, tag_text='About PyMOL APBS Tools')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = """HaaaaaaaaaaaaaAAAAAaaaaaaaaaaaaaPPPY!"""
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
        
        t = '%10s %10s %10s %10s %10s %10s %10s %10s %10s\n' % ('EL-TOT' .rjust(10),
                                                              'EL-MEA' .rjust(10),
                                                              'EL-EA'  .rjust(10),
                                                              'POL'    .rjust(10),
                                                              'REP'    .rjust(10),
                                                              'EL-MEA-C'.rjust(10),
                                                              'EL-EA-C' .rjust(10),
                                                              'TOTAL'   .rjust(10),
                                                              'RMS'     .rjust(10) )
        
        efp = EFP(elect=1,pol=1,rep=1,corr=0,freq=True,
                  ccut=None,pcut=None,ecut=None,
                  cunit=True)
        
        n_atoms         = len(xyz)
        n_solvent_mol   = (n_atoms - n_solute_atoms)/n_solvent_atoms
        
        supl = list()
        s_1 = self.supl_1.getvalue(); s_2 = self.supl_2.getvalue()
        if s_1.lower() == 'none': s_1 = None
        if s_2.lower() == 'none': s_2 = None
        supl.append(s_1)  # solute
        supl.append(s_2)  # solvent
        
        # --- build nmol, ind and bsm data
        nmol = [  frg_solute.get_natoms() , ]
        ind  = [  0                       , ]
        for i in range(n_solvent_mol):
            nmol.append( frg_solvent.get_natoms() )
            ind .append( 1                        )
            bsm  = ( frg_solute, frg_solvent )
        
        # --- evaluate the shifts
        efp.set(xyz,ind,nmol,bsm,supl)
        efp.eval(0)
            
        # --- get the shifts
        shifts = efp.get_shift()
        rms = efp.get_rms()
        # line of output frequency shifts
        log = t
        log+= "%10.2f"    % (shifts[0]+shifts[1])
        log+= 7*" %10.2f" % tuple(shifts[:7])
        log+= " %10.5f\n" % rms
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
