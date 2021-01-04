#!/usr/bin/env python

from mechelastic import EOS

eos_object = EOS()

infiles = ["Rm-3m", "Pm-3m", "Im-3m", "I4_mcm", "C12:m1"]
natoms = [6, 1, 2, 8, 4]
eos_object.plot_enthalpy_curves(infiles, natoms, au=True)
