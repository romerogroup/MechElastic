#!/usr/bin/env python

from mechelastic import EOS

eos_object = EOS()

infiles = ["Rm-3m", "Pm-3m", "Im-3m", "I4_mcm", "C2-m"]
natoms = [6, 1, 2, 9, 4]
eos_object.plot_enthalpy_curves(infiles, natoms, au=True)
