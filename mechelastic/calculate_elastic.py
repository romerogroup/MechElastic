#!/usr/bin/env python

from .comms import printer
from .core import elastic_bulk
from .core import elastic_2D
from .parsers import VaspOutcar


def calculate_elastic(infile="OUTCAR", dim="3D", crystal=None, code="vasp"):

    """
    This method calculates the elastic properties
    of a material from a DFT calculation.
    """

    # welcome message
    printer.print_mechelastic()

    # calling parser
    if code == "vasp":
        output = VaspOutcar(infile=infile, crystal_type=crystal)
        # parserclass = VASPParser(infile=infile)
        # cnew = parserclass.cnew
        # snew = parserclass.snew
        # cell = parserclass.cell
        # density = parserclass.density
        # natoms = parserclass.natoms
        # totalmass = parserclass.totalmass
        # lattconst = parserclass.lattconst

    # elastic constants calculation for 3D materials
    if dim == "3D":
        # elastic_bulk.elastic_const_bulk(
        #     cnew, snew, crystal, cell, density, natoms, totalmass
        # )
        output.elastic_properties.print_properties()

    # elastic constants calculation for 2D materials
    # elif dim == "2D": TODO
    #     elastic_2D.elastic_const_2D(cnew, lattconst)

    # other
    # else: We don't need this
    #     elastic_bulk.elastic_const_bulk(
    #         cnew, snew, crystal, cell, density, natoms, totalmass
    #     )

    print("\nThanks! See you later. ")
    return output
