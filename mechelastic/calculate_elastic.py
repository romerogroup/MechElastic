#!/usr/bin/env python

from .comms import printer
from .parsers import VaspOutcar
from .core import ElasticProperties
from .core import ElasticProperties2D


def calculate_elastic(infile="OUTCAR", dim="3D", crystal=None, code="vasp"):

    """
    This method calculates the elastic properties
    of a material from a DFT calculation.
    """

    # welcome message
    printer.print_mechelastic()

    elastic_tensor = None
    structure = None
    lattice_constant = None
    crystal_type = crystal

    # calling parser
    if code == "vasp":
        output = VaspOutcar(infile=infile)
        elastic_tensor = output.elastic_tensor
        structure = output.structure
        lattice_constant = output.lattice_constant

    # elastic constants calculation for 3D materials
    if dim == "3D":
        elastic_properties = ElasticProperties(elastic_tensor, structure, crystal_type)
        elastic_properties.print_properties()

    # elastic constants calculation for 2D materials
    elif dim == "2D":
        elastic_properties = ElasticProperties2D(elastic_tensor, lattice_constant)
        elastic_properties.print_properties()

    # other
    # else: We don't need this
    #     elastic_bulk.elastic_const_bulk(
    #         cnew, snew, crystal, cell, density, natoms, totalmass
    #     )

    print("\nThanks! See you later. ")
    return output
