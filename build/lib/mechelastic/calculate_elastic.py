#!/usr/bin/env python


from .comms import printer
from .parsers import VaspOutcar
from .parsers import AbinitOutput
from .parsers import QE_ElaStic_Parser
from .parsers import QE_thermo_pw_Parser
from .core import ElasticProperties
from .core import ElasticProperties2D
import os, sys


def calculate_elastic(
    infile="OUTCAR",
    dim="3D",
    elastic_tensor=None,
    crystal=None,
    lattice_type=None,
    code="vasp",
    anaddbfile=None,
    qe_outfile=None,
    adjust_pressure=True,
    verbose=True,
    outfile="elastic_properties.txt",
):
    """
    Calculate the elastic properties of a material from a DFT calculation.

    Parameters
    ----------
    infile : str, optional
        Path to the input file which is a DFT calculation outputfile. The default is "OUTCAR".
    dim : str, optional
        Dimension of the structure. The default is "3D".
    crystal : str, optional
        Crystal family (only used in 3D). The default is None.
    elastic_tensor : float, optional
        The elastic tensor. This option is useful if one does not want to use a DFT output.
        The default is None.
    lattice_type : TYPE, optional
        2D lattice type. The default is None.
    code : str, optional
        DFT code used to generate the outputs. The default is "vasp".
    anaddbfile : str, optional
        Path to the DDB file (applicable only in abinit). The default is None.
    qe_outfile : str, optional
        Path to the Quantum Espresso output file. The default is None.
    adjust_pressure : bool, optional
        To adjust the cell pressure according to the output file. The default is True.
    verbose : str, optional
        To print the progress of the elastic calculations. The default is True.
    outfile : TYPE, optional
        Path to the desired output file. Acceptable files are JSON, XML, TXT.
        The default is "elastic_properties.txt".

    Returns
    -------
    elastic_properties : TYPE
        DESCRIPTION.

    """
    # Check if infile is present
    if not os.path.exists(infile):
        print("%s doesn't exist. Exiting." % infile)
        sys.exit()

    # welcome message
    if verbose:
        printer.print_mechelastic()
        if code != None:
            print("\nThis matrix was computed from " + code)
    elastic_tensor = elastic_tensor
    structure = None
    lattice_constant = None
    crystal_type = crystal

    # calling parser
    if code == "vasp" and elastic_tensor is None:
        output = VaspOutcar(
            infile=infile, adjust_pressure=adjust_pressure, verbose=verbose
        )
        elastic_tensor = output.elastic_tensor
        structure = output.structure
        lattice_constant = output.lattice_constant

    elif code == "abinit":
        output = AbinitOutput(infile=infile, anaddbfile=anaddbfile)
        elastic_tensor = output.elastic_tensor
        structure = output.structure
        lattice_constant = output.lattice_constant

    elif code == "qe_ElaStic":
        output = QE_ElaStic_Parser(outfile=qe_outfile, infile=infile)
        elastic_tensor = output.elastic_tensor
        structure = output.structure
        lattice_constant = output.lattice_constant

    elif code == "qe_thermo_pw":
        output = QE_thermo_pw_Parser(outfile=qe_outfile, infile=infile)
        elastic_tensor = output.elastic_tensor
        structure = output.structure
        lattice_constant = output.lattice_constant

    # elastic constants calculation for 3D materials
    if dim == "3D":
        elastic_properties = ElasticProperties(
            elastic_tensor, structure, crystal_type, verbose=verbose
        )
        if verbose:
            print(elastic_properties)
        if "json" in outfile:
            elastic_properties.to_json(outfile)
        elif "xml" in outfile:
            elastic_properties.to_xml(outfile)
        else:
            elastic_properties.to_file(outfile)

    # elastic constants calculation for 2D materials
    elif dim == "2D":
        elastic_properties = ElasticProperties2D(
            elastic_tensor, lattice_constant, lattice_type=lattice_type
        )
        elastic_properties.print_properties()

    if verbose:
        print("\nThanks! See you later. ")
    return elastic_properties
