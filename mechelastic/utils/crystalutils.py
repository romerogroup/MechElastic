#!/usr/bin/env python

from ..tests import stability
import spglib
import numpy as np

crystallist = np.array(
    [
        "cubic",
        "hexagonal",
        "tetragonal",
        "rhombohedral-1",
        "rhombohedral-2",
        "orthorhombic",
        "monoclinic",
    ]
)

latticelist = np.array(
    [
        "hexagonal",
        "square",
        "rectangular",
        "rectangular-center",
        "oblique",
    ]
)


def crystal_select(cnew=None, cell=None, crystal_type=None, verbose=True):
    """This method selects crystal types."""

    to_print = ""
    to_print += "\n------------------------------------------------------------------\n"
    to_print += "Mechanical Stability Tests"
    to_print += "\n------------------------------------------------------------------\n"

    if crystal_type is not None:
        if verbose:
            print(to_print)
        stability.stability_test(cnew, crystal_type, verbose=verbose)

    else:
        to_print += "WARNING: crystal symmetry class  was not provided by user, it will be taken from the OUTCAR.\n"
        to_print += "One of the following was expected as the cystal_type argument: \n 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'"
        crystal_type = ""
        spg = int(spglib.get_spacegroup(cell, symprec=1e-5).split()[1][1:-1])
        if spg >= 1 and spg <= 2:
            crystal_type = "triclinic"
        if spg >= 3 and spg <= 15:
            crystal_type = "monoclinic"
        if spg >= 16 and spg <= 74:
            crystal_type = "orthorhombic"
        if spg >= 75 and spg <= 142:
            crystal_type = "tetragonal"
        if spg >= 143 and spg <= 167:
            if spg == 155 or spg == 160 or spg == 166:
                crystal_type = "rhombohedral-1"
            else:
                crystal_type = "rhombohedral-2"
        if spg >= 168 and spg <= 194:
            crystal_type = "hexagonal"
        if spg >= 195:
            crystal_type = "cubic"

        to_print += "\nFrom OUTCAR the crystal type is = %s" % crystal_type
        if verbose:
            print(to_print)
        stability.stability_test(cnew, crystal_type, verbose)


def lattice_select(cnew=None, cell=None, lattice_type=None, verbose=True):
    """This method selects crystal types."""
    to_print = ""
    if lattice_type is not None:

        to_print = ""
        to_print += (
            "\n------------------------------------------------------------------\n"
        )
        to_print += "Mechanical Stability Tests"
        to_print += (
            "\n------------------------------------------------------------------\n"
        )
        if verbose:
            print(to_print)
        stability.stability_test_2d(cnew, lattice_type, verbose)

    else:
        to_print += "\n"
        to_print += "WARNING: crystal symmetry class  was not provided by user.\n"
        to_print += "One of the following was expected as the lattice_type argument: \n 'hexagonal', 'square', 'rectangular', 'rectangular-center', 'oblique'"
        if verbose:
            print(to_print)
