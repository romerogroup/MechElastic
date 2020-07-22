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


def crystalselect(cnew=None, cell=None, crystaltype=None):

    """This method selects crystal types."""

    if crystaltype is not None:
        print("\n \t \t Mechanical Stability Test \n")
        stability.stability_test(cnew, crystaltype)

    else:
        print("\n")
        print(
            "WARNING: crystal symmetry class  was not provided by user, it will be taken from the OUTCAR"
        )
        print(
            "One of the following was expected as the second argument: \n 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'"
        )
        crystaltype = ""
        spg = int(spglib.get_spacegroup(cell, symprec=1e-5).split()[1][1:-1])
        if spg >= 1 and spg <= 2:
            crystaltype = "triclinic"
        if spg >= 3 and spg <= 15:
            crystaltype = "monoclinic"
        if spg >= 16 and spg <= 74:
            crystaltype = "orthorhombic"
        if spg >= 75 and spg <= 142:
            crystaltype = "tetragonal"
        if spg >= 143 and spg <= 167:
            if spg == 155 or spg == 160 or spg == 166:
                crystaltype = "rhombohedral-1"
            else:
                crystaltype = "rhombohedral-2"
        if spg >= 168 and spg <= 194:
            crystaltype = "hexagonal"
        if spg >= 195:
            crystaltype = "cubic"

        print("From OUTCAR the crystal type is = ", crystaltype)
        stability.stability_test(cnew, crystaltype)
