#!/usr/bin/env python

"""
This contains methods to calculate elastic constants.
"""

import numpy as np
from ..comms import printer
from ..tests import ductile


def elastic_const_2D(cnew, lattconst):
    """Convert the units from GPa or N/m2 to N/m for two-dimensional systems
        1 GPa = 10^9 N/m2

        Here, we use second Piola-Kirchhoff stress method to express the 2D forces per unit length in N/m units.
        Ref: [Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013). ]
        [Singh et al., Phys. Rev. B 95, 165444 (2017)]

        We multiply elastic tensor by the thickness of the simulation cell to consider the vacuum thickness.
        In 2D:  Cij = bulk_Cij * C_latticevector (final units N/m)

        For example: if bulk_Cij = 15 GPa and out-of-plane cell parameter c = 10 Angs.
                  Then  2D_Cij = [15*10^9 N/m2] * [10*10^(-10) m] ; i.e 15*(0.1*c) N/m """

    c2d = np.zeros((6, 6))
    for i in range(0, 6):
        for j in range(0, 6):
            c2d[i][j] = cnew[i][j] * 0.1 * float(lattconst[2])

    print("\n \n Elastic tensor for two-dimensional system in N/m units \n")
    np.set_printoptions(precision=3, suppress=True)
    printer.printMatrix(c2d)

    # Elastic properties
    # Layer modulus: represents the resistance of a 2D sheet to stretching; Lm = (1/4)*[c11 + c22 + 2*c12]  [Ref: Andrew et al.;  Phys. Rev. B 85, 125428 (2012)]
    Lm = 0.25 * (c2d[0][0] + c2d[1][1] + 2 * c2d[0][1])

    # 2D Young's modulus or in-plane stiffness: Y[10] = [c11c22 - c12^2]/[c22], Y[01] = [c11c22 - c12^2]/[c11]
    Y10 = (c2d[0][0] * c2d[1][1] - c2d[0][1] * c2d[0][1]) / (c2d[1][1])
    Y01 = (c2d[0][0] * c2d[1][1] - c2d[0][1] * c2d[0][1]) / (c2d[0][0])

    # 2D Poisson's ratio;  nu10 = c12/c22, nu01 = c12/c11
    nu10 = c2d[0][1] / c2d[1][1]
    nu01 = c2d[0][1] / c2d[0][0]

    # 2D shear modulus; G2d = C66
    G2d = c2d[5][5]

    print("\n \n Elastic properties in two-dimensions \n")
    print(
        "[Useful refs. Andrew et al.;  Phys. Rev. B 85, 125428 (2012), Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013)  ]"
    )
    print("                   ")
    print("-------------------------------------------------------")
    print("2D layer modulus (N/m)         :   %10.3f " % Lm)
    print("2D Young's modulus Y[10] (N/m) :   %10.3f " % Y10)
    print("2D Young's modulus Y[01] (N/m) :   %10.3f " % Y01)
    print("2D Shear modulus G (N/m)       :   %10.3f " % G2d)
    print("2D Poisson ratio v[10]         :   %10.3f " % nu10)
    print("2D Poisson ratio v[01]         :   %10.3f " % nu01)
    print("-------------------------------------------------------")
    print("Note:  The elastic stabilty test for 2D systems is not yet implemented. ")

    printer.print_warning_2D()

    return
