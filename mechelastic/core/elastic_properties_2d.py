# -*- coding: utf-8 -*-
import numpy as np
from ..comms import printer
from ..tests import ductile
from ..tests import eigenvals
from ..utils.constants import *
from ..utils.elements import ELEMENTS
from ..utils.crystalutils import *
from .structure import Structure


class ElasticProperties2D:
    """Convert the units from GPa or N/m2 to N/m for two-dimensional systems
        1 GPa = 10^9 N/m2

        Here, we use second Piola-Kirchhoff stress method to express the 2D forces per unit length in N/m units.
        Ref: [Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013). ]
        [Singh et al., Phys. Rev. B 95, 165444 (2017)]

        We multiply elastic tensor by the thickness of the simulation cell to consider the vacuum thickness.
        In 2D:  Cij = bulk_Cij * C_latticevector (final units N/m)

        For example: if bulk_Cij = 15 GPa and out-of-plane cell parameter c = 10 Angs.
                  Then  2D_Cij = [15*10^9 N/m2] * [10*10^(-10) m] ; i.e 15*(0.1*c) N/m """

    def __init__(self, elastic_tensor, lattice_constant):
        self.elastic_tensor = elastic_tensor
        self.lattice_constant = lattice_constant
        self.c2d = np.zeros((6, 6))

        self._c2d()

    def _c2d(self):
        for i in range(0, 6):
            for j in range(0, 6):
                self.c2d[i, j] = (
                    self.elastic_tensor[i, j] * 0.1 * float(self.lattice_constant[2])
                )

        print("\n \n Elastic tensor for two-dimensional system in N/m units \n")
        np.set_printoptions(precision=3, suppress=True)
        printer.printMatrix(self.c2d)

    # Elastic Properties

    @property
    def Lm(self):
        """
        Layer modulus: represents the resistance of a 2D sheet to stretching; Lm = (1/4)*[c11 + c22 + 2*c12]  [Ref: Andrew et al.;  Phys. Rev. B 85, 125428 (2012)]

        Returns
        -------

        Lm : float
            Layer modulus

        """
        Lm = 0.25 * (self.c2d[0][0] + self.c2d[1][1] + 2 * self.c2d[0][1])

        return Lm

    @property
    def layer_modulus(self):

        """
        Returns
        -------

        layer_modulus : float
           Layer modulus

        """
        return self.Lm

    @property
    def Y10(self):
        """

        2D Young's modulus or in-plane stiffness: Y[10] = [c11c22 - c12^2]/[c22]

        Returns
        -------

        Y10 : float
            2D Young's modulus (in-plane stiffness)

        """

        Y10 = (self.c2d[0][0] * self.c2d[1][1] - self.c2d[0][1] * self.c2d[0][1]) / (
            self.c2d[1][1]
        )

        return Y10

    @property
    def Y01(self):
        """

        2D Young's modulus or in-plane stiffness: Y[01] = [c11c22 - c12^2]/[c11]

        Returns
        -------
        Y01 : float
            2D Young's modulus (in-plane stiffness)

        """

        Y01 = (self.c2d[0][0] * self.c2d[1][1] - self.c2d[0][1] * self.c2d[0][1]) / (
            self.c2d[0][0]
        )

        return Y01

    @property
    def nu10(self):
        """

        2D Poisson's ratio;  nu10 = c12/c22

        Returns
        -------
        nu10 : float
            2D Poisson's ratio

        """

        nu10 = self.c2d[0][1] / self.c2d[1][1]

        return nu10

    @property
    def nu01(self):
        """

        2D Poisson's ratio; nu01 = c12/c11

        Returns
        -------
        nu01 : float
            2D Poisson's ratio

        """

        nu01 = self.c2d[0][1] / self.c2d[0][0]
        return nu01

    @property
    def G2d(self):
        """
        2D shear modulus; G2d = C66

        Returns
        -------
        G2d : float
            2D shear modulus

        """
        G2d = self.c2d[5][5]

        return G2d

    @property
    def shear_modulus_2d(self):
        """
        2D shear modulus; G2d = C66

        Returns
        -------
        G2d : float
            2D shear modulus

        """

        return self.G2d

    def print_properties(self):

        print("\n \n Elastic properties in two-dimensions \n")

        print(
            "[Useful refs. Andrew et al.;  Phys. Rev. B 85, 125428 (2012), Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013)  ]"
        )
        print("                   ")
        print("-------------------------------------------------------")

        print("2D layer modulus (N/m)         :   %10.3f " % self.Lm)
        print("2D Young's modulus Y[10] (N/m) :   %10.3f " % self.Y10)
        print("2D Young's modulus Y[01] (N/m) :   %10.3f " % self.Y01)
        print("2D Shear modulus G (N/m)       :   %10.3f " % self.G2d)
        print("2D Poisson ratio v[10]         :   %10.3f " % self.nu10)
        print("2D Poisson ratio v[01]         :   %10.3f " % self.nu01)
        print("-------------------------------------------------------")
        print(
            "Note:  The elastic stabilty test for 2D systems is not yet implemented. "
        )

        printer.print_warning_2D()
