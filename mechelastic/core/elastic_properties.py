# -*- coding: utf-8 -*-
import numpy as np
from ..comms import printer
from ..tests import ductile
from ..tests import eigenvals
from ..utils.constants import *
from ..utils.elements import ELEMENTS
from ..utils.crystalutils import *
from .structure import Structure


class ElasticProperties:
    def __init__(self, elastic_tensor, structure, crystal_type):
        self.elastic_tensor = np.matrix(elastic_tensor)
        self.compaliance_tensor = self.elastic_tensor.I
        self.structure = structure
        self.crystal_type = crystal_type

        crystal_select(
            cnew=self.elastic_tensor,
            cell=self.structure.spglib_cell,
            crystal_type=self.crystal_type,
        )

    @property
    def K_v(self):
        """

        Returns
        -------
        KV : float
            Bulk Modulus Voigt.

        """
        cnew = self.elastic_tensor
        K_v = (cnew[0, 0] + cnew[1, 1] + cnew[2, 2]) + 2 * (
            cnew[0, 1] + cnew[1, 2] + cnew[2, 0]
        )
        K_v = K_v / 9.0
        return K_v

    @property
    def bulk_modulus_voigt(self):
        """


        Returns
        -------
        float
            Bulk Modulus Voigt.

        """
        return self.K_v

    @property
    def G_v(self):
        """


        Returns
        -------
        G_v : float
            Shear Modulus Voigt.

        """
        cnew = self.elastic_tensor
        ## Shear: Voigt
        G_v = (
            (cnew[0, 0] + cnew[1, 1] + cnew[2, 2])
            - (cnew[0, 1] + cnew[1, 2] + cnew[2, 0])
            + 3 * (cnew[3, 3] + cnew[4, 4] + cnew[5, 5])
        )
        G_v = G_v / 15.0
        return G_v

    @property
    def shear_modulus_voight(self):
        """


        Returns
        -------
        float
            Shear Modulus Voigt.

        """
        return self.G_v

    @property
    def E_v(self):
        """


        Returns
        -------
        float
            Young's Modulus.

        """
        return (9 * self.K_v * self.G_v) / (3 * self.K_v + self.G_v)

    @property
    def youngs_modulus_voigt(self):
        """


        Returns
        -------
        float
            Young's Modulus Voigt.

        """
        return self.E_v

    @property
    def Nu_v(self):
        """


        Returns
        -------
        float
            Poisson's Ratio Voigt.

        """
        return (3 * self.K_v - self.E_v) / (6 * self.K_v)

    @property
    def poissons_ratio_voigt(self):
        """


        Returns
        -------
        float
            Poisson's Ratio Voigt.

        """
        return self.Nu_v

    @property
    def M_v(self):
        """


        Returns
        -------
        float
            P-wave modulus voigt.

        """
        return self.K_v + (4 * self.G_v / 3.0)

    @property
    def p_wave_modulus_voigt(self):
        """


        Returns
        -------
        float
            P-wave modulus voigt.

        """
        return self.M_v

    @property
    def K_r(self):
        """


        Returns
        -------
        K_r : float
            Bulk Modulus Reuss.

        """
        snew = self.compaliance_tensor
        K_r = (snew[0, 0] + snew[1, 1] + snew[2, 2]) + 2 * (
            snew[0, 1] + snew[1, 2] + snew[2, 0]
        )
        K_r = 1.0 / K_r
        return K_r

    @property
    def bulk_modulus_reuss(self):
        """


        Returns
        -------
        float
            Bulk Modulus Reuss.

        """
        return self.K_r

    @property
    def G_r(self):
        """


        Returns
        -------
        G_r : float
            Shear Modulus Reuss.

        """
        snew = self.compaliance_tensor
        G_r = (
            4 * (snew[0, 0] + snew[1, 1] + snew[2, 2])
            - 4 * (snew[0, 1] + snew[1, 2] + snew[2, 0])
            + 3 * (snew[3, 3] + snew[4, 4] + snew[5, 5])
        )
        G_r = 15.0 / G_r
        return G_r

    @property
    def shear_modulus_reuss(self):
        """


        Returns
        -------
        G_r : float
            Shear Modulus Reuss.

        """
        return self.G_r

    @property
    def E_r(self):
        """


        Returns
        -------
        float
            Young's Modulus Reuss.

        """
        return (9 * self.K_r * self.G_r) / (3 * self.K_r + self.G_r)

    @property
    def youngs_modulus_reuss(self):
        """


        Returns
        -------
        float
            Young's modulus Reuss.

        """
        return self.E_r

    @property
    def Nu_r(self):
        """


        Returns
        -------
        float
            Poisson's Ration Reuss.

        """
        return (3 * self.K_r - self.E_r) / (6 * self.K_r)

    @property
    def poissons_ratio_reuss(self):
        """


        Returns
        -------
        float
            Poisson's Ration Reuss.

        """
        return self.Nu_r

    @property
    def M_r(self):
        """


        Returns
        -------
        float
            P-Wave Modulus Reuss.

        """
        return self.K_r + (4 * self.G_r / 3.0)

    @property
    def p_wave_modulus_reuss(self):
        """


        Returns
        -------
        float
            P-Wave Modulus Reuss.

        """
        return self.M_r

    @property
    def K_vrh(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Bulk Modulus.

        """
        return (self.K_r + self.K_v) / 2

    @property
    def bulk_modulus_voigt_reuss_hill(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Bulk Modulus.

        """
        return self.K_vrh

    @property
    def G_vrh(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Shear Modulus.

        """
        return (self.G_v + self.G_r) / 2

    @property
    def shear_modulus_voight_reuss_hill(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Shear Modulus.

        """
        return self.G_vrh

    @property
    def E_vrh(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Young's Modulus.

        """
        return (self.E_v + self.E_r) / 2

    @property
    def youngs_modulus_voigt_reuss_hill(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Young's Modulus.

        """
        return self.E_vrh

    @property
    def Nu_vrh(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Poisson's Modulus.

        """
        return (self.Nu_v + self.Nu_r) / 2.0

    @property
    def poissons_ratio_voigt_reuss_hill(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation Poisson's Modulus.

        """
        return self.Nu_vrh

    @property
    def M_vrh(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation P-wave Modulus.

        """
        return (self.M_v + self.M_r) / 2.0

    @property
    def p_wave_modulus_voigt_reuss_hill(self):
        """


        Returns
        -------
        float
            Voigt-Reuss-Hill Approximation P-wave Modulus.

        """
        return self.M_vrh

    @property
    def KG_ratio_v(self):
        """


        Returns
        -------
        float
            Bulk/Shear ration voigt.

        """
        return self.K_v / self.G_v

    @property
    def KG_ratio_r(self):
        """


        Returns
        -------
        float
            Bulk/Shear ration reuss.

        """
        return self.K_r / self.G_r

    @property
    def KG_ratio_vrh(self):
        """


        Returns
        -------
        float
            Bulk/Shear ration voigt reuss hill.

        """
        return self.K_vrh / self.G_vrh

    @property
    def A_z(self):
        """


        Returns
        -------
        float
            Zenner Anisotropy only for Cubic Crystals.

        """

        cnew = self.elastic_tensor
        return 2 * cnew[3, 3] / (cnew[0, 0] - cnew[0, 1])

    @property
    def anisotropy_zenner(self):
        """


        Returns
        -------
        float
            Zenner Anisotropy only for Cubic Crystals.

        """
        return self.anisotropy_zenner

    @property
    def A_u(self):
        """


        Returns
        -------
        float
            Ranganathan and Ostoja-Starzewski method: Phys. Rev. Lett. 101, 055504 (2008).
            for any crystalline symmetry: Universal anisotropy index.
            Note: AU is a relative measure of anisotropy with respect to a limiting value.
            For example, AU does not prove that a crystal having AU = 3 has double the anisotropy
            of another crystal with AU = 1.5. I""'

        """
        return (self.K_v / self.K_r) + 5 * (self.G_v / self.G_r) - 6.0

    @property
    def anisotropy_universal(self):
        """


        Returns
        -------
        float
            Ranganathan and Ostoja-Starzewski method: Phys. Rev. Lett. 101, 055504 (2008).
            for any crystalline symmetry: Universal anisotropy index
            Note: AU is a relative measure of anisotropy with respect to a limiting value.
            For example, AU does not prove that a crystal having AU = 3 has double the anisotropy
            of another crystal with AU = 1.5. I""'

        """
        return self.A_u

    @property
    def A_l(self):
        """


        Returns
        -------
        float
               log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016)
               AL  CV , CR   is based on the distance between the averaged stiffnesses CV and
               CR , which is more appropriate. Clearly, AL  CV , CR   is zero when the crystallite is isotropic.

        """
        return np.sqrt(5) * 2.303 * np.log(1 + (self.A_u / 5))

    @property
    def ductility(self):
        """


        Returns
        -------
        string
            The ductility of the material

        """
        return ductile.ductile_test(self.KG_ratio_vrh)

    @property
    def anisotropy_log_euclidean(self):
        """


        Returns
        -------
        float
               log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016)
               AL  CV , CR   is based on the distance between the averaged stiffnesses CV and
               CR , which is more appropriate. Clearly, AL  CV , CR   is zero when the crystallite is isotropic.

        """
        return self.A_l

    @property
    def velocity_transverse(self):
        """


        Returns
        -------
        float
            Transverse Sound Velocity(m/s) from Navier's equation using Voigt-Reuss-Hill Approximation

        """
        G = self.G_vrh * 1.0e9
        return np.sqrt((G / self.structure.density))

    @property
    def velocity_logitudinal(self):
        """


        Returns
        -------
        float
            longitudinal Sound velocity(m/s) from Navier's equation using Voigt-Reuss-Hill Approximation

        """
        G = self.G_vrh * 1.0e9  # converting from GPa to Pascal units (kg/ms^2)
        K = self.K_vrh * 1.0e9
        return np.sqrt(((3 * K + 4 * G) / (3.0 * self.structure.density)))

    @property
    def velocity_average(self):
        """


        Returns
        -------
        float
            Average Sound velocity(m/s)

        """
        vt = self.velocity_transverse
        vl = self.velocity_logitudinal

        return 1.0 / (
            np.cbrt((1.0 / 3.0) * (2.0 / (vt * vt * vt) + 1.0 / (vl * vl * vl)))
        )

    @property
    def debye_temperature(self):
        """


        Returns
        -------
        theta : float
            Debye temperature calculated using  Orson Anderson's proposal [Ref- J. Phys. Chem. Solids (1963) 24, 909-917].
            WARNING: Debye model for the atomic displacement is based on a monoatomic crystal, here we consider an average mass if your crystal has several species

        """
        q = self.structure.natoms
        density = self.structure.density
        total_mass = np.sum(self.structure.masses)
        theta = (
            (h_Planck / kB)
            * self.velocity_average
            * np.cbrt((3 * q * N_avogadro * density) / (4 * (np.pi) * total_mass))
        )
        return theta

    @property
    def melting_temperature(self):
        """


        Returns
        -------
        float
                Melting temperature estimated using empirical relation from Ref: Johnston I, Keeler G, Rollins R and Spicklemire S 1996
                Solid State Physics Simulations, The Consortium for Upper-Level Physics Software (New York: Wiley)

        """
        return 607 + 9.3 * self.K_vrh

    def print_properties(self):

        print("\n \n                         Voigt     Reuss    Average")
        print("-------------------------------------------------------")
        print(
            "Bulk modulus   (GPa)  %9.3f %9.3f %9.3f "
            % (self.K_v, self.K_r, self.K_vrh)
        )
        print(
            "Shear modulus  (GPa)  %9.3f %9.3f %9.3f "
            % (self.G_v, self.G_r, self.G_vrh)
        )
        print(
            "Young modulus  (GPa)  %9.3f %9.3f %9.3f "
            % (self.E_v, self.E_r, self.E_vrh)
        )
        print(
            "Poisson ratio         %9.3f %9.3f %9.3f "
            % (self.Nu_v, self.Nu_r, self.Nu_vrh)
        )
        print(
            "P-wave modulus  (GPa) %9.3f %9.3f %9.3f "
            % (self.M_v, self.M_r, self.M_vrh)
        )
        print(
            "Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) "
            % (self.KG_ratio_v, self.KG_ratio_r, self.KG_ratio_vrh, self.ductility)
        )
        print("-------------------------------------------------------")

        print(" \n \n \t \t Elastic Anisotropy \n ")
        print("Zener anisotropy (true for cubic crystals only) Az = %10.5f" % self.A_z)
        print(
            "Universal anisotropy index (Ranganathan and Ostoja-Starzewski method; PRL 101, 055504 (2008)) Au = %10.5f"
            % self.A_u
        )
        print(
            "Log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016) AL = %10.5f "
            % self.A_l
        )
        print(
            "\n \t \t  Elastic wave velocities calculated using Navier's equation  (in m/s units) \n"
        )
        print("----------------------------------------------- ")
        print("Longitudinal wave velocity (vl) : %10.5f " % self.velocity_logitudinal)
        print("Transverse wave velocity (vt) : %10.5f " % self.velocity_transverse)
        print("Average wave velocity (vm) : %10.5f " % self.velocity_average)
        print("Debye temperature  (in K) : %10.5f " % self.debye_temperature)
        print(
            "WARNING: Debye model for the atomic displacement is based on a monoatomic crystal, here we consider an average mass in case your crystal has several species."
        )
        # print(
        #     "Atomic displacement at 150, 300 and 450 K  (in A^2) : %10.5f %10.5f %10.5f"
        #     % (
        #         u2FromDebye(mass, natoms, theta, 150.0),
        #         u2FromDebye(mass, natoms, theta, 300.0),
        #         u2FromDebye(mass, natoms, theta, 450.0),
        #     )
        # )
        print(
            "\nMelting temperature calculated from empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K)"
        )
        print("Tm (in K)=  %10.5f (plus-minus 555 K) " % self.melting_temperature)
        print("------------------------------------------------------------------ ")

        print("Hardness analysis from 6 different semi-empirical relationships.")
        H1a, H1b, H2, H3, H4, H5 = self.hardness
        print("------------------------------------------------------------------ ")
        print("Hardness (H1a) =", "{:.2f}".format(H1a), "GPa", "  Ref.[1]")
        print("Hardness (H1b) =", "{:.2f}".format(H1b), "GPa", "  Ref.[1]")
        print("Hardness (H2)  =", "{:.2f}".format(H2), "GPa", "  Ref.[2]")
        print("Hardness (H3)  =", "{:.2f}".format(H3), "GPa", "  Ref.[3]")
        print("Hardness (H4)  =", "{:.2f}".format(H4), "GPa", "  Ref.[4]")
        print("Hardness (H5)  =", "{:.2f}".format(H5), "GPa", "  Ref.[5]")
        print("References:")
        print(
            "[1] Correlation between hardness and elastic moduli of the covalent crystals. Jiang, et al. (2011)."
        )
        print(
            "[2] Computational alchemy: the search for new superhard materials. Teter (1998)."
        )
        print(
            "[3] Mechanical and electronic properties of B12-based ternary crystals of orthorhombic phase. Jiang et al. (2010)."
        )
        print(
            "[4] Theoretical investigation on the transition-metal borides with Ta3B4-type structure: A class of hard and refractory materials. Miao et al. (2011)."
        )
        print(
            "[5] Modeling hardness of polycrystalline materials and bulk metallic glasses. Chen et al. (2011)."
        )

    @property
    def elastic_stability(self):
        return eigenvals.positive_evals(self.elastic_tensor)

    @property
    def hardness(self):
        """
        Returns
        -------
        float
            The hardness calculated by 6 different methods:
            [H1a and H1b] Correlation between hardness and elastic moduli of the covalent crystals. Jiang, et al. (2011).
            [H2] Computational alchemy: the search for new superhard materials. Teter (1998).
            [H3] Mechanical and electronic properties of B12-based ternary crystals of orthorhombic phase. Jiang et al. (2010).
            [H4] Theoretical investigation on the transition-metal borides with Ta3B4-type structure: A class of hard and refractory materials. Miao et al. (2011).
            [H5] Modeling hardness of polycrystalline materials and bulk metallic glasses. Chen et al. (2011).
        """
        B = self.K_vrh
        G = self.G_vrh
        Y = self.E_vrh
        v = self.Nu_vrh
        k = G / B
        H1a = (1 / 6.78) * G
        H1b = (1 / 16.48) * Y
        H2 = (0.1769 * G) - 2.899
        H3 = (1 / 15.76) * Y
        H4 = ((1 - 2 * v) * B) / (6 * (1 + v))
        H5 = 2 * ((k * k * G) ** 0.585) - 3
        return H1a, H1b, H2, H3, H4, H5
