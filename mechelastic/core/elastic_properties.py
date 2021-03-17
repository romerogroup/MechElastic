# -*- coding: utf-8 -*-
import numpy as np
from ..comms import printer
from ..tests import ductile
from ..tests import eigenvals
from ..utils.constants import *
from ..utils.elements import ELEMENTS
from ..utils.crystalutils import *
from .structure import Structure
from numpy import linalg as LA


class ElasticProperties:
    def __init__(self, elastic_tensor, structure=None, crystal_type=None, code=None):
        """


        Parameters
        ----------
        elastic_tensor : TYPE
            DESCRIPTION.
        structure : TYPE, optional
            DESCRIPTION. The default is None.
        crystal_type : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        self.elastic_tensor = np.matrix(elastic_tensor)
        self.compaliance_tensor = self.elastic_tensor.I
        self.structure = structure
        self.crystal_type = crystal_type
        print("\n------------------------------------------------------------------")
        print("Elastic Tensor (GPa units)")
        print("------------------------------------------------------------------\n")

        printMatrix(self.elastic_tensor)
        if code != None:
            print("\nThis matrix was computed from " + code)
        print("Note: For VASP users this is the pressure-corrected matrix")

        print("\n------------------------------------------------------------------")
        print("Elastic Tensor Eigen Values (GPa units)")
        print("------------------------------------------------------------------\n")

        positive_evals(self.elastic_tensor)

        if crystal_type is not None or structure is not None:
            print(
                "\n------------------------------------------------------------------"
            )
            print("Mechanical Stability Tests")
            print(
                "------------------------------------------------------------------\n"
            )
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
        return self.A_z

    @property
    def A_cb(self):
        """


        Returns
        -------
        float
            Chung-Buessem only for Cubic Crystals.

        """

        a_cb = (self.G_v - self.G_r) / (self.G_v + self.G_r)
        return a_cb

    @property
    def anisotropy_Chung_Buessem(self):
        """


        Returns
        -------
        float
            Chung-Buessem  only for Cubic Crystals.

        """
        return self.A_cb

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
        return np.sqrt( ( np.log(self.K_v/self.K_r) )**2 + 5 * ( np.log(self.G_v/self.G_r) )**2  )
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
        if self.structure is None:
            raise Exception("No structure was provided")

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
        if self.structure is None:
            raise Exception("No structure was provided")
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
        if self.structure is None:
            raise Exception("No structure was provided")
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
        print("\n------------------------------------------------------------------")
        print("Elastic Moduli")
        print("------------------------------------------------------------------\n")

        print("                          Voigt     Reuss    Average")
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
            "Young's modulus  (GPa)  %9.3f %9.3f %9.3f "
            % (self.E_v, self.E_r, self.E_vrh)
        )
        print(
            "Poisson's ratio         %9.3f %9.3f %9.3f "
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

        print("\n------------------------------------------------------------------")
        print("Elastic parameters")
        print("------------------------------------------------------------------\n")

        print("Lame's first and second parameter; Ref.[3]")
        print("Lambda (GPa)  =  %10.3f  " % self.lambda_lame_coefficient)
        print("Mu (GPa)  =  %10.3f  " % self.mu_lame_coefficient)

        print("\n------------------------------------------------------------------")
        print("Bonding information")
        print("------------------------------------------------------------------\n ")

        print("Kleinman’s parameter; Ref.[4,5]")
        print("NOTE: K = 0 (1) bending (stretching) would dominate")
        print("K =  %10.5f  " % self.kleinman_parameter)
        print("")
        print("Cauchy's Pressure calculated from the relation : CP = C_12 - C_44")
        print("     CP > 0 (+ve) indicates that ionic bonding dominates")
        print("     CP < 0 (-ve) indicates that covalent bonding dominates")

        print("CP (GPa) =  %10.3f  " % self.cauchy_pressure)
        print("Bonding is mainly " + self.bonding_type)

        print("\n------------------------------------------------------------------")
        print("Elastic Anisotropy")
        print("------------------------------------------------------------------\n")

        print(
            "Zener's anisotropy (true for cubic crystals only); Az = %6.3f; Ref.[6]"
            % self.A_z
        )
        print(
            "Chung-Buessem's anisotropy (true for cubic crystals only); Acb = %6.3f; Ref.[7]"
            % self.A_cb
        )
        print("Universal anisotropy index; Au = %6.3f; Ref.[8]" % self.A_u)
        print("Log-Euclidean's anisotropy; AL = %6.3f; Ref.[9]" % self.A_l)

        print("\n------------------------------------------------------------------")
        print("Elastic Wave Velocities and Debye Temperature")
        print("------------------------------------------------------------------\n")

        print(
            "Longitudinal wave velocity (vl) : %10.3f m/s; Ref.[10]"
            % self.velocity_logitudinal,
        )
        print(
            "Transverse wave velocity (vt) : %10.3f m/s; Ref.[10]"
            % self.velocity_transverse,
        )
        print(
            "Average wave velocity (vm) : %10.3f m/s; Ref.[10]" % self.velocity_average,
        )
        print("Debye temperature  : %10.3f K; Ref.[10]" % self.debye_temperature)
        print("")
        print(
            "WARNING: The  Debye model for the atomic displacement is based on a monoatomic crystal approximation.\n Here we consider an averaged mass, in case your crystal has several species."
        )

        # print(
        #     "Atomic displacement at 150, 300 and 450 K  (in A^2) : %10.5f %10.5f %10.5f"
        #     % (
        #         u2FromDebye(mass, natoms, theta, 150.0),
        #         u2FromDebye(mass, natoms, theta, 300.0),
        #         u2FromDebye(mass, natoms, theta, 450.0),
        #     )
        # )
        print("\n------------------------------------------------------------------")
        print("Melting Temperature")
        print("------------------------------------------------------------------\n")

        print(
            "Melting temperature calculated from the empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K); Ref.[11]"
        )
        print("Tm =  %10.3f K (plus-minus 555 K) " % self.melting_temperature)
        print("\n")
        print(
            "WARNING: This is a crude approximation and its validity needs to be checked! "
        )

        print("\n------------------------------------------------------------------")
        print("Hardness Analysis")
        print("------------------------------------------------------------------\n")

        H1a, H1b, H2, H3, H4, H5 = self.hardness
        print("Hardness (H1a) =", "{:.2f}".format(H1a), "GPa;", "  Ref.[12]")
        print("Hardness (H1b) =", "{:.2f}".format(H1b), "GPa;", "  Ref.[12]")
        print("Hardness (H2)  =", "{:.2f}".format(H2), "GPa;", "  Ref.[13]")
        print("Hardness (H3)  =", "{:.2f}".format(H3), "GPa;", "  Ref.[14]")
        print("Hardness (H4)  =", "{:.2f}".format(H4), "GPa;", "  Ref.[15]")
        print("Hardness (H5)  =", "{:.2f}".format(H5), "GPa;", "  Ref.[16]")
        print("")
        print(
            """Hardness recommendation model:
        ********************************************************************
                     Cubic  Hexagonal  Orthorhombic  Rhombohedral  General
        ********************************************************************
        Insulator      H2      H1b          H2            H2          H2
        Semiconductor  H5    H1b, H3                      H2          H5
        Metal          H1a     H4           H4            H4          H4
        ********************************************************************
        Insulator     : bandgap > 2 eV
        Semiconductor : bandgap < 2 eV
        Metal         : bandgap = 0 """
        )

        print("\n------------------------------------------------------------------")
        print("References")
        print("------------------------------------------------------------------\n")
        print(
            "[1] Necessary and Sufficient Elastic Stability Conditions in Various Crystal Systems. Félix Mouhat and François-Xavier Coudert. Phys. Rev. B (2014)"
        )
        print(
            "[2] Crystal Structures and Elastic Properties of Superhard IrN2 and IrN3 from First Principles. Zhi-jian Wu et al. Phys. Rev. B (2007)"
        )
        print(
            "[3]  The rock physics handbook, Cam-bridge university press. G. Mavko, T. Mukerji, J. Dvorkin. Cambridge University Press. (2020)"
        )
        print(
            "[4]  Deformation potentials in silicon. i. uniaxial strain. Leonard Kleinman. Phys.Rev. 128. (1962)"
        )
        print(
            "[5] Electronic Structure and the Properties of Solids: The Physics of the Chemical Bond. Walter A. Harrison.  (2012)"
        )
        print(
            "[6] Elasticity and Anelasticity of Metals. Clarence M. Zener et al. The Journal of Physical Chemistry. (1949)"
        )
        print(
            "[7] The Elastic Anisotropy of Crystals. D. H. Chung and W. R. Buessem. Journal of Applied Physics (1967)"
        )
        print(
            "[8] Universal Elastic Anisotropy Index. Shivakumar I. Ranganathan et al. Phys. Rev. Lett. (2008)"
        )
        print(
            "[9] Elastic Anisotropy of Crystals. Christopher M. Kube. AIP Advances. (2016)"
        )
        print(
            "[10] A Simplified Method for Calculating the Debye Temperature from Elastic Constants. Orson L.Anderson. Journal of Physics and Chemistry of Solids. (1963)"
        )
        print(
            "[11] Elastic constants versus melting temperature in metals. M.E.Fine et al. Scripta Metallurgica. (1984)"
        )
        print(
            "[12] Correlation between hardness and elastic moduli of the covalent crystals. Jiang, et al. (2011)."
        )
        print(
            "[13] Computational alchemy: the search for new superhard materials. Teter (1998)."
        )
        print(
            "[14] Mechanical and electronic properties of B12-based ternary crystals of orthorhombic phase. Jiang et al. (2010)."
        )
        print(
            "[15] Theoretical investigation on the transition-metal borides with Ta3B4-type structure: A class of hard and refractory materials. Miao et al. (2011)."
        )
        print(
            "[16] Modeling hardness of polycrystalline materials and bulk metallic glasses. Chen et al. (2011)."
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

    @property
    def cauchy_pressure(self):
        """
        This parameter desceibes the nature of bonding
        CP > 0 (+ve) indicates that ionic bonding dominates
        CP < 0 (-ve) indicates that covalent bonding dominates
        Returns
        -------
        None.

        """
        return self.elastic_tensor[0, 1] - self.elastic_tensor[3, 3]

    @property
    def bonding_type(self):
        """
        This parameter desceibes the nature of bonding
        CP > 0 (+ve) indicates that ionic bonding dominates
        CP < 0 (-ve) indicates that covalent bonding dominates

        Returns
        -------
        str
            DESCRIPTION.

        """
        cauchy_pressure = self.cauchy_pressure
        if cauchy_pressure > 0:
            return "ionic"
        elif cauchy_pressure < 0:
            return "covalent"

    @property
    def kleinman_parameter(self):
        c = self.elastic_tensor
        return (c[0, 0] + 8 * c[0, 1]) / (7 * c[0, 0] - 2 * c[0, 1])

    @property
    def bond_bending_vs_streching(self):
        return

    @property
    def lambda_lame_coefficient(self):
        """


        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.E_vrh* self.Nu_vrh/ ((1 + self.Nu_vrh) * (1 - 2 * self.Nu_vrh))

    @property
    def mu_lame_coefficient(self):
        return self.E_vrh / (2 * (1 + self.Nu_vrh))


def positive_evals(cnew):
    """This method checks the postivity of the eigenvalues
    of a matrix."""

    # print("Eigen Values of the matrix:")
    evals = list(LA.eigvals(cnew))
    evalsPrint = list(np.around(np.array(evals), 3))
    print("%s" % evalsPrint)
    check = 0
    for i in range(len(evals)):
        if evals[i] > 0.0:
            pass
        else:
            check = 1
    # if check == 1:
    #     print(
    #         "ATTENTION: One or more eigen values are negative indicating elastic instability."
    #     )
    # if check == 0:
    #     print("All eigen values are positive indicating elastic stability.")

    return not (bool(check))


def printMatrix(c):
    row = c.shape[0]
    col = c.shape[1]
    for i in range(row):
        for j in range(col):
            print("{:>10.3f} ".format(c[i, j]), end=" ")
            if j == (col - 1):
                print(" ")

    return
