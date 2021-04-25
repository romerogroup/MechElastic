# -*- coding: utf-8 -*-
import numpy as np
import json
from dicttoxml import dicttoxml
from xml.dom.minidom import parseString
from ..comms import printer
from ..tests import ductile
from ..tests import eigenvals
from ..utils.constants import *
from ..utils.elements import ELEMENTS
from ..utils.crystalutils import *
from .structure import Structure


class ElasticProperties:
    def __init__(self, elastic_tensor, structure=None, crystal_type=None, verbose=True):
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
        self.verbose = verbose
        # if verbose:
        #     print_matrix(self.elastic_tensor) alredy printed

        eigenvals.positive_evals(self.elastic_tensor, verbose)

        if crystal_type is not None or structure is not None:

            crystal_select(
                cnew=self.elastic_tensor,
                cell=self.structure.spglib_cell,
                crystal_type=self.crystal_type,
                verbose=verbose,
            )
        if verbose:
            print(self.__str__)

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
            Bulk/Shear ratio voigt.

        """
        return self.K_v / self.G_v

    @property
    def KG_ratio_r(self):
        """


        Returns
        -------
        float
            Bulk/Shear ratio reuss.

        """
        return self.K_r / self.G_r

    @property
    def KG_ratio_vrh(self):
        """


        Returns
        -------
        float
            Bulk/Shear ratio voigt reuss hill.

        """
        return self.K_vrh / self.G_vrh

    @property
    def bulk_shear_ratio_voigt(self):
        """


        Returns
        -------
        float
            Bulk/Shear ratio voigt.

        """
        return self.K_v / self.G_v

    @property
    def bulk_shear_ratio_reuss(self):
        """


        Returns
        -------
        float
            Bulk/Shear ratio reuss.

        """
        return self.K_r / self.G_r

    @property
    def bulk_shear_ratio_voigt_reuss_hill(self):
        """


        Returns
        -------
        float
            Bulk/Shear ratio voigt reuss hill.

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
        return np.sqrt(
            (np.log(self.K_v / self.K_r)) ** 2 + 5 * (np.log(self.G_v / self.G_r)) ** 2
        )

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
        return self.E_vrh * self.Nu_vrh / ((1 + self.Nu_vrh) * (1 - 2 * self.Nu_vrh))

    @property
    def mu_lame_coefficient(self):
        return self.E_vrh / (2 * (1 + self.Nu_vrh))

    def to_dict(self, symprec=1e-5):
        """


        Parameters
        ----------
        symprec : float
            Precision used in calculating the space group in angstroms. The default is 1e-5.

        Returns
        -------
        dict
            DESCRIPTION.

        """
        return {
            "anisotropy_Chung_Buessem": self.anisotropy_Chung_Buessem,
            "anisotropy_log_euclidean": self.anisotropy_log_euclidean,
            "anisotropy_universal": self.anisotropy_universal,
            "anisotropy_zenner": self.anisotropy_zenner,
            "bond_bending_vs_streching": self.bond_bending_vs_streching,
            "bonding_type": self.bonding_type,
            "bulk_modulus_reuss": self.bulk_modulus_reuss,
            "bulk_modulus_voigt": self.bulk_modulus_voigt,
            "bulk_modulus_voigt_reuss_hill": self.bulk_modulus_voigt_reuss_hill,
            "bulk_shear_ratio_reuss": self.bulk_shear_ratio_reuss,
            "bulk_shear_ratio_voigt": self.bulk_shear_ratio_voigt,
            "bulk_shear_ratio_voigt_reuss_hill": self.bulk_shear_ratio_voigt_reuss_hill,
            "cauchy_pressure": self.cauchy_pressure,
            "compaliance_tensor": self.compaliance_tensor.tolist(),
            "crystal_type": self.crystal_type,
            "debye_temperature": self.debye_temperature,
            "ductility": self.ductility,
            "elastic_stability": self.elastic_stability,
            "elastic_tensor": self.elastic_tensor.tolist(),
            "hardness": self.hardness,
            "kleinman_parameter": self.kleinman_parameter,
            "lambda_lame_coefficient": self.lambda_lame_coefficient,
            "melting_temperature": self.melting_temperature,
            "mu_lame_coefficient": self.mu_lame_coefficient,
            "p_wave_modulus_reuss": self.p_wave_modulus_reuss,
            "p_wave_modulus_voigt": self.p_wave_modulus_voigt,
            "p_wave_modulus_voigt_reuss_hill": self.p_wave_modulus_voigt_reuss_hill,
            "poissons_ratio_reuss": self.poissons_ratio_reuss,
            "poissons_ratio_voigt": self.poissons_ratio_voigt,
            "poissons_ratio_voigt_reuss_hill": self.poissons_ratio_voigt_reuss_hill,
            "shear_modulus_reuss": self.shear_modulus_reuss,
            "shear_modulus_voight": self.shear_modulus_voight,
            "shear_modulus_voight_reuss_hill": self.shear_modulus_voight_reuss_hill,
            "structure": self.structure.to_dict(symprec),
            "velocity_average": self.velocity_average,
            "velocity_logitudinal": self.velocity_logitudinal,
            "velocity_transverse": self.velocity_transverse,
            "youngs_modulus_reuss": self.youngs_modulus_reuss,
            "youngs_modulus_voigt": self.youngs_modulus_voigt,
            "youngs_modulus_voigt_reuss_hill": self.youngs_modulus_voigt_reuss_hill,
        }

    def to_json(self, outfile="elastic_properties.json", symprec=1e-5):
        """


        Parameters
        ----------
        outfile : str, optional
            file path to the output file. The default is "elastic_properties.json".
        symprec : float, optional
            Precision used in calculating the space group in angstroms. The default is 1e-5.

        Returns
        -------
        None.

        """
        wf = open(outfile, "w")
        json.dump(
            self.to_dict(symprec), wf, sort_keys=True, indent=4, separators=(",", ": ")
        )
        wf.close()

    def to_xml(self, outfile="elastic_properties.xml", symprec=1e-5):
        """


        Parameters
        ----------
        outfile : str, optional
            file path to the output file. The default is "elastic_properties.xml".
        symprec : float, optional
            Precision used in calculating the space group in angstroms. The default is 1e-5.

        Returns
        -------
        None.

        """
        wf = open(outfile, "w")
        xml = dicttoxml(self.to_dict(symprec))
        dom = parseString(xml)
        wf.write(dom.toprettyxml())
        wf.close()

    def to_file(self, outfile="elastic_properties.txt"):
        """


        Parameters
        ----------
        outfile : str, optional
            Path to the output file. The default is "elastic_properties.txt".

        Returns
        -------
        None.

        """
        wf = open(outfile, "w")
        wf.write(self.__str__())
        wf.close()

    def __str__(self):
        """


        Returns
        -------
        None.

        """

        ret = ""
        ret += "\n------------------------------------------------------------------\n"
        ret += "Elastic Moduli\n"
        ret += "------------------------------------------------------------------\n\n"

        ret += "                          Voigt     Reuss    Average\n"
        ret += "-------------------------------------------------------\n"
        ret += "Bulk modulus   (GPa)  %9.3f %9.3f %9.3f \n" % (
            self.K_v,
            self.K_r,
            self.K_vrh,
        )
        ret += "Shear modulus  (GPa)  %9.3f %9.3f %9.3f \n" % (
            self.G_v,
            self.G_r,
            self.G_vrh,
        )
        ret += "Young's modulus  (GPa)  %9.3f %9.3f %9.3f \n" % (
            self.E_v,
            self.E_r,
            self.E_vrh,
        )
        ret += "Poisson's ratio         %9.3f %9.3f %9.3f \n" % (
            self.Nu_v,
            self.Nu_r,
            self.Nu_vrh,
        )
        ret += "P-wave modulus  (GPa) %9.3f %9.3f %9.3f \n" % (
            self.M_v,
            self.M_r,
            self.M_vrh,
        )
        ret += "Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) \n" % (
            self.KG_ratio_v,
            self.KG_ratio_r,
            self.KG_ratio_vrh,
            self.ductility,
        )

        ret += "\n------------------------------------------------------------------\n"
        ret += "Elastic parameters\n"
        ret += "------------------------------------------------------------------\n\n"

        ret += "Lame's first and second parameter; Ref.[3]\n"
        ret += "Lambda (GPa)  =  %10.3f  \n" % self.lambda_lame_coefficient
        ret += "Mu (GPa)  =  %10.3f  \n" % self.mu_lame_coefficient

        ret += "\n------------------------------------------------------------------\n"
        ret += "Bonding information\n"
        ret += "------------------------------------------------------------------\n \n"

        ret += "Kleinman’s parameter; Ref.[4,5]\n"
        ret += "NOTE: K = 0 (1) bending (stretching) would dominate\n"
        ret += "K =  %10.5f  \n" % self.kleinman_parameter
        ret += "\n"
        ret += "Cauchy's Pressure calculated from the relation : CP = C_12 - C_44\n"
        ret += "     CP > 0 (+ve) indicates that ionic bonding dominates\n"
        ret += "     CP < 0 (-ve) indicates that covalent bonding dominates\n"

        ret += "CP (GPa) =  %10.3f  \n" % self.cauchy_pressure
        ret += "Bonding is mainly " + self.bonding_type + "\n"

        ret += "\n------------------------------------------------------------------\n"
        ret += "Elastic Anisotropy\n"
        ret += "------------------------------------------------------------------\n\n"

        ret += (
            "Zener's anisotropy (true for cubic crystals only); Az = %6.3f; Ref.[6]\n"
            % self.A_z
        )
        ret += (
            "Chung-Buessem's anisotropy (true for cubic crystals only); Acb = %6.3f; Ref.[7]\n"
            % self.A_cb
        )
        ret += "Universal anisotropy index; Au = %6.3f; Ref.[8]\n" % self.A_u
        ret += "Log-Euclidean's anisotropy; AL = %6.3f; Ref.[9]\n" % self.A_l

        ret += "\n------------------------------------------------------------------\n"
        ret += "Elastic Wave Velocities and Debye Temperature\n"
        ret += "------------------------------------------------------------------\n\n"

        ret += (
            "Longitudinal wave velocity (vl) : %10.3f m/s; Ref.[10]\n"
            % self.velocity_logitudinal
        )
        ret += (
            "Transverse wave velocity (vt) : %10.3f m/s; Ref.[10]\n"
            % self.velocity_transverse
        )
        ret += (
            "Average wave velocity (vm) : %10.3f m/s; Ref.[10]\n"
            % self.velocity_average
        )
        ret += "Debye temperature  : %10.3f K; Ref.[10]\n" % self.debye_temperature
        ret += "\n"
        ret += "WARNING: The  Debye model for the atomic displacement is based on a monoatomic crystal approximation.\n"
        ret += "Here we consider an averaged mass, in case your crystal has several species.\n"

        ret += "\n------------------------------------------------------------------\n"
        ret += "Melting Temperature\n"
        ret += "------------------------------------------------------------------\n\n"

        ret += "Melting temperature calculated from the empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K); Ref.[11]"
        ret += "Tm =  %10.3f K (plus-minus 555 K) \n" % self.melting_temperature
        ret += "\n\n"
        ret += "WARNING: This is a crude approximation and its validity needs to be checked! \n"

        ret += "\n------------------------------------------------------------------\n"
        ret += "Hardness Analysis\n"
        ret += "------------------------------------------------------------------\n\n"

        H1a, H1b, H2, H3, H4, H5 = self.hardness
        ret += "Hardness (H1a) = {:.2f} GPa;  Ref.[12]\n".format(H1a)
        ret += "Hardness (H1b) = {:.2f} GPa;  Ref.[12]\n".format(H1b)
        ret += "Hardness (H2)  = {:.2f} GPa;  Ref.[13]\n".format(H2)
        ret += "Hardness (H3)  = {:.2f} GPa;  Ref.[14]\n".format(H3)
        ret += "Hardness (H4)  = {:.2f} GPa;  Ref.[15]\n".format(H4)
        ret += "Hardness (H5)  = {:.2f} GPa;  Ref.[16]\n".format(H5)
        ret += "\n"
        ret += """Hardness recommendation model:
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

        ret += "\n------------------------------------------------------------------\n"
        ret += "References\n"
        ret += "------------------------------------------------------------------\n\n"
        ret += "[1] Necessary and Sufficient Elastic Stability Conditions in Various Crystal Systems. Félix Mouhat and François-Xavier Coudert. Phys. Rev. B (2014)\n"
        ret += "[2] Crystal Structures and Elastic Properties of Superhard IrN2 and IrN3 from First Principles. Zhi-jian Wu et al. Phys. Rev. B (2007)\n"
        ret += "[3]  The rock physics handbook, Cam-bridge university press. G. Mavko, T. Mukerji, J. Dvorkin. Cambridge University Press. (2020)\n"
        ret += "[4]  Deformation potentials in silicon. i. uniaxial strain. Leonard Kleinman. Phys.Rev. 128. (1962)\n"
        ret += "[5] Electronic Structure and the Properties of Solids: The Physics of the Chemical Bond. Walter A. Harrison.  (2012)\n"
        ret += "[6] Elasticity and Anelasticity of Metals. Clarence M. Zener et al. The Journal of Physical Chemistry. (1949)\n"
        ret += "[7] The Elastic Anisotropy of Crystals. D. H. Chung and W. R. Buessem. Journal of Applied Physics (1967)\n"
        ret += "[8] Universal Elastic Anisotropy Index. Shivakumar I. Ranganathan et al. Phys. Rev. Lett. (2008)\n"
        ret += "[9] Elastic Anisotropy of Crystals. Christopher M. Kube. AIP Advances. (2016)\n"
        ret += "[10] A Simplified Method for Calculating the Debye Temperature from Elastic Constants. Orson L.Anderson. Journal of Physics and Chemistry of Solids. (1963)\n"
        ret += "[11] Elastic constants versus melting temperature in metals. M.E.Fine et al. Scripta Metallurgica. (1984)\n"
        ret += "[12] Correlation between hardness and elastic moduli of the covalent crystals. Jiang, et al. (2011).\n"
        ret += "[13] Computational alchemy: the search for new superhard materials. Teter (1998).\n"
        ret += "[14] Mechanical and electronic properties of B12-based ternary crystals of orthorhombic phase. Jiang et al. (2010).\n"
        ret += "[15] Theoretical investigation on the transition-metal borides with Ta3B4-type structure: A class of hard and refractory materials. Miao et al. (2011).\n"
        ret += "[16] Modeling hardness of polycrystalline materials and bulk metallic glasses. Chen et al. (2011).\n"
        return ret

    @property
    def elastic_stability(self):
        return eigenvals.positive_evals(self.elastic_tensor, verbose=self.verbose)

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


# def positive_evals(cnew, verbose=True):
#     """This method checks the postivity of the eigenvalues
#     of a matrix."""

#     # print("Eigen Values of the matrix:")
#     evals = list(np.linalg.eigvals(cnew))
#     evalsPrint = list(np.around(np.array(evals), 3))
#     if verbose:
#         print("%s" % evalsPrint)
#     check = 0
#     for i in range(len(evals)):
#         if evals[i] > 0.0:
#             pass
#         else:
#             check = 1
#     # if check == 1:
#     #     print(
#     #         "ATTENTION: One or more eigen values are negative indicating elastic instability."
#     #     )
#     # if check == 0:
#     #     print("All eigen values are positive indicating tic stability.")

#     return not (bool(check))


def print_matrix(c):
    print("\n------------------------------------------------------------------")
    print("Elastic Tensor (GPa units)")
    print("------------------------------------------------------------------\n")
    row = c.shape[0]
    col = c.shape[1]
    for i in range(row):
        for j in range(col):
            print("{:>10.3f} ".format(c[i, j]), end=" ")
            if j == (col - 1):
                print(" ")
    print("Note: For VASP users this is the pressure-corrected matrix")

    print("\n------------------------------------------------------------------")
    print("Elastic Tensor Eigen Values (GPa units)")
    print("------------------------------------------------------------------\n")

    return
