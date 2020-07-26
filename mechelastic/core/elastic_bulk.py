#!/usr/bin/env python

"""
This contains the method to calculate elastic constants
for bulk systems.
"""

import numpy as np
from ..comms import printer
from ..tests import ductile
from ..tests import eigenvals
from ..utils.constants import *
from ..utils.elements import ELEMENTS
from ..utils.crystalutils import *


def elastic_const_bulk(
    cnew, snew, crystal=None, cell=None, density=None, natoms=None, totalmass=None,
):
    """Calculate elastic constants for bulk."""

    ## Bulk: Voigt
    KV = (cnew[0][0] + cnew[1][1] + cnew[2][2]) + 2 * (
        cnew[0][1] + cnew[1][2] + cnew[2][0]
    )
    KV = KV / 9.0

    ## Shear: Voigt
    GV = (
        (cnew[0][0] + cnew[1][1] + cnew[2][2])
        - (cnew[0][1] + cnew[1][2] + cnew[2][0])
        + 3 * (cnew[3][3] + cnew[4][4] + cnew[5][5])
    )
    GV = GV / 15.0

    # Young's: Voigt
    EV = (9 * KV * GV) / (3 * KV + GV)

    # Poisson's ratio: Voigt
    Nu_V = (3 * KV - EV) / (6 * KV)

    # P-wave modulus, M: Voigt
    MV = KV + (4 * GV / 3.0)

    # Reuss Method
    mc_new = np.matrix(cnew)
    mci_new = mc_new.I

    for i in range(0, 6):
        for j in range(0, 6):
            snew[i][j] = mci_new[i, j]

    ## bulk: Reuss
    KR = (snew[0][0] + snew[1][1] + snew[2][2]) + 2 * (
        snew[0][1] + snew[1][2] + snew[2][0]
    )
    KR = 1.0 / KR

    ## Shear: Reuss
    GR = (
        4 * (snew[0][0] + snew[1][1] + snew[2][2])
        - 4 * (snew[0][1] + snew[1][2] + snew[2][0])
        + 3 * (snew[3][3] + snew[4][4] + snew[5][5])
    )
    GR = 15.0 / GR

    # Young's: Reuss
    ER = (9 * KR * GR) / (3 * KR + GR)

    # Poisson's ratio: Reuss
    Nu_R = (3 * KR - ER) / (6 * KR)

    # P-wave modulus, M: Reuss
    MR = KR + (4 * GR / 3.0)

    # Voigt-Reuss-Hill Approximation: average of both methods
    # VRH
    Kvrh = (KV + KR) / 2.0
    Gvrh = (GV + GR) / 2.0
    Evrh = (EV + ER) / 2.0
    Nu_vrh = (Nu_V + Nu_R) / 2.0
    Mvrh = (MV + MR) / 2.0
    KG_ratio_V = KV / GV
    KG_ratio_R = KR / GR
    KG_ratio_vrh = Kvrh / Gvrh
    # Elastic Anisotropy
    # Zener anisotropy for cubic crystals only
    Az = 2 * cnew[3][3] / (cnew[0][0] - cnew[0][1])

    # Ranganathan and Ostoja-Starzewski method: Phys. Rev. Lett. 101, 055504 (2008).
    # for any crystalline symmetry: Universal anisotropy index
    AU = (KV / KR) + 5 * (GV / GR) - 6.0

    # ""'Note: AU is a relative measure of anisotropy with respect to a limiting value. For example, AU does not prove that a crystal having AU = 3 has double the anisotropy of another crystal with AU = 1.5. I""'

    # log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016)
    # ""'AL  CV , CR   is based on the distance between the averaged stiffnesses CV and
    # CR , which is more appropriate. Clearly, AL  CV , CR   is zero when the crystallite is isotropic.

    AL = np.sqrt(5) * 2.303 * np.log(1 + (AU / 5))

    print("\n \n                         Voigt     Reuss    Average")
    print("-------------------------------------------------------")
    print("Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (KV, KR, Kvrh))
    print("Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (GV, GR, Gvrh))
    print("Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (EV, ER, Evrh))
    print("Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh))
    print("P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh))
    print(
        "Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) "
        % (KG_ratio_V, KG_ratio_R, KG_ratio_vrh, ductile.ductile_test(KG_ratio_vrh))
    )
    print("-------------------------------------------------------")

    print(" \n \n \t \t Elastic Anisotropy \n ")
    print("Zener anisotropy (true for cubic crystals only) Az = %10.5f" % Az)
    print(
        "Universal anisotropy index (Ranganathan and Ostoja-Starzewski method; PRL 101, 055504 (2008)) Au = %10.5f"
        % AU
    )
    print(
        "Log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016) AL = %10.5f "
        % AL
    )

    # Calculation of elastic wave velocities and Debye temperature using the values obtained from Voigt-Reuss-Hill Approximation
    G = Gvrh * 1.0e9  # converting from GPa to Pascal units (kg/ms^2)
    K = Kvrh * 1.0e9

    # transverse velocity: Navier's equation
    vt = np.sqrt((G / density))

    # longitudinal velocity: Navier's equation
    vl = np.sqrt(((3 * K + 4 * G) / (3.0 * density)))

    # average
    vm = 1.0 / (np.cbrt((1.0 / 3.0) * (2.0 / (vt * vt * vt) + 1.0 / (vl * vl * vl))))

    # Debye temperature calculated using  Orson Anderson's proposal [Ref- J. Phys. Chem. Solids (1963) 24, 909-917]
    q = np.sum(natoms)
    theta = (
        (h_Planck / kB)
        * vm
        * np.cbrt((3 * q * N_avogadro * density) / (4 * (np.pi) * totalmass))
    )

    # Melting temperature estimated using empirical relation from Ref: Johnston I, Keeler G, Rollins R and Spicklemire S 1996
    # Solid State Physics Simulations, The Consortium for Upper-Level Physics Software (New York: Wiley)
    Tm = 607 + 9.3 * Kvrh

    print(
        "\n \t \t  Elastic wave velocities calculated using Navier's equation  (in m/s units) \n"
    )
    print("----------------------------------------------- ")
    print("Longitudinal wave velocity (vl) : %10.5f " % vl)
    print("Transverse wave velocity (vt) : %10.5f " % vt)
    print("Average wave velocity (vm) : %10.5f " % vm)
    print("Debye temperature  (in K) : %10.5f " % theta)
    print(
        "WARNING: Debye model for the atomic displacement is based on a monoatomic crystal, here we consider an average mass if your crystal has several species"
    )
    #    print "Atomic displacement at 150, 300 and 450 K  (in A^2) : %10.5f %10.5f %10.5f" %(u2FromDebye(mass,natoms,theta,150.),u2FromDebye(mass,natoms,theta,300.), u2FromDebye(mass,natoms,theta,450.))
    print(
        "\nMelting temperature calculated from empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K)"
    )
    print("Tm (in K)=  %10.5f (plus-minus 555 K) " % Tm)
    print("----------------------------------------------- ")

    eigenvals.positive_evals(cnew)

    crystalselect(cnew=cnew, cell=cell, crystaltype=crystal)

    return
