#!/usr/bin/env python

import numpy as np


class EOS:
    """EOS.

    This class contains methods to calculate the
    Equation of State for a selection of different
    models.

    """

    def __init__(self, volume, energy, order=2):
        """Initialize with obtaining fitting parameters
        for provided volume and energy data."""

        a, b, c = np.polyfit(volume, energy, order)
        V0 = -b / (2 * a)
        E0 = a * V0 ** 2 + b * V0 + c
        B0 = 2 * a * V0
        Bp = 4.0

        self.x0 = [V0, E0, B0, Bp]

    def fitting(self, volume, energy):
        """Fits using the least square method"""

        # Murnaghan
        target = lambda params, y, x: y - eos_murnaghan(params, x)
        eos_murnaghan_fitted = leastsq(target, self.x0, args=(energy, volume))

    def get_peos(self, params, v, dv=1e-3, eostype=1):
        """get_peos.

        This method uses a central difference method to
        differentiate EOS in energy to get pressure.

        central difference : f(v+dv) - f(v-dv))/2dv

        P = dE/dV

        Parameters
        ----------
        params : float
            params
        v : float
            v
        dv : float
            dv
        eostyype : int
            eostype

        """
        v1 = v - dv
        v2 = v + dv

        # Murnaghan
        if eostype == 1:
            p1 = self.eos_murnaghan(params, v1)
            p2 = self.eos_murnaghan(params, v2)

        # Birch-Murnaghan
        elif eostype == 2:
            p1 = self.eos_birchmurnaghan(params, v1)
            p2 = self.eos_birchmurnaghan(params, v2)

        # Birch
        elif eostype == 3:
            p1 = self.eos_birch(params, v1)
            p2 = self.eos_birch(params, v2)

        P = -1 * (p2 - p1) / (2.00 * dv)

        return P

    # Murnaghan equation of state
    def eos_murnaghan(self, params, vol):
        "From Phys. Rev. B 28, 5480 (1983)"
        E0, B0, Bp, V0 = params
        E = (
            E0
            + B0 / Bp * vol * ((V0 / vol) ** Bp / (Bp - 1.0) + 1.0)
            - V0 * B0 / (Bp - 1.0)
        )
        return E

    # Birch-Murnaghan equation of state
    def eos_birch_murnaghan(self, params, vol):
        "From Phys. Rev. B 70, 224107"
        E0, B0, Bp, V0 = params
        eta = (vol / V0) ** (1.0 / 3.0)
        E = E0 + 9.0 * B0 * V0 / 16.0 * (eta ** 2 - 1.0) ** 2 * (
            6.0 + Bp * (eta ** 2 - 1.0) - 4.0 * eta ** 2
        )
        return E

    # Birch equation of state
    def eos_birch(self, params, vol):
        """
        From Intermetallic compounds: Principles and Practice, Vol. I: Princples
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
        """
        E0, B0, Bp, V0 = params
        E = (
            E0
            + 9.0 / 8.0 * B0 * V0 * ((V0 / vol) ** (2.0 / 3.0) - 1.0) ** 2
            + 9.0 / 16.0 * B0 * V0 * (Bp - 4.0) * ((V0 / vol) ** (2.0 / 3.0) - 1.0) ** 3
        )
        return E

    # Vinet equation of state
    def eos_vinet(self, params, vol):
        "From Phys. Rev. B 70, 224107"
        E0, B0, Bp, V0 = params
        eta = (vol / V0) ** (1.0 / 3.0)
        E = E0 + 2.0 * B0 * V0 / (Bp - 1.0) ** 2 * (
            2.0
            - (5.0 + 3.0 * Bp * (eta - 1.0) - 3.0 * eta)
            * numpy.exp(-3.0 * (Bp - 1.0) * (eta - 1.0) / 2.0)
        )
        return E
