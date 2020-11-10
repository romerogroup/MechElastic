#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq, least_squares
import matplotlib.pyplot as plt

# Setting up plotting class
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams["font.family"] = "Georgia"
plt.rc("font", size=22)  # controls default text sizes
plt.rc("axes", titlesize=22)  # fontsize of the axes title
plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
plt.rc("ytick", labelsize=22)  # fontsize of the tick labels


class EOS:

    """EOS.

    This class contains methods to calculate the
    Equation of State for a selection of different
    models.

    """

    def __init__(self, infile=None):
        """Initialize with obtaining fitting parameters
        for provided volume and energy data."""

        # reading the files with volume and energy
        # columns.

        fi = open(infile, "r")
        data = fi.readlines()
        fi.close()

        self.volume = []
        self.energy = []

        for i in range(len(data)):
            self.volume.append(float(data[i].split()[0]))
            self.energy.append(float(data[i].split()[1]))

        self.volume = np.array(self.volume)
        self.energy = np.array(self.energy)

        a, b, c = np.polyfit(self.volume, self.energy, 2)

        V0 = -b / (2 * a)
        self.E0 = a * V0 ** 2 + b * V0 + c
        B0 = 2 * a * V0
        Bp = 4.0

        self.coeffs0 = [self.E0, B0, Bp, V0]

        return

    def _get_peos(self, dv=1e-3):
        """get_peos.

        This method uses a central difference method to
        differentiate EOS in energy to get pressure.

        central difference : (f(v+dv) - f(v-dv))/2dv

        P = dE/dV

        Parameters
        ----------
        v : float
            v
        dv : float
            dv

        """

        v = self.volume

        v1 = v - dv
        v2 = v + dv

        def P(e1, e2, dv):
            return -1 * (e2 - e1) / (2.00 * dv)

        # Murnaghan
        e1 = self.eos_murnaghan(self.eos_murnaghan_fitted.x, v1)
        e2 = self.eos_murnaghan(self.eos_murnaghan_fitted.x, v2)
        self.P_murnaghan = P(e1, e2, dv)

        # Birch-Murnaghan
        e1 = self.eos_birch_murnaghan(self.eos_birch_murnaghan_fitted.x, v1)
        e2 = self.eos_birch_murnaghan(self.eos_birch_murnaghan_fitted.x, v2)
        self.P_birch_murnaghan = P(e1, e2, dv)

        # Birch
        e1 = self.eos_birch(self.eos_birch_fitted.x, v1)
        e2 = self.eos_birch(self.eos_birch_fitted.x, v2)
        self.P_birch = P(e1, e2, dv)

        # Vinet
        e1 = self.eos_vinet(self.eos_vinet_fitted.x, v1)
        e2 = self.eos_vinet(self.eos_vinet_fitted.x, v2)
        self.P_vinet = P(e1, e2, dv)

    def plot_eos(self):
        """plot_eos.

        This funciton plots the fitted EOS models.

        """

        # Call fitting function
        self.fitting()

        # Call Pressure calculation funciton
        self._get_peos()

        vol_array = np.linspace(min(self.volume), max(self.volume), 100)

        fig = plt.figure(figsize=(13, 9))
        axs = fig.add_subplot(111)

        # Energy vs Volume plot

        axs.plot(
            vol_array,
            self.eos_vinet(self.eos_vinet_fitted.x, vol_array),
            "orangered",
            linewidth=11,
            label="Vinet",
        )

        axs.plot(
            vol_array,
            self.eos_birch(self.eos_birch_fitted.x, vol_array),
            "lime",
            linestyle="solid",
            linewidth=7,
            label="Birch",
        )

        axs.plot(
            vol_array,
            self.eos_murnaghan(self.eos_murnaghan_fitted.x, vol_array),
            color="magenta",
            linestyle="dashed",
            linewidth=5,
            label="Murnaghan",
        )

        axs.plot(
            vol_array,
            self.eos_birch_murnaghan(self.eos_birch_murnaghan_fitted.x, vol_array),
            "k--",
            linewidth=2,
            label="Birch-Murnaghan",
        )

        axs.scatter(
            self.volume,
            self.energy,
            s=600,
            facecolors="none",
            edgecolors="black",
            label="Raw Data",
        )

        axs.set_xlabel("Volume ($\AA^3$)")
        axs.set_ylabel("Energy (eV)")
        axs.set_title("Energy vs. Volume")
        plt.legend(loc="best")
        plt.show()

        # Pressure vs Volume

        fig = plt.figure(figsize=(13, 9))
        axs2 = fig.add_subplot(111)

        axs2.plot(
            self.volume, self.P_murnaghan, "orangered", linewidth=11, label="Vinet",
        )

        axs2.plot(
            self.volume,
            self.P_murnaghan,
            "lime",
            linestyle="solid",
            linewidth=7,
            label="Birch",
        )

        axs2.plot(
            self.volume,
            self.P_murnaghan,
            color="magenta",
            linestyle="dashed",
            linewidth=5,
            label="Murnaghan",
        )

        axs2.plot(
            self.volume, self.P_murnaghan, "k--", linewidth=2, label="Birch-Murnaghan",
        )

        axs2.set_xlabel("Volume ($\AA^3$)")
        axs2.set_ylabel("Pressure (eV/$\AA^3$)")
        axs2.set_title("Pressure vs. Volume")
        plt.legend(loc="best")
        plt.show()

    ############################# FITTING #############################

    def fitting(self):
        """fitting.

        Fits coefficients using the least square method.
        returns an array of fitted coefficients [E0,B0,Bp,V0].

        """

        print("Fitting coefficients : [E0,B0,Bp,V0]\n")

        # Vinet
        def residuals_vinet(coeffs, y, x):
            return y - self.eos_vinet(coeffs, x)

        self.eos_vinet_fitted = least_squares(
            residuals_vinet, self.coeffs0, args=(self.energy, self.volume)
        )
        print("Vinet : %s" % self.eos_vinet_fitted.x)
        print("Vinet Cost (MSE) : %s\n" % self.eos_vinet_fitted.cost)

        # Birch
        def residuals_birch(coeffs, y, x):
            return y - self.eos_birch(coeffs, x)

        self.eos_birch_fitted = least_squares(
            residuals_birch, self.coeffs0, args=(self.energy, self.volume)
        )
        print("Birch : %s" % self.eos_birch_fitted.x)
        print("Birch Cost (MSE) : %s\n" % self.eos_birch_fitted.cost)

        # Murnaghan
        def residuals_murnaghan(coeffs, y, x):
            return y - self.eos_murnaghan(coeffs, x)

        self.eos_murnaghan_fitted = least_squares(
            residuals_murnaghan, self.coeffs0, args=(self.energy, self.volume)
        )

        print("Murnaghan : %s" % self.eos_murnaghan_fitted.x)
        print("Murnaghan Cost (MSE) : %s\n" % self.eos_murnaghan_fitted.cost)

        # Birch-Murnaghan
        def residuals_birch_murnaghan(coeffs, y, x):
            return y - self.eos_birch_murnaghan(coeffs, x)

        self.eos_birch_murnaghan_fitted = least_squares(
            residuals_birch_murnaghan, self.coeffs0, args=(self.energy, self.volume)
        )
        print("Birch-Murnaghan : %s" % self.eos_birch_murnaghan_fitted.x)
        print(
            "Birch-Murnaghan Cost (MSE) : %s\n" % self.eos_birch_murnaghan_fitted.cost
        )

        cost_array = [
            self.eos_vinet_fitted.cost,
            self.eos_birch_fitted.cost,
            self.eos_murnaghan_fitted.cost,
            self.eos_birch_murnaghan_fitted.cost,
        ]
        eos_name = ["Vinet", "Birch", "Murnaghan", "Birch-Murnaghan"]

        print(
            "Based on MSE %s is the best EOS model for this dataset."
            % (eos_name[cost_array.index(min(cost_array))])
        )

    ####################### EOS MODELS ##############################

    # Murnaghan equation of state
    def eos_murnaghan(self, coeffs, vol):
        "Ref: Phys. Rev. B 28, 5480 (1983)"
        E0, B0, Bp, V0 = coeffs
        E = (
            E0
            + B0 / Bp * vol * ((V0 / vol) ** Bp / (Bp - 1.0) + 1.0)
            - V0 * B0 / (Bp - 1.0)
        )
        return E

    # Birch-Murnaghan equation of state
    def eos_birch_murnaghan(self, coeffs, vol):
        "Ref: Phys. Rev. B 70, 224107"
        E0, B0, Bp, V0 = coeffs
        eta = (vol / V0) ** (1.0 / 3.0)
        E = E0 + 9.0 * B0 * V0 / 16.0 * (eta ** 2 - 1.0) ** 2 * (
            6.0 + Bp * (eta ** 2 - 1.0) - 4.0 * eta ** 2
        )
        return E

    # Birch equation of state
    def eos_birch(self, coeffs, vol):
        """
        Ref: Intermetallic compounds: Principles and Practice, Vol. I: Princples
        Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
        """
        E0, B0, Bp, V0 = coeffs
        E = (
            E0
            + 9.0 / 8.0 * B0 * V0 * ((V0 / vol) ** (2.0 / 3.0) - 1.0) ** 2
            + 9.0 / 16.0 * B0 * V0 * (Bp - 4.0) * ((V0 / vol) ** (2.0 / 3.0) - 1.0) ** 3
        )
        return E

    # Vinet equation of state
    def eos_vinet(self, coeffs, vol):
        "Ref: Phys. Rev. B 70, 224107"
        E0, B0, Bp, V0 = coeffs
        eta = (vol / V0) ** (1.0 / 3.0)
        E = E0 + 2.0 * B0 * V0 / (Bp - 1.0) ** 2 * (
            2.0
            - (5.0 + 3.0 * Bp * (eta - 1.0) - 3.0 * eta)
            * np.exp(-3.0 * (Bp - 1.0) * (eta - 1.0) / 2.0)
        )
        return E


if __name__ == "__main__":
    eos = EOS(infile="EvsV.dat")
    eos.plot_eos()
