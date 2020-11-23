#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq, least_squares
from scipy import integrate
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

    def __init__(self, infile=None, eostype="energy"):
        """Initialize with obtaining fitting parameters
        for provided volume and energy data.
        If eostype=="energy" then it will take in data of
        energy and volume and differentiate for pressure.
        If eostype=="volume" then it will do similarly for
        pressure.

        """

        # reading the files with volume and energy
        # columns.

        fi = open(infile, "r")
        data = fi.readlines()
        fi.close()

        self.eostype = eostype

        if self.eostype == "energy":

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

        elif self.eostype == "pressure":
            self.volume = []
            self.pressure = []

            for i in range(len(data)):
                self.volume.append(float(data[i].split()[0]))
                self.pressure.append(float(data[i].split()[1]))

            self.volume = np.array(self.volume)
            self.pressure = np.array(self.pressure)

            a, b, c = np.polyfit(self.volume, self.pressure, 2)

            V0 = -b / (2 * a)
            self.P0 = a * V0 ** 2 + b * V0 + c
            K0 = 2 * a * V0
            Kp = 4.0

            self.coeffs0 = [self.P0, K0, Kp, V0]

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
        self.P_murnaghan = P(e1, e2, dv) * 1.60217e02  # GPa

        # Birch-Murnaghan
        e1 = self.eos_birch_murnaghan(self.eos_birch_murnaghan_fitted.x, v1)
        e2 = self.eos_birch_murnaghan(self.eos_birch_murnaghan_fitted.x, v2)
        self.P_birch_murnaghan = P(e1, e2, dv) * 1.60217e02  # GPa

        # Birch
        e1 = self.eos_birch(self.eos_birch_fitted.x, v1)
        e2 = self.eos_birch(self.eos_birch_fitted.x, v2)
        self.P_birch = P(e1, e2, dv) * 1.60217e02  # GPa

        # Vinet
        e1 = self.eos_vinet(self.eos_vinet_fitted.x, v1)
        e2 = self.eos_vinet(self.eos_vinet_fitted.x, v2)
        self.P_vinet = P(e1, e2, dv) * 1.60217e02  # GPa

    def _get_eeos(self):
        """get_eeos.

        This method uses scipy.integrate to integrate
        pressure to obtain energy.


        """

        v = self.volume
        P_birch_murnaghan = self.eos_birch_murnaghan_pressure(
            self.eos_birch_murnaghan_pressure_fitted.x, v
        )

        # def P(v):
        #     if np.abs(v) < 1e-10:
        #         res = v
        #     else:
        #         res = self.eos_birch_murnaghan_pressure(
        #             self.eos_birch_murnaghan_pressure_fitted.x, v
        #         )
        #     return res

        # def E(v):
        #     res = np.zeros_like(v)
        #     for i, val in enumerate(v):
        #         print(val)
        #         y, err = quad(P, 0, val)
        #         res[i] = y

        #     print(res * 6.241509 * 1e-03)
        #     return res

        self.E_birch_murnaghan = (
            integrate.cumtrapz(v, P_birch_murnaghan, initial=0) * 6.2415e-03
        )  # eV

    def plot_eos(self):
        """plot_eos.

        This funciton plots the fitted EOS models.

        """

        # Call fitting function
        self.fitting()

        if self.eostype == "energy":
            # Call Pressure calculation function
            self._get_peos()
        elif self.eostype == "pressure":
            # Call Energy calculation function
            self._get_eeos()

        vol_array = np.linspace(min(self.volume), max(self.volume), 100)

        fig = plt.figure(figsize=(13, 9))
        axs = fig.add_subplot(111)

        if self.eostype == "energy":

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
                self.volume, self.P_vinet, "orangered", linewidth=11, label="Vinet",
            )

            axs2.plot(
                self.volume,
                self.P_birch,
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
                self.volume,
                self.P_birch_murnaghan,
                "k--",
                linewidth=2,
                label="Birch-Murnaghan",
            )

            axs2.set_xlabel("Volume ($\AA^3$)")
            axs2.set_ylabel("Pressure (GPa)")
            axs2.set_title("Pressure vs. Volume")
            plt.legend(loc="best")
            plt.show()

        elif self.eostype == "pressure":
            # Pressure vs Volume plot

            axs.plot(
                vol_array,
                self.eos_birch_murnaghan_pressure(
                    self.eos_birch_murnaghan_pressure_fitted.x, vol_array
                ),
                "r-",
                linewidth=2,
                label="Birch-Murnaghan",
            )

            axs.scatter(
                self.volume,
                self.pressure,
                s=600,
                facecolors="none",
                edgecolors="black",
                label="Raw Data",
            )

            axs.set_xlabel("Volume ($\AA^3$)")
            axs.set_ylabel("Pressure (GPa)")
            axs.set_title("Pressure vs. Volume")
            plt.legend(loc="best")
            plt.show()

            # Energy vs Volume

            fig = plt.figure(figsize=(13, 9))
            axs2 = fig.add_subplot(111)

            axs2.plot(
                self.volume,
                self.E_birch_murnaghan,
                "r-",
                linewidth=2,
                label="Birch-Murnaghan",
            )

            axs2.set_xlabel("Volume ($\AA^3$)")
            axs2.set_ylabel("Energy (eV)")
            axs2.set_title("Energy vs. Volume")
            plt.legend(loc="best")
            plt.show()

    ############################# FITTING #############################

    def fitting(self):
        """fitting.

        Fits coefficients using the least square method.
        returns an array of fitted coefficients [E0,B0,Bp,V0] or
        [P0,K0,Kp,V0].

        """

        if self.eostype == "energy":

            print(
                "Fitting coefficients [E0 (eV), B0 (GPa), Bp (GPa), V0 (Anstroms^3)]\nalong with Mean Squared Error (MSE)\n"
            )
            np.set_printoptions(formatter={"float": "{: 0.3f}".format})

            # Vinet
            def residuals_vinet(coeffs, y, x):
                return y - self.eos_vinet(coeffs, x)

            self.eos_vinet_fitted = least_squares(
                residuals_vinet, self.coeffs0, args=(self.energy, self.volume)
            )
            print("Vinet : %s" % self.eos_vinet_fitted.x)
            print(
                "Vinet (MSE) : {:.3e}\n".format(np.mean(self.eos_vinet_fitted.fun ** 2))
            )

            # Birch
            def residuals_birch(coeffs, y, x):
                return y - self.eos_birch(coeffs, x)

            self.eos_birch_fitted = least_squares(
                residuals_birch, self.coeffs0, args=(self.energy, self.volume)
            )
            print("Birch : %s" % self.eos_birch_fitted.x)
            print(
                "Birch (MSE) : {:.3e}\n".format(np.mean(self.eos_birch_fitted.fun ** 2))
            )

            # Murnaghan
            def residuals_murnaghan(coeffs, y, x):
                return y - self.eos_murnaghan(coeffs, x)

            self.eos_murnaghan_fitted = least_squares(
                residuals_murnaghan, self.coeffs0, args=(self.energy, self.volume)
            )

            print("Murnaghan : %s" % self.eos_murnaghan_fitted.x)
            print(
                "Murnaghan (MSE) : {:.3e}\n".format(
                    np.mean(self.eos_murnaghan_fitted.fun ** 2)
                )
            )

            # Birch-Murnaghan
            def residuals_birch_murnaghan(coeffs, y, x):
                return y - self.eos_birch_murnaghan(coeffs, x)

            self.eos_birch_murnaghan_fitted = least_squares(
                residuals_birch_murnaghan, self.coeffs0, args=(self.energy, self.volume)
            )
            print("Birch-Murnaghan : %s" % self.eos_birch_murnaghan_fitted.x)
            print(
                "Birch-Murnaghan (MSE) : {:.3e}\n".format(
                    np.mean(self.eos_birch_murnaghan_fitted.fun ** 2)
                )
            )

            cost_array = [
                np.mean(self.eos_vinet_fitted.fun ** 2),
                np.mean(self.eos_birch_fitted.fun ** 2),
                np.mean(self.eos_murnaghan_fitted.fun ** 2),
                np.mean(self.eos_birch_murnaghan_fitted.fun ** 2),
            ]
            eos_name = ["Vinet", "Birch", "Murnaghan", "Birch-Murnaghan"]

            print(
                "Based on MSE %s is the best EOS model for this dataset."
                % (eos_name[cost_array.index(min(cost_array))])
            )

        elif self.eostype == "pressure":
            print(
                "Fitting coefficients [P0 (GPa), K0 (GPa), Kp (GPa), V0 (Anstroms^3)]\nalong with Mean Squared Error (MSE)\n"
            )
            np.set_printoptions(formatter={"float": "{: 0.3f}".format})

            # Birch-Murnaghan Pressure
            def residuals_birch_murnaghan_pressure(coeffs, y, x):
                return y - self.eos_birch_murnaghan_pressure(coeffs, x)

            self.eos_birch_murnaghan_pressure_fitted = least_squares(
                residuals_birch_murnaghan_pressure,
                self.coeffs0,
                args=(self.pressure, self.volume),
            )
            print("Birch-Murnaghan : %s" % self.eos_birch_murnaghan_pressure_fitted.x)
            print(
                "Birch-Murnaghan (MSE) : {:.3e}\n".format(
                    np.mean(self.eos_birch_murnaghan_pressure_fitted.fun ** 2)
                )
            )

            cost_array = [
                np.mean(self.eos_birch_murnaghan_pressure_fitted.fun ** 2),
            ]
            eos_name = ["Birch-Murnaghan"]

            print(
                "Based on MSE %s is the best EOS model for this dataset."
                % (eos_name[cost_array.index(min(cost_array))])
            )

    ####################### EOS MODELS for Energy ##############################

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
        Ref: Michael J. Mehl; Barry M. Klein; Dimitris A. Papaconstantopoulos. First-Principles Calculation of Elastic Properties. In Intermetallic Compounds; John Wiley & Sons Ltd, 1994; Vol. 1.

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

    ####################### EOS MODELS for Pressure ##############################

    # Birch-Murnaghan Pressure EOS
    def eos_birch_murnaghan_pressure(self, coeffs, vol):
        """
        Ref: doi:10.3390/min9120745
        """
        P0, K0, Kp, V0 = coeffs
        eta = (V0 / vol) ** (1.0 / 3.0)

        P = (
            (3 / 2)
            * K0
            * (eta ** 7 - eta ** 5)
            * (1 + 0.75 * (Kp - 4) * (eta ** 2 - 1))
        )
        return P


if __name__ == "__main__":
    eos = EOS(infile="EvsV.dat")
    eos.plot_eos()
