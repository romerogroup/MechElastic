#!/usr/bin/env python

import numpy as np
from scipy.optimize import leastsq, least_squares
from scipy import integrate
import matplotlib.pyplot as plt
import sys
from intersect import intersection
import networkx as nx
import matplotlib as mpl
import json
from dicttoxml import dicttoxml
from xml.dom.minidom import parseString

# Setting up plotting class
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams["font.family"] = "Arial"
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

    def __init__(self):
        """Initialize with obtaining fitting parameters
        for provided volume and energy data.
        If eostype=="energy" then it will take in data of
        energy and volume and differentiate for pressure.
        If eostype=="volume" then it will do similarly for
        pressure.
        If eostype=="phase_transition_pressure" it would
        plot the H vs P curves for each phase provided
        by the user as a list for infile.

        """

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

        def P(e1, e2, dv):
            return -1 * (e2 - e1) / (2.00 * dv)

        if self.eostype != "enthalpy":

            v = self.vol_array

            v1 = v - dv
            v2 = v + dv

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

        else:
            self.pressure = []
            for i, icoeffs in enumerate(self.selected_coeffs):

                v = self.vol_array[i]

                v1 = v - dv
                v2 = v + dv

                if self.model[i] == "Vinet":
                    e1 = self.eos_vinet(icoeffs, v1)
                    e2 = self.eos_vinet(icoeffs, v2)
                    self.pressure.append(P(e1, e2, dv))

                elif self.model[i] == "Birch":
                    e1 = self.eos_birch(icoeffs, v1)
                    e2 = self.eos_birch(icoeffs, v2)
                    self.pressure.append(P(e1, e2, dv))

                elif self.model[i] == "Murnaghan":
                    e1 = self.eos_murnaghan(icoeffs, v1)
                    e2 = self.eos_murnaghan(icoeffs, v2)
                    self.pressure.append(P(e1, e2, dv))

                elif self.model[i] == "Birch-Murnaghan":
                    e1 = self.eos_birch_murnaghan(icoeffs, v1)
                    e2 = self.eos_birch_murnaghan(icoeffs, v2)
                    self.pressure.append(P(e1, e2, dv))

            # Converting self.pressure into numpy array and converting to GPa
            self.pressure = 1.60217e02 * np.array(self.pressure, dtype="float64")

    def _get_eeos(self):
        """get_eeos.

        This method uses scipy.integrate to integrate
        pressure to obtain energy.


        """

        v = self.volume
        P_birch_murnaghan = self.eos_birch_murnaghan_pressure(
            self.eos_birch_murnaghan_pressure_fitted.x, v
        )

        self.E_birch_murnaghan = (
            integrate.cumtrapz(v, P_birch_murnaghan, initial=0) * 6.2415e-03
        )  # eV

    def plot_eos(
        self,
        infile=None,
        eostype="energy",
        natoms=1,
        au=False,
        vlim=None,
        model=None,
        raw_data=True,
        export=True,
        savefig=False,
    ):
        """plot_eos.

        This funciton plots the fitted EOS models.
        Initialize with obtaining fitting parameters
        for provided volume and energy data.
        If eostype=="energy" then it will take in data of
        energy and volume and differentiate for pressure.
        If eostype=="volume" then it will do similarly for
        pressure.

        Parameters
        ----------
            infile : list, optional (default ``None``)
                List of input files.

            eostype : str, optional (default ``"energy"``)
                Type of EOS input data; energy or pressure.

            natoms : int, optional (default ``1``)
                Number of atoms

            au : bool, optional (default ``False``)
                Output in atomic units.

            vlim : list, optional (default ``None``)
               Minimum and maximum for volume range.

            model : str, optional (default ``None``)
                Name of model; Vinet, Birch, Murnaghan, Birch-Murnaghan

            raw_data : bool, optional (default ``True``)
                Also include the raw input data in the plot.

            export : bool, optional (default ``True``)
                Export data as .json and .xml.

            savefig : bool, optional (default ``False``)
                Save figures in the pdf format.


        Returns
        -------
            param : None
        """

        self.eostype = eostype
        self.natoms = natoms
        self.au = au
        self.vlim = vlim
        self.model = model
        self.raw_data = raw_data
        self.export = export
        self.savefig = savefig

        if self.eostype == "energy":

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
            E0 = a * V0 ** 2 + b * V0 + c
            B0 = 2 * a * V0
            Bp = 4.0

            self.coeffs0 = [E0, B0, Bp, V0]

        elif self.eostype == "pressure":

            # reading the files with volume and energy
            # columns.

            fi = open(infile, "r")
            data = fi.readlines()
            fi.close()

            self.volume = []
            self.pressure = []

            for i in range(len(data)):
                self.volume.append(float(data[i].split()[0]))
                self.pressure.append(float(data[i].split()[1]))

            self.volume = np.array(self.volume)
            self.pressure = np.array(self.pressure)

            a, b, c = np.polyfit(self.volume, self.pressure, 2)

            V0 = -b / (2 * a)
            P0 = a * V0 ** 2 + b * V0 + c
            K0 = 2 * a * V0
            Kp = 4.0

            self.coeffs0 = [P0, K0, Kp, V0]

        # Interpolating volume array
        if self.vlim:
            self.vol_array = np.linspace(self.vlim[0], self.vlim[1], 1000)
        else:
            self.vol_array = np.linspace(np.min(self.volume), np.max(self.volume), 1000)

        # Call fitting function
        self.fitting()

        if self.eostype == "energy":
            # Call Pressure calculation function
            self._get_peos()
        elif self.eostype == "pressure":
            # Call Energy calculation function
            self._get_eeos()

        fig = plt.figure(figsize=(13, 9))
        axs = fig.add_subplot(111)

        # eV to Ha and A^3 to Bohr^3 conversion if au=True
        if self.au:
            au_converter_e = 0.0367493
            au_converter_v = 6.748376012323316
        else:
            au_converter_e = 1.0
            au_converter_v = 1.0

        if self.eostype == "energy":

            # Energy vs Volume plot

            if self.model is None:

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_vinet(self.eos_vinet_fitted.x, self.vol_array)
                    / self.natoms,
                    "orangered",
                    linewidth=11,
                    label="Vinet",
                )

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_birch(self.eos_birch_fitted.x, self.vol_array)
                    / self.natoms,
                    "lime",
                    linestyle="solid",
                    linewidth=7,
                    label="Birch",
                )

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_murnaghan(self.eos_murnaghan_fitted.x, self.vol_array)
                    / self.natoms,
                    color="magenta",
                    linestyle="dashed",
                    linewidth=5,
                    label="Murnaghan",
                )

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_birch_murnaghan(
                        self.eos_birch_murnaghan_fitted.x, self.vol_array
                    )
                    / self.natoms,
                    "k--",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

                # axs.scatter(
                #     au_converter_v * self.volume,
                #     au_converter_e * self.energy / self.natoms,
                #     s=600,
                #     facecolors="none",
                #     edgecolors="black",
                #     label="Raw Data",
                # )

            elif self.model == "Vinet":
                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_vinet(self.eos_vinet_fitted.x, self.vol_array)
                    / self.natoms,
                    "r",
                    linewidth=2,
                    label="Vinet",
                )

            elif self.model == "Birch":

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_birch(self.eos_birch_fitted.x, self.vol_array)
                    / self.natoms,
                    "r",
                    linewidth=2,
                    label="Birch",
                )

            elif self.model == "Murnaghan":
                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_murnaghan(self.eos_murnaghan_fitted.x, self.vol_array)
                    / self.natoms,
                    color="r",
                    linewidth=2,
                    label="Murnaghan",
                )

            elif self.model == "Birch-Murnaghan":
                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    au_converter_e
                    * self.eos_birch_murnaghan(
                        self.eos_birch_murnaghan_fitted.x, self.vol_array
                    )
                    / self.natoms,
                    "r",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            if self.raw_data:
                axs.scatter(
                    au_converter_v * self.volume / self.natoms,
                    au_converter_e * self.energy / self.natoms,
                    s=600,
                    facecolors="none",
                    edgecolors="black",
                    label="Raw Data",
                )

            if self.au:
                axs.set_xlabel("Volume ($Bohr^3$/atom)")
                axs.set_ylabel("Energy (Ha/atom)")
            else:
                axs.set_xlabel("Volume ($\AA^3$/atom)")
                axs.set_ylabel("Energy (eV/atom)")
            axs.set_title("Energy vs. Volume")
            plt.legend(loc="best")
            if self.savefig:
                plt.savefig("EvsV.pdf")
            plt.show()

            # Pressure vs Volume

            fig = plt.figure(figsize=(13, 9))
            axs2 = fig.add_subplot(111)

            if self.model is None:

                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_vinet,
                    "orangered",
                    linewidth=11,
                    label="Vinet",
                )

                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_birch,
                    "lime",
                    linestyle="solid",
                    linewidth=7,
                    label="Birch",
                )

                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_murnaghan,
                    color="magenta",
                    linestyle="dashed",
                    linewidth=5,
                    label="Murnaghan",
                )

                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_birch_murnaghan,
                    "k--",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            elif self.model == "Vinet":
                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_vinet,
                    "r",
                    linewidth=2,
                    label="Vinet",
                )

            elif self.model == "Birch":
                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_birch,
                    "r",
                    linewidth=2,
                    label="Birch",
                )

            elif self.model == "Murnaghan":
                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_murnaghan,
                    color="r",
                    linewidth=2,
                    label="Murnaghan",
                )

            elif self.model == "Birch-Murnaghan":
                axs2.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.P_birch_murnaghan,
                    "r",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            if self.au:
                axs2.set_xlabel("Volume ($Bohr^3$/atom)")
            else:
                axs2.set_xlabel("Volume ($\AA^3$/atom)")
            axs2.set_ylabel("Pressure (GPa)")
            axs2.set_title("Pressure vs. Volume")
            plt.legend(loc="best")
            if self.savefig:
                plt.savefig("PvsV.pdf")
            plt.show()

        elif self.eostype == "pressure":
            # Pressure vs Volume plot

            ## TODO: Add other models

            if self.model is None:

                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.eos_birch_murnaghan_pressure(
                        self.eos_birch_murnaghan_pressure_fitted.x, self.vol_array
                    ),
                    "r-",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

                # axs.scatter(
                #     au_converter_v * self.volume,
                #     self.pressure / self.natoms,
                #     s=600,
                #     facecolors="none",
                #     edgecolors="black",
                #     label="Raw Data",
                # )

            elif self.model == "Birch-Murnaghan":
                axs.plot(
                    au_converter_v * self.vol_array / self.natoms,
                    self.eos_birch_murnaghan_pressure(
                        self.eos_birch_murnaghan_pressure_fitted.x, self.vol_array
                    ),
                    "r-",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            if self.raw_data:
                axs.scatter(
                    au_converter_v * self.volume / self.natoms,
                    self.pressure / self.natoms,
                    s=600,
                    facecolors="none",
                    edgecolors="black",
                    label="Raw Data",
                )

            if self.au:
                axs.set_xlabel("Volume ($Bohr^3$/atom)")
            else:
                axs.set_xlabel("Volume ($\AA^3$/atom)")
            axs.set_ylabel("Pressure (GPa)")
            axs.set_title("Pressure vs. Volume")
            plt.legend(loc="best")
            if self.savefig:
                plt.savefig("PvsV.pdf")
            plt.show()

            # Energy vs Volume

            fig = plt.figure(figsize=(13, 9))
            axs2 = fig.add_subplot(111)

            if self.model is None:

                axs2.plot(
                    au_converter_v * self.volume / self.natoms,
                    au_converter_e * self.E_birch_murnaghan / self.natoms,
                    "r-",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            elif self.model == "Birch-Murnaghan":
                axs2.plot(
                    au_converter_v * self.volume / self.natoms,
                    au_converter_e * self.E_birch_murnaghan / self.natoms,
                    "r-",
                    linewidth=2,
                    label="Birch-Murnaghan",
                )

            if self.au:
                axs2.set_ylabel("Energy (Ha/atom)")
                axs2.set_xlabel("Volume ($Bohr^3$/atom)")
            else:
                axs2.set_ylabel("Energy (eV/atom)")
                axs2.set_xlabel("Volume ($\AA^3$/atom)")
            axs2.set_title("Energy vs. Volume")
            plt.legend(loc="best")
            if self.savefig:
                plt.savefig("EvsV.pdf")
            plt.show()

        # Saving as .json and .xml
        if self.export:
            self.to_json()
            self.to_xml()

    def plot_enthalpy_curves(
        self,
        infiles=None,
        natoms=1,
        au=False,
        vlim=None,
        vlim_list=None,
        deltaH_index=None,
        labels=None,
        export=True,
        savefig=False,
        intersect_tol=1E-03
    ):
        """plot_enthalpy_curves.

        Initialize with obtaining fitting parameters
        for provided volume and energy data.
        This function plots Enthalpy vs Pressure curves
        for multiple datasets and calculates the phase transition pressure.
        Additionally it plots enthalpy differences with respect to a selected
        phase.


        Parameters
        ----------
            infiles : list, optional (default ``None``)
                List of input files.

            natoms : int, optional (default ``1``)
                Number of atoms.

            labels : str, optional (default ``None``)
                List of names of phases. Supports Latex.
                If not provided infiles will be used.

            au : bool, optional (default ``False``)
                Output in atomic units?

            vlim : list, optional (default ``None``)
                Minimum and maximum volume.

            vlim_list : list, optional (default ``None``)
                List of minimum and maximum volumes for each file.

            deltaH_index : int, optional (default ``None``)
                Index of file to calculate enthalpy difference with respect to.

            export : bool, optional (default ``True``)
                Export data as .json and .xml.

            savefig : bool, optional (default ``False``)
                Save figures in the pdf format.

            intersect_tol : float, optional (default ``1E-03``)
                Tolerance for finding intersection points between enthalpy curves.

        Returns
        -------
            param : None
        """

        self.eostype = "enthalpy"
        self.infiles = infiles
        self.natoms = natoms
        self.au = au
        self.vlim = vlim
        self.vlim_list = vlim_list
        self.deltaH_index = deltaH_index
        self.labels = labels
        self.export = export
        self.savefig = savefig
        self.intersect_tol = intersect_tol

        nfiles = len(self.infiles)
        self.volume = []
        self.energy = []
        self.vol_array = []

        # read a list of files for different phases.
        for ifile, fname in enumerate(self.infiles):
            volume_temp = []
            energy_temp = []
            with open(fname) as fi:
                rawdata = fi.readlines()
                for i in range(len(rawdata)):
                    volume_temp.append(float(rawdata[i].split()[0]))
                    energy_temp.append(float(rawdata[i].split()[1]))

            self.volume.append(volume_temp)
            self.energy.append(energy_temp)

        # Converting 2D lists to numpy arrays
        self.volume = np.array(self.volume, dtype="float64")
        self.energy = np.array(self.energy, dtype="float64")

        # initial fitting coefficients
        self.coeffs0 = []
        for i in range(nfiles):
            a, b, c = np.polyfit(self.volume[i], self.energy[i], 2)

            V0 = -b / (2 * a)
            E0 = a * V0 ** 2 + b * V0 + c
            B0 = 2 * a * V0
            Bp = 4.0

            self.coeffs0.append([E0, B0, Bp, V0])

        # Converting to a numpy array
        self.coeffs0 = np.array(self.coeffs0, dtype="float64")

        # Interpolating volume arrays and storing in a 2D self.vol_array
        for i in range(nfiles):
            if self.vlim:
                self.vol_array.append(np.linspace(self.vlim[0], self.vlim[1], 2000))
            elif self.vlim_list:
                self.vol_array.append(
                    np.linspace(self.vlim_list[i][0], self.vlim_list[i][1], 2000)
                )
            else:
                self.vol_array.append(
                    np.linspace(np.min(self.volume[i]), np.max(self.volume[i]), 2000)
                )

        # converting self.vol_array to numpy array
        self.vol_array = np.array(self.vol_array, dtype="float64")

        # Call fitting function
        self.fitting()

        # Call Pressure calculation function
        self._get_peos()

        # Calculating enthalpy for each phase
        # Converting PV from GPaA3 -> eV

        self.H = []
        self.E = []


        for i in range(nfiles):

            if self.model[i] == "Vinet":
                self.H.append(
                    self.eos_vinet(self.selected_coeffs[i], self.vol_array[i])
                    + np.multiply(self.pressure[i], self.vol_array[i])
                    * (1 / 160.2176621)
                )
                self.E.append(
                    self.eos_vinet(self.selected_coeffs[i], self.vol_array[i])
                )

            elif self.model[i] == "Birch":
                self.H.append(
                    self.eos_birch(self.selected_coeffs[i], self.vol_array[i])
                    + np.multiply(self.pressure[i], self.vol_array[i])
                    * (1 / 160.2176621)
                )
                self.E.append(
                    self.eos_birch(self.selected_coeffs[i], self.vol_array[i])
                )

            elif self.model[i] == "Murnaghan":
                self.H.append(
                    self.eos_murnaghan(self.selected_coeffs[i], self.vol_array[i])
                    + np.multiply(self.pressure[i], self.vol_array[i])
                    * (1 / 160.2176621)
                )
                self.E.append(
                    self.eos_murnaghan(self.selected_coeffs[i], self.vol_array[i])
                )

            elif self.model[i] == "Birch-Murnaghan":
                self.H.append(
                    self.eos_birch_murnaghan(self.selected_coeffs[i], self.vol_array[i])
                    + np.multiply(self.pressure[i], self.vol_array[i])
                    * (1 / 160.2176621)
                )
                self.E.append(
                    self.eos_birch_murnaghan(self.selected_coeffs[i], self.vol_array[i])
                )


        # Converting self.H and self.E to numpy array
        self.H = np.array(self.H, dtype="float64")
        self.E = np.array(self.E, dtype="float64")

        # Normalize H for number of atoms
        for i, iatom in enumerate(self.natoms):
            self.H[i] = self.H[i] / iatom
        # Normalize E,V for number of atoms
        for i, iatom in enumerate(self.natoms):
            self.E[i] = self.E[i] / iatom
            self.vol_array[i] = self.vol_array[i] / iatom


        # Convert eV to Ha if au=True
        if self.au:
            self.H = 0.0367493 * self.H
            self.E = 0.0367493 * self.E

        # Enthalpy vs Pressure plots
        fig = plt.figure(figsize=(13, 9))
        axs = fig.add_subplot(111)

        if self.labels is None:
            self.labels = self.infiles

        for i, filename in enumerate(self.infiles):
            axs.plot(self.pressure[i], self.H[i], label=self.labels[i])
        axs.set_xlabel("Pressure (GPa)")
        if self.au:
            axs.set_ylabel("Enthalpy (Ha/atom)")
        else:
            axs.set_ylabel("Enthalpy (eV/atom)")
        axs.set_title("Enthalpy vs. Pressure")
        plt.legend(loc="best")
        axs.grid(color="gainsboro", ls="--", lw=0.6)
        if self.savefig:
            plt.savefig("Enthalpy.pdf")
        plt.show()

        # Energy vs. Volume plots
        fig = plt.figure(figsize=(13, 9))
        axs = fig.add_subplot(111)

        if self.labels is None:
            self.labels = self.infiles

        for i, filename in enumerate(self.infiles):
            axs.plot(self.vol_array[i], self.E[i], label=self.labels[i])
        axs.set_xlabel("Volume ($\AA^3/atom$)")
        if self.au:
            axs.set_ylabel("Energy (Ha/atom)")
        else:
            axs.set_ylabel("Energy (eV/atom)")
        axs.set_title("Energy vs. Volume")
        plt.legend(loc="best")
        axs.grid(color="gainsboro", ls="--", lw=0.6)
        if self.savefig:
            plt.savefig("Energy.pdf")
        plt.show()


        # Matrix to store transition information
        self.trans_mat = []

        # Intersection points
        if nfiles > 1:
            # print("\nEquilibrium enthalpy points between phases : ")
            for i in range(nfiles):
                for j in range(nfiles):
                    if i < j:
                        # print("\n%s and %s:" % (self.infiles[i], self.infiles[j]))
                        self.int_x, self.int_y = intersection(
                            self.pressure[i], self.H[i], self.pressure[j], self.H[j]
                        )

                        # Check if there are multiple intersection points
                        if self.int_x.size:
                            if self.int_x.size < 2:
                                if self.au:
                                    # print(
                                    #     "%.2f GPa at %.2f Ha "
                                    #     % (self.int_x, self.int_y)
                                    # )
                                    self.trans_mat.append(
                                        [
                                            self.infiles[i],
                                            self.infiles[j],
                                            self.int_x,
                                            self.int_y,
                                        ]
                                    )
                                else:
                                    # print(
                                    #     "%.2f GPa at %.2f eV "
                                    #     % (self.int_x, self.int_y)
                                    # )
                                    self.trans_mat.append(
                                        [
                                            self.infiles[i],
                                            self.infiles[j],
                                            self.int_x,
                                            self.int_y,
                                        ]
                                    )

                            else:
                                for jj in range(self.int_x.size):
                                    if self.au:
                                        # print(
                                        #     "%.2f GPa at %.2f Ha "
                                        #     % (self.int_x[jj], self.int_y[jj])
                                        # )
                                        self.trans_mat.append(
                                            [
                                                self.infiles[i],
                                                self.infiles[j],
                                                self.int_x,
                                                self.int_y,
                                            ]
                                        )

                                    else:
                                        # print(
                                        #     "%.2f GPa at %.2f eV "
                                        #     % (self.int_x[jj], self.int_y[jj])
                                        # )
                                        self.trans_mat.append(
                                            [
                                                self.infiles[i],
                                                self.infiles[j],
                                                self.int_x,
                                                self.int_y,
                                            ]
                                        )

                        # else:
                        #     print("NOT FOUND FOR THIS VOLUME RANGE!")
                        # # plt.plot(x, y, "*k")
        self.trans_mat = np.array(self.trans_mat, dtype="object")

        # reordered trans_mat to contain minimum enthalpy values for multiple
        # intersections and their corresponding pressure.
        self.trans_mat_reordered = np.zeros((self.trans_mat.shape), dtype="object")
        self.trans_mat_reordered[:, 0:2] = self.trans_mat[:, 0:2]
        self.trans_mat_reordered[:, 3] = [
            min(self.trans_mat[i, 3]).item() for i in range(len(self.trans_mat))
        ]
        indx_arr = [np.argmin(self.trans_mat[i, 3]) for i in range(len(self.trans_mat))]
        self.trans_mat_reordered[:, 2] = [
            self.trans_mat[i, 2][j] for i, j in enumerate(indx_arr)
        ]
        self.trans_mat_reordered = self.trans_mat_reordered[
            np.argsort(self.trans_mat_reordered[:, 3])
        ]
        self.trans_mat_reordered_tmp = self.trans_mat_reordered.copy()

        # Find lowest curve of the two intersecting curves for each
        # row in trans_mat_reordered. c1, c2 are the indexes of the pressure
        # arrays.
        for i in range(len(self.trans_mat_reordered)):
            c1 = np.where(
                np.isclose(
                    self.pressure[self.infiles.index(self.trans_mat_reordered[i][0])],
                    self.trans_mat_reordered[i][2],
                    self.intersect_tol,
                )
            )
            c2 = np.where(
                np.isclose(
                    self.pressure[self.infiles.index(self.trans_mat_reordered[i][1])],
                    self.trans_mat_reordered[i][2],
                    self.intersect_tol,
                )
            )

            # corresponding enthalpy values for the two intersecting curves
            # Adding 500 to see a point a bit down along the interseciton point.
            h1 = self.H[self.infiles.index(self.trans_mat_reordered[i][0])][
                c1[0][-1] + 500
            ]
            h2 = self.H[self.infiles.index(self.trans_mat_reordered[i][1])][
                c1[0][-1] + 500
            ]

            # reorder the phases such that lowest energy phase appears first
            if h1 > h2:
                self.trans_mat_reordered_tmp[i][0] = self.trans_mat_reordered[i][1]
                self.trans_mat_reordered_tmp[i][1] = self.trans_mat_reordered[i][0]
            else:
                self.trans_mat_reordered_tmp[i][0] = self.trans_mat_reordered[i][0]
                self.trans_mat_reordered_tmp[i][1] = self.trans_mat_reordered[i][1]

        # setting to original array
        self.trans_mat_reordered = self.trans_mat_reordered_tmp
        # Adding additional phase for testing
        # tm = np.append(self.trans_mat_reordered, ["Im-3m", "XXX", 8.0, -2.0])
        # self.trans_mat_reordered = tm.reshape(7, 4)

        print("\nPossible transition paths for the provided volume range:\n")
        print(
            "NOTE: Ordered with ascending enthalpy. i.e. topmost transition is the most probable.\n"
        )
        selected_array = []
        selected_array2 = []
        for i in range(len(self.trans_mat_reordered)):
            for j in range(len(self.trans_mat_reordered)):
                if i < j:
                    if self.trans_mat_reordered[i][1] == self.trans_mat_reordered[j][0]:
                        selected_array.append(i)

                        print(
                            "{} -> {} (at {:0.3f} GPa) -> {} (at {:0.3f} GPa)".format(
                                self.trans_mat_reordered[i][0],
                                self.trans_mat_reordered[i][1],
                                self.trans_mat_reordered[i][2],
                                self.trans_mat_reordered[j][1],
                                self.trans_mat_reordered[j][2],
                            ),
                        )

                        for k in range(len(self.trans_mat_reordered)):
                            if j < k:
                                if (
                                    self.trans_mat_reordered[j][1]
                                    == self.trans_mat_reordered[k][0]
                                ):
                                    selected_array2.append(j)

                        if j in selected_array2:
                            print(
                                "-> {} (at {:0.3f} GPa)".format(
                                    self.trans_mat_reordered[k][1],
                                    self.trans_mat_reordered[k][2],
                                ),
                            )

            if i not in selected_array:
                print(
                    "{} -> {} (at {:0.3f} GPa)".format(
                        self.trans_mat_reordered[i][0],
                        self.trans_mat_reordered[i][1],
                        self.trans_mat_reordered[i][2],
                    ),
                )

        # print("\nPossible transition paths for the provided volume range:\n")
        # selected_array = []
        # selected_array2 = []
        # for i in range(len(self.trans_mat_reordered)):
        #     for j in range(len(self.trans_mat_reordered)):

        #         def recursion(i, j):
        #             if i < j:
        #                 if (
        #                     self.trans_mat_reordered[i][1]
        #                     == self.trans_mat_reordered[j][0]
        #                 ):
        #                     selected_array.append([])
        #                     selected_array.append(i)

        #                     print(
        #                         "{} -> {} -> {}".format(
        #                             self.trans_mat_reordered[i][0],
        #                             self.trans_mat_reordered[i][1],
        #                             self.trans_mat_reordered[j][1],
        #                         ),
        #                         end="\n",
        #                     )

        #                     for k in range(len(self.trans_mat_reordered)):
        #                         if j < k:
        #                             if (
        #                                 self.trans_mat_reordered[j][1]
        #                                 == self.trans_mat_reordered[k][0]
        #                             ):
        #                                 selected_array2.append(j)

        #                     if j in selected_array:
        #                         print(
        #                             "-> {}".format(self.trans_mat_reordered[k][1]),
        #                         )

        #     if i not in selected_array:
        #         print(
        #             "{} -> {}".format(
        #                 self.trans_mat_reordered[i][0], self.trans_mat_reordered[i][1]
        #             ),
        #         )



        # Network map for phase transitions
        # Added by Pedram
        cmap = plt.get_cmap("nipy_spectral")
        DG = nx.DiGraph()
        names = np.unique(self.trans_mat_reordered[:, 0:2])

        DG.add_nodes_from(names)
        for trans in self.trans_mat_reordered:
            DG.add_edge(
                trans[0],
                trans[1],
                weight=trans[3],
            )
        pos = nx.layout.spring_layout(DG)

        node_sizes = [3500] * DG.number_of_nodes()
        edge_colors = [
            DG.get_edge_data(edge[0], edge[1])["weight"] for edge in DG.edges
        ]  # self.tras_mat_reordered[:, 3]

        plt.figure(figsize=(13, 9))
        nodes = nx.draw_networkx_nodes(
            DG, pos, node_size=node_sizes, node_color="gainsboro"
        )
        edges = nx.draw_networkx_edges(
            DG,
            pos,
            node_size=node_sizes,
            # arrowstyle="->",
            arrowsize=16,
            edge_color=edge_colors,
            edge_cmap=cmap,
            width=4,
        )

        nx.draw_networkx_labels(DG, pos, font_size=16)

        pc = mpl.collections.PatchCollection(edges, cmap=cmap)
        pc.set_array(edge_colors)
        ax = plt.gca()
        ax.set_axis_off()
        cb = plt.colorbar(pc, ax=ax,orientation="vertical")
        if self.au:
            cb.set_label("Enthalpy (Ha/atom)")
        else:
            cb.set_label("Enthalpy (eV/atom)")
        if self.savefig:
            plt.savefig("Enthalpy-network.pdf")
        plt.show()

        # Enthalpy differences with respect to a selected phase
        if self.deltaH_index:

            # Converting human input index to python 0 index
            self.deltaH_index = int(self.deltaH_index)

            # name of phase
            phase_name = self.infiles[self.deltaH_index]

            # removing selected phase from infiles list
            # self.infiles.remove(phase_name)

            # max and min pressure of selected phase
            self.plim = [
                min(self.pressure[self.deltaH_index]),
                max(self.pressure[self.deltaH_index]),
            ]

            # Pressure matrix with removed row
            # self.pressure = np.delete(self.pressure, self.deltaH_index, axis=0)

            # Enthalpy matrix with removed row
            # self.deltaH = np.delete(
            #     (self.H - self.H[self.deltaH_index]), self.deltaH_index, axis=0
            # )

            # Enthalpy difference
            self.deltaH = self.H - self.H[self.deltaH_index]

            fig = plt.figure(figsize=(13, 9))
            axs2 = fig.add_subplot(111)

            for i, filename in enumerate(self.infiles):
                axs2.plot(self.pressure[i], self.deltaH[i], label=self.labels[i])
            axs2.set_xlabel("Pressure (GPa)")
            if self.au:
                axs2.set_ylabel("$\Delta$H (Ha/atom)")
            else:
                axs2.set_ylabel("$\Delta$H (eV/atom)")
            title_string = "$\Delta$H vs. Pressure wrt " + str(phase_name)
            axs2.set_title(title_string)
            # axs2.axhline(color="black")
            axs2.set_xlim(self.plim)
            plt.legend(loc="best")
            if self.savefig:
                plt.save("Enthalpy-difference.pdf")
            plt.show()

        # Saving as .json and .xml
        if self.export:
            self.to_json()
            self.to_xml()

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

        elif self.eostype == "enthalpy":

            print("\nDetermining best model for each phase...")

            self.selected_coeffs = []
            self.model = []

            for i, ifile in enumerate(self.infiles):

                # Vinet
                def residuals_vinet(coeffs, y, x):
                    return y - self.eos_vinet(coeffs, x)

                self.eos_vinet_fitted = least_squares(
                    residuals_vinet,
                    self.coeffs0[i],
                    args=(self.energy[i], self.volume[i]),
                )

                # Birch
                def residuals_birch(coeffs, y, x):
                    return y - self.eos_birch(coeffs, x)

                self.eos_birch_fitted = least_squares(
                    residuals_birch,
                    self.coeffs0[i],
                    args=(self.energy[i], self.volume[i]),
                )

                # Murnaghan
                def residuals_murnaghan(coeffs, y, x):
                    return y - self.eos_murnaghan(coeffs, x)

                self.eos_murnaghan_fitted = least_squares(
                    residuals_murnaghan,
                    self.coeffs0[i],
                    args=(self.energy[i], self.volume[i]),
                )

                # Birch-Murnaghan
                def residuals_birch_murnaghan(coeffs, y, x):
                    return y - self.eos_birch_murnaghan(coeffs, x)

                self.eos_birch_murnaghan_fitted = least_squares(
                    residuals_birch_murnaghan,
                    self.coeffs0[i],
                    args=(self.energy[i], self.volume[i]),
                )

                eos_name = ["Vinet", "Birch", "Murnaghan", "Birch-Murnaghan"]

                cost_array = [
                    np.mean(self.eos_vinet_fitted.fun ** 2),
                    np.mean(self.eos_birch_fitted.fun ** 2),
                    np.mean(self.eos_murnaghan_fitted.fun ** 2),
                    np.mean(self.eos_birch_murnaghan_fitted.fun ** 2),
                ]

                coeffs_array = [
                    self.eos_vinet_fitted.x,
                    self.eos_birch_fitted.x,
                    self.eos_murnaghan_fitted.x,
                    self.eos_birch_murnaghan_fitted.x,
                ]

                print(
                    "%s :  %s " % (ifile, eos_name[cost_array.index(min(cost_array))])
                )

                self.selected_coeffs.append(
                    coeffs_array[cost_array.index(min(cost_array))]
                )
                self.model.append(eos_name[cost_array.index(min(cost_array))])

    ####################### EOS MODELS for Energy ##############################

    # Murnaghan
    def eos_murnaghan(self, coeffs, vol):
        "Ref: Phys. Rev. B 28, 5480 (1983)"
        E0, B0, Bp, V0 = coeffs
        E = (
            E0
            + B0 / Bp * vol * ((V0 / vol) ** Bp / (Bp - 1.0) + 1.0)
            - V0 * B0 / (Bp - 1.0)
        )
        return E

    # Birch-Murnaghan
    def eos_birch_murnaghan(self, coeffs, vol):
        "Ref: Phys. Rev. B 70, 224107"
        E0, B0, Bp, V0 = coeffs
        eta = (vol / V0) ** (1.0 / 3.0)
        E = E0 + 9.0 * B0 * V0 / 16.0 * (eta ** 2 - 1.0) ** 2 * (
            6.0 + Bp * (eta ** 2 - 1.0) - 4.0 * eta ** 2
        )
        return E

    # Birch
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

    # Vinet
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

    # Outputting data to files

    def to_dict(self):
        if self.eostype == "enthalpy":
            return {
                "H": self.H.tolist(),
                "au": self.au,
                "coeffs0": self.coeffs0.tolist(),
                "deltaH_index": self.deltaH_index,
                "energy": self.energy.tolist(),
                "eos_birch_fitted_cost": self.eos_birch_fitted["cost"],
                "eos_birch_murnaghan_fitted_cost": self.eos_birch_murnaghan_fitted[
                    "cost"
                ],
                "eos_murnaghan_fitted_cost": self.eos_murnaghan_fitted["cost"],
                "eos_vinet_fitted_cost": self.eos_vinet_fitted["cost"],
                "eostype": self.eostype,
                "infiles": self.infiles,
                "labels": self.labels,
                "model": self.model,
                "natoms": self.natoms,
                "pressure": self.pressure.tolist(),
                "selected_coeffs": [i.tolist() for i in self.selected_coeffs],
                "trans_mat_reordered": self.trans_mat_reordered.tolist(),
                "vlim": self.vlim,
                "vlim_list": self.vlim_list,
                "vol_array": self.vol_array.tolist(),
                "volume": self.volume.tolist(),
            }
        elif self.eostype == "energy":
            return {
                "P_birch": self.P_birch.tolist(),
                "P_birch_murnaghan": self.P_birch_murnaghan.tolist(),
                "P_murnaghan": self.P_murnaghan.tolist(),
                "P_vinet": self.P_vinet.tolist(),
                "au": self.au,
                "coeffs0": self.coeffs0,
                "energy": self.energy.tolist(),
                "eos_birch_fitted_cost": self.eos_birch_fitted["cost"],
                "eos_birch_murnaghan_fitted_cost": self.eos_birch_murnaghan_fitted[
                    "cost"
                ],
                "eos_murnaghan_fitted_cost": self.eos_murnaghan_fitted["cost"],
                "eos_vinet_fitted_cost": self.eos_vinet_fitted["cost"],
                "eostype": self.eostype,
                "model": self.model,
                "natoms": self.natoms,
                "raw_data": self.raw_data,
                "vlim": self.vlim,
                "vol_array": self.vol_array.tolist(),
                "volume": self.volume.tolist(),
            }
        else:
            return {
                "E_birch_murnaghan": self.E_birch_murnaghan.tolist(),
                "au": self.au,
                "coeffs0": self.coeffs0,
                # "energy": self.energy.tolist(),
                "eos_birch_murnaghan_fitted_cost": self.eos_birch_murnaghan_pressure_fitted[
                    "cost"
                ],
                "eostype": self.eostype,
                "model": self.model,
                "natoms": self.natoms,
                "raw_data": self.raw_data,
                "vlim": self.vlim,
                "vol_array": self.vol_array.tolist(),
                "volume": self.volume.tolist(),
            }

    def to_json(self, outfile="eos.json"):

        wf = open(outfile, "w")
        json.dump(self.to_dict(), wf, sort_keys=True, indent=4, separators=(",", ": "))
        wf.close()

    def to_xml(self, outfile="eos.xml"):
        wf = open(outfile, "w")
        xml = dicttoxml(self.to_dict())
        dom = parseString(xml)
        wf.write(dom.toprettyxml())
        wf.close()
