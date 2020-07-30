#!/usr/bin/env python

"""
This file contains messages outputed to the screen.
"""

from ..version import version as __version__
from ..version import date as __date__
from prettytable import PrettyTable


def print_mechelastic():
    # created at http://www.network-science.de/ascii/.
    mechelastic_str = """
      /\/\   ___  ___| |__   /__\ | __ _ ___| |_(_) ___
     /    \ / _ \/ __| '_ \ /_\ | |/ _` / __| __| |/ __|
    / /\/\ \  __/ (__| | | //__ | | (_| \__ \ |_| | (__
    \/    \/\___|\___|_| |_\__/ |_|\__,_|___/\__|_|\___| """
    print(("  \n" + mechelastic_str))
    print("\nA Python library to calculate elastic properties of materials.\n")
    print("Version %s created on %s\n" % (__version__, __date__))
    print("Please cite:")
    print(
        "Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018)."
    )
    print("\nDisclaimer:")
    print(
        "Please check the authenticity of your results before publishing. AUTHORS of this script do not guarantee the quality and/or accuracy of results generated using this script.\n"
    )
    print(
        "------------------------------------------------------------------------------------------"
    )

    return


def print_warning_2D():
    warning_str = """    The conversion of elastic constants from GPa to N/m units was done
              by multiplying the Cij matrix elements by the thickness of the simulation cell,
              i.e. by the out-of-plane 'c' lattice parameters. Of course, this assumes that
              vacuum is along the c-axis and this is a good approximation only for atomically
              thin monolayers. One needs to account for the thickness of the 2D system, if
              the system is considerably thick.
              I HOPE YOU UNDERSTAND WHAT YOU ARE DOING HERE. """
    print(("WARNING:  " + warning_str))

    return


# def printMatrix(c):
#     p = PrettyTable()
#     for row in c:
#         p.add_row(row)
#     print(p.get_string(header=False, border=False))


def printMatrix(c):
    row = c.shape[0]
    col = c.shape[1]
    for i in range(row):
        for j in range(col):
            print("{:>10.4f} ".format(c[i, j]), end=" ")
            if j == (col - 1):
                print(" ")

    return
