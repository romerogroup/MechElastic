import numpy as np
import re
from ..utils.constants import *
from ..utils.elements import *

from ..comms import printer
from ..tests import symmetry


class AbinitParser:
    """
    This class contains methods to parse the Abinit output file.
    """

    def __init__(self, infile="abinit.out"):

        self.infile = infile

        self.cnew = np.zeros((6, 6))
        self.snew = np.zeros((6, 6))

        self._readOUTPUT()

        return

    def _readOUTPUT(self):

        # s is the complaiance matrix and c is the elastic constant matrix
        s = np.zeros((6, 6))
        c = np.zeros((6, 6))

        # parsing OUTCAR

        f = open(self.infile, "r")
        lines = f.readlines()
        f.close()

        nions = -1
        mass = []
        foundcell = 0
        icell = -1
        iatom = -1
        lattice = []
        positions = []
        atomic_numbers = []
        symbols = []
        for i in lines:
            word = i.split()
            if "TOTAL ELASTIC MODULI (kBar)" in i:
                ii = lines.index(i)
            if "POMASS" in i and "ZVAL" in i:
                m = i.split()[2]
                mass.append(float(m[:-1]))
            if "volume of cell" in i:
                vol = float(i.split()[4])
            if "NIONS =" in i:
                nions = int(i.split()[11])
            if "ions per type =" in i:
                self.natoms = list(map(int, i.split()[4:]))
            if "LATTYP" in i:
                lattyp = i
            if icell > -1 and not foundcell:
                lattice.append(list(map(float, i.split()))[:3])
                icell = icell + 1
                if icell == 3:
                    foundcell = 1
            if not foundcell and "direct lattice vectors" in i:
                icell = icell + 1
            if nions > 0 and iatom > -1:
                positions.append(list(map(float, word)))
                iatom = iatom + 1
                if iatom == nions:
                    iatom = -1
            if nions > 0 and "position of ions in fractional coordinates" in i:
                iatom = iatom + 1
            if "ions per type = " in i:
                iontype = list(map(int, word[4:]))
            if "VRHFIN" in i:
                # symbols.append(i.split()[1][1:-1])
                # Modified by Uthpala to detect symbols in OUTCAR when there is a space between
                # the element symbol and ':'. E.g.- La : instead of La:
                symbols.append(re.findall(r"=([a-zA-Z\s]*):", i)[0].split()[0])
            if "length of vectors" in i:
                n = lines.index(i)
                l = lines[n + 1]
                self.lattconst = l.split()

        # cell parameters
        A = float(self.lattconst[0])
        B = float(self.lattconst[1])
        C = float(self.lattconst[2])

        print(
            "Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f"
            % (A, B, C)
        )

        atomic_numbers = np.zeros(nions, dtype=np.int32)
        k = 0
        for i in range(len(iontype)):
            for j in range(iontype[i]):
                atomic_numbers[k] = ELEMENTS[symbols[i]]
                k = k + 1

        print("Atomic numbers")
        print(atomic_numbers)
        lattice = np.array(lattice)
        positions = np.array(positions)
        self.cell = (lattice, positions, atomic_numbers)

        mass = np.array(mass)
        print("Mass of atoms (in g/mol units): ")
        print((np.array(mass)))
        print("Number of atoms: %d" % np.sum(self.natoms))

        self.totalmass = 0.0
        for i in range(len(self.natoms)):
            self.totalmass = self.totalmass + self.natoms[i] * mass[i]
        print("Total mass (in g/mol): %10.4f " % self.totalmass)
        print("Volume of the cell (in Ang^3 units): %10.4f " % vol)

        # converting the units
        vol = vol * 1.0e-30  # from Angstrom to meters
        self.totalmass = self.totalmass * 1.0e-3  # from gram to kg
        self.density = self.totalmass / (vol * N_avogadro)

        print("\nDensity (in kg/m^3 units ): %10.5f" % self.density)

        for i in range(0, 6):
            l = lines[ii + 3 + i]
            word = l.split()
            s[i][:] = word[1:7]
            c[i] = list(map(float, s[i][:]))

        # print c
        mc = np.matrix(c)
        mci = mc.I

        for i in range(0, 6):
            for j in range(0, 6):
                s[i][j] = mci[i, j]

        print("Printing Cij matrix as read from OUTCAR\n")
        printer.printMatrix(c)

        # Redefining the Cij matrix into the correct Voigt notation since VASP's OUTCAR has a different order
        # In VASP: Columns and rows are listed as: 1, 2, 3, 6, 4, 5
        # In this format OUTCAR's C44 values would be actually C66, C55 would be C44, and C66 would be C55.
        # OUTCAR follows the below order:
        # [C11 C12 C13 C16 C14 C15]
        # [C21 C22 C23 C26 C24 C25]
        # [C31 C32 C33 C36 C34 C35]
        # [C61 C62 C63 C66 C64 C65]
        # [C41 C42 C43 C46 C44 C45]
        # [C51 C52 C53 C56 C54 C55]

        self.cnew = np.copy(c)

        # for i in range(0,3):
        #   for j in range(0,6):
        #     cnew[i][j]=c[i][j]

        for j in range(0, 6):
            self.cnew[3][j] = c[4][j]
            self.cnew[4][j] = c[5][j]
            self.cnew[5][j] = c[3][j]

        ctemp = np.zeros((6, 6))
        ctemp = np.copy(self.cnew)

        for i in range(0, 6):
            self.cnew[i][3] = self.cnew[i][4]
            self.cnew[i][4] = self.cnew[i][5]
            self.cnew[i][5] = ctemp[i][3]

        # Change the units of Cij from kBar to GPa
        for i in range(0, 6):
            for j in range(0, 6):
                self.cnew[i][j] = self.cnew[i][j] / 10.0

        print(
            "\n \n printing CNEW: Modified matrix in correct order (in GPa units)... \n For example- to generate input for the ELATE code [https://github.com/fxcoudert/elate] \n"
        )

        np.set_printoptions(precision=3, suppress=True)

        printer.printMatrix(self.cnew)

        print(
            (
                "\n Checking if the modified matrix CNEW is symmetric: i.e. Cij = Cji:  %10s"
                % symmetry.check_symmetric(self.cnew)
            )
        )
