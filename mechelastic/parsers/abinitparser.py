import numpy as np
import re
from itertools import groupby
from ..utils.constants import *
from ..utils.elements import ELEMENTS, elements_inversed

from ..comms import printer
from ..tests import symmetry

from ..core import ElasticProperties
from ..core import Structure


class AbinitOutput:

    """
    This class contains methods to parse the outputs from
    Abinit. It requires the scf/RF output and the ddb analysis
    output.
    """

    def __init__(self, infile="abinit.out", ddbfile="abinit2.out"):

        self.infile = infile
        self.ddbfile = ddbfile

        self.elastic_tensor = None
        self.compaliance_tensor = None
        self.structure = None
        self.text = None

        self._parse_output()

        # self.elastic_properties = ElasticProperties(
        #     self.elastic_tensor, self.structure)#, self.crystal_type)

        return

    def _parse_ddb(self):
        """
        This funcion reads the elastic tensor from the abinit output file
        after running anaddb.
        """

        rf = open(self.ddbfile, "r")
        datamat = rf.read()
        rf.close()

        ET = re.findall(
            r"Elastic\s*Tensor\s*\(relaxed\s*ion\)\s*[\sA-Za-z:\d\^()]*\n([-\s0-9.]*)",
            datamat,
        )
        ET2 = np.array([float(x) for x in ET[0].split()])
        elastic_tensor = ET2.reshape(6, 6)

        return np.matrix(elastic_tensor)

    def _parse_output(self):
        """
        This function reads the abinit output for scf and RF calculations
        and initializes the class variables.

        """

        rf = open(self.infile)
        self.text = rf.read()
        rf.close()
        data = self.text
        mass = np.array(
            [float(x) for x in re.findall(r"amu\s*([E+=0-9.\s]*)\n", data)[0].split()]
        )
        volume = float(
            re.findall("Unit\s*cell\s*volume\s*ucvol\s*=\s*([-+0-9.E\s]*)", data)[-1]
        )
        # converting from Bohr^3 to Angstrom^3
        volume = volume * (0.529177249) ** 3

        nions = int(re.findall(r"\bnatom\b\s*([0-9]*)", data)[-1])

        # counting the number of atoms per type
        typat = [int(x) for x in re.findall(r"\btypat\b\s*([0-9\s]*)", data)[0].split()]
        species_grouped = [i[0] for i in groupby(typat)]
        iontype = [typat.count(i) for i in species_grouped]
        natoms = np.array([typat.count(i) for i in species_grouped])

        # Sometimes LATTYP is not present.
        # Only parse if available.
        lattice_type = re.findall("LATTYP.*", data)
        if lattice_type:
            lattyp = lattice_type[-1]

        full_latt = re.findall(r"\(Bohr,Bohr\^-1\):*\n([-+0-9.RG=\(\)\s\n]*)\n", data)[
            -1
        ].split()
        r1 = [float(x) for x in full_latt[1:4]]
        r2 = [float(x) for x in full_latt[9:12]]
        r3 = [float(x) for x in full_latt[17:20]]
        lattice = np.array((r1, r2, r3), dtype="float64")
        # converting from Bohr to Angstroms
        lattice = 0.529177249 * lattice

        positions = np.array(
            [
                float(x)
                for x in re.findall(r"\bxred\b\s*([0-9-+.\sE]*)\n", data)[-1].split()
            ]
        )
        positions = positions.reshape(nions, 3)

        self.lattice_constant = [
            0.529177249 * float(x)
            for x in re.findall(r"acell\d\s*([-+0-9.E\s]*)Bohr", data)[0].split()
        ]

        # cell parameters
        A = float(self.lattice_constant[0])
        B = float(self.lattice_constant[1])
        C = float(self.lattice_constant[2])

        print(
            "Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f"
            % (A, B, C)
        )

        # external pressure in GPa
        self.pressure = float(re.findall(r"Pressure=\s*([0-9-+.E\s]*)\s*GPa", data)[-1])

        # atomic numbers
        znucl = [
            int(float(x)) for x in re.findall(r"znucl\s*([0-9.\s]*)\n", data)[0].split()
        ]

        species = [elements_inversed[x] for x in znucl]

        atomic_numbers = np.zeros(nions, dtype=np.int32)
        k = 0
        for i in range(len(iontype)):
            for j in range(iontype[i]):
                atomic_numbers[k] = ELEMENTS[species[i]]
                k = k + 1

        print("Atomic numbers")
        print(atomic_numbers)
        atomic_numbers = atomic_numbers
        symbols = [elements_inversed[x] for x in atomic_numbers]
        cell = (lattice, positions, atomic_numbers)

        print("Mass of atoms (in g/mol units): ")
        print((np.array(mass)))
        print("Number of atoms: %d" % np.sum(natoms))

        total_mass = 0.0
        for i in range(len(natoms)):
            total_mass = total_mass + natoms[i] * mass[i]
        print("Total mass (in g/mol): %10.4f " % total_mass)
        print("Volume of the cell (in Ang^3 units): %10.4f " % volume)

        # converting the units
        volume *= 1.0e-30  # from Angstrom to meters
        total_mass *= 1.0e-3  # from gram to kg
        density = total_mass / (volume * N_avogadro)

        print("\nDensity (in kg/m^3 units ): %10.5f" % density)
        print("External Pressure (in GPa units ): %10.5f" % self.pressure)

        c = self._parse_ddb()

        print("\nPrinting Cij matrix as read from Abinit DDB output.\n")
        np.set_printoptions(precision=4, suppress=True)
        printer.printMatrix(c)

        self.elastic_tensor = c.copy()

        self.compaliance_tensor = self.elastic_tensor.I
        # print(self.elastic_tensor)
        print(
            "\n \n printing CNEW: Modified matrix in correct order (in GPa units)... \n For example- to generate input for the ELATE code [https://github.com/fxcoudert/elate] \n"
        )

        # in GPa
        self.elastic_tensor = self.elastic_tensor * 100

        np.set_printoptions(precision=3, suppress=True)

        printer.printMatrix(self.elastic_tensor)

        print(
            (
                "\n Checking if the modified matrix CNEW is symmetric: i.e. Cij = Cji:  %10s"
                % symmetry.check_symmetric(self.elastic_tensor)
            )
        )
        self.structure = Structure(symbols, positions, lattice)

        return
