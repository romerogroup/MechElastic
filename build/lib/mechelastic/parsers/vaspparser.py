import numpy as np
import re
from ..utils.constants import *
from ..utils.elements import ELEMENTS, elements_inversed

from ..comms import printer
from ..tests import symmetry

from ..core import ElasticProperties
from ..core import Structure


class VaspOutcar:

    """
    This class contains methods to parse the OUTCAR file.
    """

    def __init__(self, infile="OUTCAR", adjust_pressure=True, verbose=True):

        self.infile = infile
        self.adjust_pressure = adjust_pressure

        self.elastic_tensor = None
        self.compaliance_tensor = None
        self.structure = None
        self.density = None
        self.text = None

        self.mass = None
        self.volume = None
        self.nions = None

        self.lattice_type = None
        self.lattice = None
        self.positions = None
        self.iontype = None
        self.species = None
        self.lattice_constant = None
        self.pressure = None
        self.total_mass = None

        self._parse_outcar(verbose)
        if verbose:
            print(self)
        # self.elastic_properties = ElasticProperties(
        #     self.elastic_tensor, self.structure)#, self.crystal_type)

        return

    def _parse_outcar(self, verbose):
        """
        This function reads outcar and initializes the class variables

        """
        if verbose:
            print("Parsing %s" % self.infile)
        rf = open(self.infile)
        self.text = rf.read()
        rf.close()
        data = self.text
        self.mass = np.array(
            [float(x)
             for x in re.findall("POMASS\s*=\s*([0-9.]*);\s*ZVAL", data)]
        )
        self.volume = float(re.findall(
            "volume of cell\s*:\s*([0-9.]*)", data)[-1])
        self.nions = int(re.findall("NIONS\s*=\s*([0-9]*)", data)[0])
        self.natoms = np.array(
            [
                int(x)
                for x in re.findall("ions per type\s*=\s*([\s0-9]*)", data)[0].split()
            ]
        )

        # Sometimes LATTYP is not present.
        # Only parse if available.
        self.lattice_type = re.findall("LATTYP.*", data)
        if self.lattice_type:
            lattyp = self.lattice_type[-1]
        self.lattice = np.array(
            [
                x.split()[:3]
                for x in re.findall(
                    "direct lattice vectors.*\n([-+0-9.\s\n]*)length of vectors", data
                )[-1].split("\n")[:3]
            ]
        ).astype(float)
        self.positions = np.array(
            [
                x.split()[:3]
                for x in re.findall(
                    "position of ions in fractional coordinates.*\n([-+0-9.\s\n]*)",
                    data,
                )[-1].split("\n")[:self.nions]
            ]
        ).astype(float)
        self.iontype = [
            int(x)
            for x in re.findall("ions per type\s*=\s*([0-9\s]*)", data)[0].split()
        ]
        self.species = re.findall("VRHFIN\s*=([a-zA-Z\s]*):", data)
        self.lattice_constant = [
            float(x)
            for x in re.findall("length of vectors.*\n([0-9.\s]*)", data)[-1].split()[
                :3
            ]
        ]

        # external pressure in kB -> GPa
        self.pressure = float(
            re.findall(r"external\s*pressure\s=\s*([-0-9.]*)\s*kB", data)[-1]
        )
        self.pressure = self.pressure / 10

        atomic_numbers = np.zeros(self.nions, dtype=np.int32)
        k = 0
        for i in range(len(self.iontype)):
            for j in range(self.iontype[i]):
                atomic_numbers[k] = ELEMENTS[self.species[i].strip()]
                k = k + 1

        self.atomic_numbers = atomic_numbers
        self.symbols = [elements_inversed[x] for x in atomic_numbers]

        self.total_mass = 0.0
        for i in range(len(self.natoms)):
            self.total_mass += self.natoms[i] * self.mass[i]

        # converting the units
        self.volume *= 1.0e-30  # from Angstrom to meters
        self.total_mass *= 1.0e-3  # from gram to kg
        self.density = self.total_mass / (self.volume * N_avogadro)
        self.density = self.density

        c = np.matrix(
            [
                x.split()[1:]
                for x in re.findall(
                    "TOTAL ELASTIC MODULI \(kBar\).*\n.*\n.*\n([XYZ0-9.\s-]*)\n\s*-",
                    data,
                )[0].split("\n")
            ][:6]
        ).astype(float)

        # compaliance_tensor = c.I
        if verbose:
            print("\nPrinting Elastic Tensor, Cij as read from OUTCAR\n")

            np.set_printoptions(precision=4, suppress=True)
            printer.print_matrix(c)

        self.elastic_tensor = c.copy()

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

        for j in range(0, 6):
            self.elastic_tensor[3, j] = c[4, j]
            self.elastic_tensor[4, j] = c[5, j]
            self.elastic_tensor[5, j] = c[3, j]

        ctemp = self.elastic_tensor.copy()

        for i in range(0, 6):
            self.elastic_tensor[i, 3] = self.elastic_tensor[i, 4]
            self.elastic_tensor[i, 4] = self.elastic_tensor[i, 5]
            self.elastic_tensor[i, 5] = ctemp[i, 3]

        # Change the units of Cij from kBar to GPa
        self.elastic_tensor /= 10.0

        # Print elastic tensor in GPa units
        if verbose:
            print(
                "\n \n printing Elastic Tensor: corrected order (in GPa units)... \n For example- to generate input for the ELATE code [https://github.com/fxcoudert/elate] \n \n"
            )
            np.set_printoptions(precision=3, suppress=True)
            printer.print_matrix(self.elastic_tensor)

        if self.adjust_pressure and self.pressure != 0.0:
            # In VASP the pressure needs to be adjusted in the elastic tensor.
            # Subtract P in the diagonal and add P in C12, C21, C13, C31, C23, C32.
            if verbose:
                print(
                    "\nAdjusted for pressure since non-zero hydrostatic pressure exists.")
                print("External Pressure : %10.5f GPa" % self.pressure)
                print(
                    "Set flag adjust_pressure=False to disable. -ap 0 in stand-alone mode."
                )

            for i in range(0, 6):
                self.elastic_tensor[i,
                                    i] = self.elastic_tensor[i, i] - self.pressure
            self.elastic_tensor[0,
                                1] = self.elastic_tensor[0, 1] + self.pressure
            self.elastic_tensor[1,
                                0] = self.elastic_tensor[1, 0] + self.pressure
            self.elastic_tensor[0,
                                2] = self.elastic_tensor[0, 2] + self.pressure
            self.elastic_tensor[2,
                                0] = self.elastic_tensor[2, 0] + self.pressure
            self.elastic_tensor[1,
                                2] = self.elastic_tensor[1, 2] + self.pressure
            self.elastic_tensor[2,
                                1] = self.elastic_tensor[2, 1] + self.pressure
            if verbose:
                print("\nPressure adjusted Elastic Tensor in GPa units:\n")
                np.set_printoptions(precision=3, suppress=True)
                printer.print_matrix(self.elastic_tensor)

        self.compaliance_tensor = self.elastic_tensor.I

        if verbose:
            print(
                (
                    "\n Checking if the modified Elastic Tensor is symmetric: i.e. Cij = Cji:  %10s"
                    % symmetry.check_symmetric(self.elastic_tensor)
                )
            )
        self.structure = Structure(self.symbols, self.positions, self.lattice)

        return

    def __str__(self):
        ret = ""
        # cell parameters
        A = float(self.lattice_constant[0])
        B = float(self.lattice_constant[1])
        C = float(self.lattice_constant[2])

        ret += (
            "Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f\n"
            % (A, B, C)
        )

        ret += ("Mass of atoms (in g/mol units): \n" +
                "  ".join(["%s " for x in np.array(self.mass)]) % tuple(self.mass))
        ret += ("Number of atoms: %d\n" % sum(self.natoms))
        ret += ("Total mass (in g/mol): %10.4f \n" % self.total_mass)
        ret += ("Volume of the cell (in Ang^3 units): %10.4f \n" % self.volume)
        ret += ("\nDensity (in kg/m^3 units ): %10.5f" % self.density)
        ret += ("External Pressure (in GPa units ): %10.5f" % self.pressure)
        return ret
