# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 00:02:30 2020

@author: lllan
"""

#import numpy as np
#import re
#import math
#from elements import ELEMENTS, elements_inversed
#from constants import N_avogadro
#import symmetry
import numpy as np
import re
from ..utils.constants import *
from ..utils.elements import ELEMENTS, elements_inversed

from ..comms import printer
from ..tests import symmetry

from ..core import ElasticProperties
from ..core import Structure

class QEParser:

    """
    This class contains methods to parse the OUTCAR file.
    """

    def __init__(self, outfile = 'ElaStic_2nd.out',scfFile = 'scf.in'):

        self.infile = "ElaStic_PW.in"
        self.outfile = outfile
        self.scfFile = scfFile
        self.infoElaStic = 'INFO_ElaStic'
        
        self.elastic_tensor = None
        self.compaliance_tensor = None
        self.structure = None
        
        self.outText = None
        self.inText = None
        self.scfText = None
        self.infoElaSticText = None
        
        self.test = None
        self._parse_files()

        # self.elastic_properties = ElasticProperties(
        #     self.elastic_tensor, self.structure)#, self.crystal_type)

        return

    def _parse_files(self):
        """
        This function reads outcar and initializes the class variables

        """


        rf = open(self.outfile)
        self.outText = rf.read()
        rf.close()
        dataOut = self.outText
        
        rf = open(self.infile)
        self.inText = rf.read()
        rf.close()
        data = self.inText
        
        rf = open(self.scfFile)
        self.scfText = rf.read()
        rf.close()
        dataScf = self.scfText
        
        rf = open(self.infoElaStic)
        self.infoElaSticText = rf.read()
        rf.close()
        dataInfo = self.infoElaSticText
        
        nspecies = int(re.findall("\s*ntyp\s*=\s*(\d)",data)[0])
#        mass = np.array(
#            [float(x) for x in re.findall("ATOMIC_SPECIES\n"+ nspecies *"\s*\w*\s*([\.\d]*).*\n", data)]
#        )
        
        nions = int(re.findall("\s*nat\s*=\s*(\d)",data)[0])
        raw_ions = re.findall("ATOMIC_POSITIONS.*\n" + nions * "(.*)\n",data)[0]
        
        composition = {}
        species_list = []
        mass = []
        raw_species = re.findall("ATOMIC_SPECIES.*\n" + nspecies * "(.*).*\n",data)[0]
        
        if(nspecies == 1):
            composition[raw_species.split()[0]] = {'numAtom':0, 'mass':raw_species.split()[1]}
            species_list.append(raw_species.split()[0])
            
        else:
            for nspec in range(nspecies):
                species_list.append(raw_species[nspec].split()[0])
                composition[raw_species.split()[0]] = {'numAtom':0, 'mass':raw_species.split()[1]}
                
        if(nions == 1):
            composition[raw_species.split()[0]]['numAtom'] = 1
        else:
            for ions in range(nions):
                for species in range(len(species_list)):
                    if(raw_ions[ions].split()[0] == species_list[species]):
                        composition[raw_ions[ions].split()[0]]['numAtom'] += 1
                
        volume = float(re.findall("Volume of equilibrium unit cell\s*=\s*([0-9.]*)[\s\S]*", dataInfo)[-1])
        
        
        
#
#        # Sometimes LATTYP is not present.
#        # Only parse if available.
#        lattice_type = re.findall("LATTYP.*", data)
#        if lattice_type:
#            lattyp = lattice_type[-1]
#        cell_parameters = np.array(
#            [
#                x.strip().split()
#                for x in re.findall("CELL_PARAMETERS \(alat\)\n" + 3 * "(.*)\n", data)[0]
#            ]
#        ).astype(float)
      
        positions = np.array(
            [
                x.strip().split()[1:]
                for x in raw_ions
            ]
        ).astype(float)
        
        
        iontype = [int(composition[x]['numAtom']) for x in list(composition.keys())]
        
        ions = [x.strip().split()[0] for x in raw_ions]
        mass  = [float(composition[x]['mass']) for x in ions]
        
        species = list(composition.keys())
        
        ibrav = int(re.findall("\s*ibrav\s*=\s*(\d*)",dataScf)[0])
        
        
        if (ibrav == 1 or ibrav == 2 or ibrav == 3):
            "cubic systems"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = 1. 
            celldm3 = 1. 
            celldm4 = 0.
            celldm5 = 0.
            celldm6 = 0.
            celldm2 = 1.
            
            if(ibrav == 1):
                a = celldm1
                b = celldm1
                c = celldm1
                lattice = np.array([[celldm1,0,0],[0,celldm1,0],[0,0,celldm1]])
            elif(ibrav == 2):
                a = celldm1*(2**0.5)/2
                b = celldm1*(2**0.5)/2
                c = celldm1*(2**0.5)/2
                lattice = np.array([[-celldm1/2,0,celldm1/2],[0,celldm1/2,celldm1/2],[-celldm1/2,celldm1/2,0]])
            elif(ibrav == 3):
                a = celldm1*(3**0.5)/2
                b = celldm1*(3**0.5)/2
                c = celldm1*(3**0.5)/2
                lattice = np.array([[celldm1/2,0,celldm1/2],[-celldm1/2,celldm1/2,celldm1/2],[-celldm1/2,-celldm1/2,celldm1/2]])
            
        if (ibrav == 4):
            "Hexagonal system's must have celldm(1) and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = 1. 
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",data)[0])
            celldm4 = 0.
            celldm5 = 0.
            celldm6 =-.5
            
            
            a = celldm1
            b = celldm1
            c = a*celldm3
            
        if (ibrav == 5):
            "Trigonal or Rhombohedragonal structures must have celldm(1) and celldm(4)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = 1. 
            celldm3 = 1.
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",data)[0])
            if (celldm4 < -1 or celldm4 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(4)" is WRONG !?!?!?    \n')
            celldm5 = celldm4 
            celldm6 = celldm4
            
            a = celldm1
            b = celldm1*((2+celldm4)/3)**0.5
            c = celldm1 
            
            
            
        if (ibrav == 6 or ibrav == 7):
            "Tetragonal structures must have celldm(1) and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = 1.
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",data)[0])
            celldm4 = .0 
            celldm5 = .0
            celldm6 = .0
            
            if(ibrav == 6):
                a = celldm1
                b = celldm1
                c = a*celldm3
            elif(ibrav == 7):
                a = (2*(celldm1/2)**2 + (celldm1*celldm3/2)**2)**0.5 
                b = (2*(celldm1/2)**2 + (celldm1*celldm3/2)**2)**0.5 
                c = (2*(celldm1/2)**2 + (celldm1*celldm3/2)**2)**0.5 
        
        if (ibrav == 8 or ibrav == 9 or ibrav == 10 or ibrav == 11):
            "Orthorhombic structures must have celldm(1), celldm(2), and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",data)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",data)[0])
            celldm4 = .0
            celldm5 = .0
            celldm6 = .0
            
            if(ibrav == 8):
                a = celldm1
                b = a*celldm2
                c = a*celldm3
            elif(ibrav == 9):
                a = ((celldm1**2+(celldm1*celldm2)**2)**0.5)/2
                b = ((celldm1**2+(celldm1*celldm2)**2)**0.5)/2
                c = celldm1*celldm3
            elif(ibrav == 10):
                a = ((celldm1**2+(celldm1*celldm3)**2)**0.5)/2
                b = ((celldm1**2+(celldm1*celldm2)**2)**0.5)/2
                c = (((celldm1*celldm2)**2+(celldm1*celldm3)**2)**0.5)/2
            elif(ibrav == 10):
                a = ((celldm1**2 + (celldm1*celldm2)**2 + (celldm1*celldm3)**2 )**0.5 )/2
                b = ((celldm1**2 + (celldm1*celldm2)**2 + (celldm1*celldm3)**2 )**0.5 )/2
                c = ((celldm1**2 + (celldm1*celldm2)**2 + (celldm1*celldm3)**2 )**0.5 )/2
            
        if (ibrav == 12 or ibrav == 13):
            "Monoclinic structures must have celldm(1), celldm(2), and celldm(3),celldm(4)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",data)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",data)[0])
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",data)[0])
            celldm5 = .0
            celldm6 = .0
            
            
            if(ibrav == 12):
                a = celldm1
                b = celldm1 * celldm2
                c = celldm1 * celldm3
            elif(ibrav == 13):
                a = ((celldm1**2+(celldm1*celldm3)**2)**0.5)/2
                b = celldm1 * celldm2
                c = ((celldm1**2+(celldm1*celldm3)**2)**0.5)/2

        if (ibrav == 14):
            "Triclinic structure must have celldm(1), celldm(2), and celldm(3),celldm(4),celldm(5), and celldm(6)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",data)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",data)[0])
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",data)[0])
            celldm5 = float(re.findall("\s*celldm\(5\)\s*=\s*([\.\d]*)",data)[0])
            celldm6 = float(re.findall("\s*celldm\(6\)\s*=\s*([\.\d]*)",data)[0])
            
            alpha = math.acos(celldm4)
            beta  = math.acos(celldm5)
            gamma = math.acos(celldm6)
            a = celldm1
            b = celldm1*celldm2
            c = celldm3*celldm1*((math.cos(beta))**2+ ((1/math.sin(gamma)) * (math.cos(alpha)-math.cos(beta)*math.cos(gamma)))**2 + ((1/math.sin(gamma)**2)*(1+ 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2)))**0.5

            lattice = np.array([[a,0,0],[b*math.cos(gamma), b*math.sin(gamma), 0],[c*math.cos(beta),  c*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma),
            c*sqrt( 1 + 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
                     - math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2 )/math.sin(gamma)]])
            
        self.lattice_constant = [a,b,c]            
        self.test = lattice 
#
        # cell parameters
        A = float(self.lattice_constant[0])
        B = float(self.lattice_constant[1])
        C = float(self.lattice_constant[2])

        print(
            "Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f"
            % (A, B, C)
        )

#        # external pressure in kB -> GPa
#        self.pressure = float(
#            re.findall(r"external\s*pressure\s=\s*([-0-9.]*)\s*kB", data)[-1]
#        )
#        self.pressure = self.pressure / 10

        atomic_numbers = np.zeros(nions, dtype=np.int32)
        k = 0
        for i in range(len(iontype)):
            for j in range(iontype[i]):
                atomic_numbers[k] = ELEMENTS[species[i]]
                k = k + 1
        self.test = atomic_numbers

        print("Atomic numbers")
        print(atomic_numbers)
        atomic_numbers = atomic_numbers
        symbols = [elements_inversed[x] for x in atomic_numbers]
        cell = (lattice, positions, atomic_numbers)

        print("Mass of atoms (in g/mol units): ")
        print((np.array(mass)))
        print("Number of atoms: %d" % np.sum(nions))

        total_mass = 0.0
        for i in range(len(ions)):
            total_mass = total_mass +  mass[i]
        print("Total mass (in g/mol): %10.4f " % total_mass)
        volAng= volume*0.529177249**3
        print("Volume of the cell (in Ang^3 units): %10.4f " % volAng)

        # converting the units
        volume *= (5.29177249**-11)**3  # from Angstrom to meters
        total_mass *= 1.0e-3  # from gram to kg
        density = total_mass / (volume * N_avogadro)

        print("\nDensity (in kg/m^3 units ): %10.5f" % density)
#        print("External Pressure (in GPa units ): %10.5f" % self.pressure)


        raw_stiffness = re.findall("(?<=Elastic constant \(stiffness\) matrix in GPa\:)([\s\S]*?)(?=Elastic compliance)",dataOut)[0].strip().split()
        c = [raw_stiffness[x:x+6] for x in range(0, len(raw_stiffness), 6)]
        c = np.array(c).astype(float)
    
        print("\nPrinting Cij matrix as read from outfile (GPa)\n")
        self.elastic_tensor = c.copy()
        np.set_printoptions(precision=3, suppress=True)
        printer.printMatrix(self.elastic_tensor)
        
        self.compaliance_tensor = np.linalg.inv(self.elastic_tensor)


        print(
            (
                "\n Checking if the modified matrix CNEW is symmetric: i.e. Cij = Cji:  %10s"
                % symmetry.check_symmetric(self.elastic_tensor)
            )
        )
        self.structure = Structure(symbols, positions, lattice)

        return
