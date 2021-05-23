# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 00:02:30 2020

@author: lllan
"""


import numpy as np
import re
import math
from ..utils.constants import *
from ..utils.elements import ELEMENTS, elements_inversed
from ..comms import printer
from ..tests import symmetry
from ..core import ElasticProperties
from ..core import Structure

def mag(v):
    return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2])**0.5
def vol(a,b,c):
    return a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0])

class QE_thermo_pw_Parser:

    """
    This class contains methods to parse the OUTCAR file.
    """

    def __init__(self, outfile = 'scf.out', infile = 'scf.in'):


        self.outfile = outfile
        self.infile  = infile

        
        self.elastic_tensor = None
        self.compaliance_tensor = None
        self.structure = None
        self.density = None
                
        self.outText = None
        self.inText  = None
        
        self._parse_files()


        return

    def _parse_files(self):
        """
        This function reads outcar and initializes the class variables

        """


        rf = open(self.outfile)
        self.outText = rf.read()
        rf.close()
        dataOut= self.outText
        
        rf = open(self.infile)
        self.inText = rf.read()
        rf.close()
        dataIn= self.inText
        
        nspecies = int(re.findall("\s*ntyp\s*=\s*(\d)",dataIn)[0])

        

        
        nions = int(re.findall("\s*nat\s*=\s*(\d)",dataIn)[0])
        raw_ions = re.findall("ATOMIC_POSITIONS.*\n" + nions * "(.*)\n",dataIn)[0]
        
        composition = {}
        species_list = []
        mass = []
        raw_species = re.findall("ATOMIC_SPECIES.*\n" + nspecies * "(.*).*\n",dataIn)[0]
        
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
        
        ibrav = int(re.findall("\s*ibrav\s*=\s*(\d*)",dataIn)[0])
        
        if (ibrav == 0):
            num_latticeVectors = 3
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",data)[0])
            raw_cell = re.findall("CELL_PARAMETERS.*\n" + num_latticeVectors * "(.*)\n",data)[0]
            lattice = np.array([[float(component) for component in row.split()] for row in raw_cell])
            lattice = lattice*celldm1
            a = mag(lattice[0,:])
            b = mag(lattice[1,:])
            c = mag(lattice[2,:])
        if (ibrav == 1 or ibrav == 2 or ibrav == 3):
            "cubic systems"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = 1. 
            celldm3 = 1. 
            celldm4 = 0.
            celldm5 = 0.
            celldm6 = 0.
            celldm2 = 1.
            
            if(ibrav == 1):
                lattice = np.array([[celldm1,0,0],
                                    [0,celldm1,0],
                                    [0,0,celldm1]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 2):
                lattice = np.array([[-celldm1/2,0,celldm1/2],
                                    [0,celldm1/2,celldm1/2],
                                    [-celldm1/2,celldm1/2,0]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 3):
                lattice = np.array([[celldm1/2,0,celldm1/2],
                                    [-celldm1/2,celldm1/2,celldm1/2],
                                    [-celldm1/2,-celldm1/2,celldm1/2]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
        if (ibrav == 4):
            "Hexagonal system's must have celldm(1) and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = 1. 
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm4 = 0.
            celldm5 = 0.
            celldm6 =-.5
            
            lattice = np.array([[celldm1,0,0],
                                [-celldm1/2,celldm1*sqrt(3)/2,0],
                                [0,0,celldm1*celldm3]])
            a = mag(lattice[0,:])
            b = mag(lattice[1,:])
            c = mag(lattice[2,:])
            
        if (ibrav == 5):
            "Trigonal or Rhombohedragonal structures must have celldm(1) and celldm(4)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = 1. 
            celldm3 = 1.
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",dataIn)[0])
            if (celldm4 < -1 or celldm4 > 1.):
                sys.exit('\n.... Oops ERROR: "celldm(4)" is WRONG !?!?!?    \n')
            celldm5 = celldm4 
            celldm6 = celldm4
            
            tx = ((1-celldm4)/2)**0.5
            ty = ((1-celldm4)/6)**0.5
            tz = ((1+2*celldm4)/3)**0.5
            
            lattice = np.array([[celldm1*tx,-celldm1*ty,celldm*tz],
                                [0,2*celldm1*ty,tz],
                                [-celldm1*tx,-celldm1*ty,celldm1*tz]])
            a = mag(lattice[0,:])
            b = mag(lattice[1,:])
            c = mag(lattice[2,:])
            
            
            
        if (ibrav == 6 or ibrav == 7):
            "Tetragonal structures must have celldm(1) and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = 1.
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm4 = .0 
            celldm5 = .0
            celldm6 = .0
            
            if(ibrav == 6):
                lattice = np.array([[celldm1,0,0],
                                    [0,celldm1,0],
                                    [0,0,celldm1*celldm3]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])

            elif(ibrav == 7):
                lattice = np.array([[celldm1/2,-celldm1/2,celldm1*celldm3/2],
                                    [celldm1/2,celldm1/2,celldm1*celldm3/2],
                                    [-celldm1/2,-celldm1/2,celldm1*celldm3/2]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
        
        if (ibrav == 8 or ibrav == 9 or ibrav == 10 or ibrav == 11):
            "Orthorhombic structures must have celldm(1), celldm(2), and celldm(3)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm4 = .0
            celldm5 = .0
            celldm6 = .0
            
            if(ibrav == 8):
                lattice = np.array([[celldm1,0,0],
                                    [0,celldm1*celldm2,0],
                                    [0,0,celldm1*celldm3]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 9):
                lattice = np.array([[celldm1/2,celldm1*celldm2/2,0],
                                    [-celldm1/2,celldm1*celldm2/2,0],
                                    [0,0,celldm1*celldm3]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 10):
                lattice = np.array([[celldm1/2,0,celldm1*celldm3/2],
                                    [celldm1/2,celldm1*celldm2/2,0],
                                    [0,celldm1*celldm2/2,celldm1*celldm3/2]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 11):
                lattice = np.array([[celldm1/2,celldm1*celldm2/2,celldm1*celldm3/2],
                                    [-celldm1/2,celldm1*celldm2/2,celldm1*celldm3/2],
                                    [-celldm1/2,-celldm1*celldm2/2,celldm1*celldm3/2]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            
        if (ibrav == 12 or ibrav == 13):
            "Monoclinic structures must have celldm(1), celldm(2), and celldm(3),celldm(4)"
            
            
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm5 = .0
            celldm6 = .0
            
            if(ibrav == 12):
                lattice = np.array([[celldm1/2,0,0],
                                    [celldm1*celldm2*celldm4,celldm1*celldm2*celldm4,0],
                                    [0,0,celldm1*celldm3]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])
            elif(ibrav == 13):
                lattice = np.array([[celldm1/2,0,-celldm1*celldm3/2],
                                    [celldm1*celldm2*celldm4,celldm1*celldm2*celldm4,0],
                                    [celldm1/2,0,celldm1*celldm3/2]])
                a = mag(lattice[0,:])
                b = mag(lattice[1,:])
                c = mag(lattice[2,:])

        if (ibrav == 14):
            "Triclinic structure must have celldm(1), celldm(2), and celldm(3),celldm(4),celldm(5), and celldm(6)"
            celldm1 = float(re.findall("\s*celldm\(1\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm2 = float(re.findall("\s*celldm\(2\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm3 = float(re.findall("\s*celldm\(3\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm4 = float(re.findall("\s*celldm\(4\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm5 = float(re.findall("\s*celldm\(5\)\s*=\s*([\.\d]*)",dataIn)[0])
            celldm6 = float(re.findall("\s*celldm\(6\)\s*=\s*([\.\d]*)",dataIn)[0])
            
            alpha = math.acos(celldm4)
            beta  = math.acos(celldm5)
            gamma = math.acos(celldm6)
           
            lattice = np.array([[celldm1,0,0],
                                [celldm1*celldm2*math.cos(gamma), celldm1*celldm2*math.sin(gamma), 0],
                                [celldm1*celldm3*math.cos(beta),  celldm1*celldm3*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma),
            celldm1*celldm3*math.sqrt( 1 + 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
                     - math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2 )/math.sin(gamma)]])
            
            a = mag(lattice[0,:])
            b = mag(lattice[1,:])
            c = mag(lattice[2,:])
            
        
        volume = vol(lattice[0,:],lattice[1,:],lattice[2,:])
        
        # Convert to a.u. (bohr) to Angtroms
        lattice = lattice*(0.529177249)   
        self.lattice_constant = [a*0.529177249  ,b*0.529177249  ,c*0.529177249  ]         


        # cell parameters
        A = float(self.lattice_constant[0])
        B = float(self.lattice_constant[1])
        C = float(self.lattice_constant[2])

        print(
            "Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f"
            % (A, B, C)
        )

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
        print("Number of atoms: %d" % np.sum(nions))

        total_mass = 0.0
        for i in range(len(ions)):
            total_mass = total_mass +  mass[i]
        print("Total mass (in g/mol): %10.4f " % total_mass)
        volAng= volume*0.529177249**3
        print("Volume of the cell (in Ang^3 units): %10.4f " % volAng)

        # converting the units
        volume *= (5.29177249**-11)**3  # from atomic units (bohr)  to meters
        total_mass *= 1.0e-3  # from gram to kg
        density = total_mass / (volume * N_avogadro)
        self.density = density 
        print("\nDensity (in kg/m^3 units ): %10.5f" % density)


        num_rows = 6
        raw_stiffness = re.findall("\s*Elastic\sconstants.*\n.*\n" + num_rows * "\s*(.*)\n",dataOut)[0]
        c = np.array([[float(x)/10.0 for x in row.split()[1:]] for row in raw_stiffness])

    
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
