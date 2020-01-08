#!/usr/bin/env python

'''This script reads elastic constants from OUTCAR, calculates elastic moduli, and performs mechanical stability test: 
   Authors: Sobhit Singh (1,2) and Aldo Romero (2)
   Email: smsingh@mix.wvu.edu, alromero@mail.wvu.edu
   (1) Rutgers University, Piscataway, NJ, USA
   (2) West Virginia University, Morgantown, WV, USA
Version: 12.2019   #Dec. 2019

Please cite the below paper if you use this code for your research: 
     Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018).


To RUN this code for a 3D system assuming you know the crystal type before hand: 
 >> python MechElastic.py -i OUTCAR-file -c hexagonal --dim 3D

Note: The crystal type is needed only to perform the elastic stabilty test. 
      This test is not yet implemented for 2D systems, so you can ignore the crystal type for 2D systems.  

If the crystal symmetry is not provided by user, then the code will rely on spglib to find it
 >> python MechElastic.py -i OUTCAR-file --dim 3D 

To RUN this code for a 2D monolayer system 
 >> python MechElastic.py -i OUTCAR-file --dim 2D
      (please pay attention to the warning for 2D systems)

OUTCAR file, if present in the current directory, is read as default unless a filename is specified by the user.


Disclaimer: Please check the authenticity of your results before publishing. 
            AUTHORS of this script do not guarantee the quality and/or accuracy 
            of results generated using this script. 

'''

import numpy as np
import spglib
import sys
from numpy import linalg as LA
from prettytable import PrettyTable
import argparse


mechelastic_version = 12.2019  #Dec. 2019
def print_version(mechelastic_version):
    print(('Version: %s \n' % mechelastic_version))


#created at http://www.network-science.de/ascii/.
def print_mechelastic():
    mechelastic_str = """ 
      /\/\   ___  ___| |__   /__\ | __ _ ___| |_(_) ___
     /    \ / _ \/ __| '_ \ /_\ | |/ _` / __| __| |/ __|
    / /\/\ \  __/ (__| | | //__ | | (_| \__ \ |_| | (__
    \/    \/\___|\___|_| |_\__/ |_|\__,_|___/\__|_|\___| """
    print(("  \n" + mechelastic_str))

    return


def printMatrix(c):
    p = PrettyTable()
    for row in c:
        p.add_row(row)
    print(p.get_string(header=False, border=False))
    
    return



def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)    
   
    
def positive_evals(cnew):
    print("\nEigen Values of the matrix CNEW:")
    evals = LA.eigvals(cnew)
    print(evals)
    check = 0
    for i in range(len(evals)):
        if evals[i] > 0.0:
            pass
        else:
            check = 1
    if check == 1:
        print("ATTENTION: One or more eigen values are negative indicating elastic instability.")
    if check == 0:
        print("All eigen values are positive indicating elastic stability.")  
        
    return


            
print_mechelastic()
print_version(mechelastic_version)


## Defining useful constants
## Avogadro number
N_avogadro = 6.022140857e+23
## Planck's constant
h_Planck = 6.62607004e-34
## Boltzmann constant J/K
kB = 1.38064852e-23
## Atomic mass unit. kg
amu = 1.66053886e-27 

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, 
                    help="input the OUTCAR", default='OUTCAR')
parser.add_argument("-c", "--crystal", type=str,
                    help="Provide the crystal type. Otherwise it would be determined from OUTCAR") 
parser.add_argument("-d", "--dim", type=str,
                    help="Enter the dimension, 2D or 3D: For example: '-d 2D' or '--dim 2D'. Default is '3D' ", default='3D')
args = parser.parse_args()
args.input

print("-----------------------------")
print("List of arguments entered:") 
print("Input file name:", args.input) 
print("Crystal type: ", args.crystal)
print("Dimensions:", args.dim)
print("-----------------------------")


# following info about elements has been taken from https://www.lfd.uci.edu/~gohlke/code/elements.py.html

ELEMENTS = {'H' : 1 , 'He' : 2 , 'Li' : 3 , 'Be' : 4 , 'B' : 5 , 'C' : 6 , 'N' : 7 , 'O' : 8 , 'F' : 9 , 'Ne' : 10 , 'Na' : 11 , 'Mg' : 12 , 'Al' : 13 , 'Si' : 14 ,
'P' : 15 , 'S' : 16 , 'Cl' : 17 , 'Ar' : 18 , 'K' : 19 , 'Ca' : 20 , 'Sc' : 21 , 'Ti' : 22 , 'V' : 23 , 'Cr' : 24 , 'Mn' : 25 , 'Fe' : 26 , 'Co' : 27 , 'Ni' : 28 ,
'Cu' : 29 , 'Zn' : 30 , 'Ga' : 31 , 'Ge' : 32 , 'As' : 33 , 'Se' : 34 , 'Br' : 35 , 'Kr' : 36 , 'Rb' : 37 , 'Sr' : 38 , 'Y' : 39 , 'Zr' : 40 , 'Nb' : 41 , 'Mo' : 42 ,
'Tc' : 43 , 'Ru' : 44 , 'Rh' : 45 , 'Pd' : 46 , 'Ag' : 47 , 'Cd' : 48 , 'In' : 49 , 'Sn' : 50 , 'Sb' : 51 , 'Te' : 52 , 'I' : 53 , 'Xe' : 54 , 'Cs' : 55 , 'Ba' : 56 ,
'La' : 57 , 'Ce' : 58 , 'Pr' : 59 , 'Nd' : 60 , 'Pm' : 61 , 'Sm' : 62 , 'Eu' : 63 , 'Gd' : 64 , 'Tb' : 65 , 'Dy' : 66 , 'Ho' : 67 , 'Er' : 68 , 'Tm' : 69 , 'Yb' : 70 ,
'Lu' : 71 , 'Hf' : 72 , 'Ta' : 73 , 'W' : 74 , 'Re' : 75 , 'Os' : 76 , 'Ir' : 77 , 'Pt' : 78 , 'Au' : 79 , 'Hg' : 80 , 'Tl' : 81 , 'Pb' : 82 , 'Bi' : 83 , 'Po' : 84 ,
'At' : 85 , 'Rn' : 86 , 'Fr' : 87 , 'Ra' : 88 , 'Ac' : 89 , 'Th' : 90 , 'Pa' : 91 , 'U' : 92 , 'Np' : 93 , 'Pu' : 94 , 'Am' : 95 , 'Cm' : 96 , 'Bk' : 97 , 'Cf' : 98 ,
'Es' : 99 , 'Fm' : 100 , 'Md' : 101 , 'No' : 102 , 'Lr' : 103 , 'Rf' : 104 , 'Db' : 105 , 'Sg' : 106 , 'Bh' : 107 , 'Hs' : 108 , 'Mt' : 109 } 




f = open(args.input,'r')
lines = f.readlines()
f.close()

# s is the complaiance matrix and c is the elastic constant matrix
s=np.zeros((6,6))
c=np.zeros((6,6))

# parsing OUTCAR
nions=-1
mass=[]
foundcell=0
icell=-1
iatom=-1
lattice=[]
positions=[]
atomic_numbers=[]
symbols=[]
for i in lines:
    word = i.split()
    if "TOTAL ELASTIC MODULI (kBar)" in i:
        ii=lines.index(i)
    if "POMASS" in i and "ZVAL" in i:
        m=i.split()[2]
        mass.append(float(m[:-1]))
    if "volume of cell" in i:
        vol=float(i.split()[4])
    if "NIONS =" in i:
        nions=int(i.split()[11])
    if "ions per type =" in i:
        natoms=list(map(int,i.split()[4:]))
    if "LATTYP" in i:
        lattyp=i
    if (icell > -1 and not foundcell):
        lattice.append(list(map(float,i.split()))[:3])
        icell=icell+1
        if (icell == 3):
            foundcell=1
    if not foundcell and "direct lattice vectors" in i:
        icell=icell+1
    if (nions>0 and iatom> -1):
        positions.append(list(map(float,word)))
        iatom=iatom+1
        if (iatom==nions):
            iatom=-1
    if nions>0 and "position of ions in fractional coordinates" in i:
        iatom=iatom+1
    if "ions per type = " in i:
        iontype=list(map(int,word[4:]))
    if "VRHFIN" in i:
        symbols.append(i.split()[1][1:-1])
    if "length of vectors" in i:
        n=lines.index(i)
        l=lines[n+1]
        lattconst = l.split()

    #   break

##print lattice
##print positions

## cell parameters
A = float(lattconst[0])
B = float(lattconst[1])
C = float(lattconst[2])


print("Lattice parameters (in Angs.): a = %10.5f      b = %10.5f     c = %10.5f" %(A, B, C))
atomic_numbers=np.zeros(nions,dtype=np.int32)
k=0
for i in range(len(iontype)):
    for j in range(iontype[i]):
        atomic_numbers[k]=ELEMENTS[symbols[i]]
        k=k+1

print("Atomic numbers")
print(atomic_numbers)
lattice=np.array(lattice)
positions=np.array(positions)
cell=(lattice,positions,atomic_numbers)

mass=np.array(mass)
print("Mass of atoms (in g/mol units): ")
print((np.array(mass)))
print("Number of atoms: %d" %np.sum(natoms))

totalmass=0.0
for i in range(len(natoms)):
        totalmass = totalmass + natoms[i]*mass[i]
print("Total mass (in g/mol): %10.4f " %totalmass)
print("Volume of the cell (in Ang^3 units): %10.4f " %vol)

## converting the units
vol = vol*1.0e-30     ## from Angstrom to meters
totalmass = totalmass*1.0e-3 ## from gram to kg
density = totalmass/(vol*N_avogadro)

print('\n Density (in kg/m^3 units ): %10.5f' % density)


for i in range(0,6):
    l=lines[ii+3+i]
    word = l.split()
    s[i][:]=word[1:7]
    c[i]=list(map(float,s[i][:]))

#print c

mc=np.matrix(c)
mci=mc.I

for i in range(0,6):
    for j in range(0,6):
        s[i][j]=mci[i,j]

print("Printing Cij matrix as read from OUTCAR\n")
printMatrix(c)

## Redefining the Cij matrix into the correct Voigt notation since VASP's OUTCAR has a different order
## In VASP: Columns and rows are listed as: 1, 2, 3, 6, 4, 5
## In this format OUTCAR's C44 values would be actually C66, C55 would be C44, and C66 would be C55. 
## OUTCAR follows the below order: 
## [C11 C12 C13 C16 C14 C15]
## [C21 C22 C23 C26 C24 C25]
## [C31 C32 C33 C36 C34 C35]
## [C61 C62 C63 C66 C64 C65]
## [C41 C42 C43 C46 C44 C45]
## [C51 C52 C53 C56 C54 C55] 


cnew = np.zeros((6,6))
snew = np.zeros((6,6))

cnew = np.copy(c)

#for i in range(0,3):
#   for j in range(0,6):
#     cnew[i][j]=c[i][j]

for j in range(0,6):
    cnew[3][j] = c[4][j]
    cnew[4][j] = c[5][j]
    cnew[5][j] = c[3][j]

#print "CNEW"
#print(matrix(cnew))

ctemp=np.zeros((6,6))
ctemp = np.copy(cnew)

for i in range(0,6):
    cnew[i][3] = cnew[i][4]
    cnew[i][4] = cnew[i][5]
    cnew[i][5] = ctemp[i][3]


### Change the units of Cij from kBar to GPa
for i in range(0,6):
    for j in range(0,6):
        cnew[i][j] = cnew[i][j]/10.0

#print "CTEMP "
#print(matrix(ctemp))

print("\n \n printing CNEW: Modified matrix in correct order (in GPa units)... \n For example- to generate input for the ELATE code [https://github.com/fxcoudert/elate] \n")
#np.set_printoptions(precision=3)
np.set_printoptions(precision=3, suppress=True)
#print(np.matrix(cnew))


printMatrix(cnew)
print(('\n Checking if the modified matrix CNEW is symmetric: i.e. Cij = Cji:  %10s' % check_symmetric(cnew)))



        
def elastic_const_bulk(cnew, snew):
    ##Bulk: Voigt
    KV = (cnew[0][0] + cnew[1][1] + cnew[2][2]) + 2*(cnew[0][1] + cnew[1][2] + cnew[2][0])
    KV = KV/9.0

    ## Shear: Voigt
    GV = (cnew[0][0] + cnew[1][1] + cnew[2][2]) - (cnew[0][1] + cnew[1][2] + cnew[2][0]) + 3*(cnew[3][3] + cnew[4][4] + cnew[5][5])
    GV = GV/15.0

    ## Young's: Voigt
    EV = (9*KV*GV)/(3*KV + GV)

    ## Poisson's ratio: Voigt
    Nu_V = (3*KV - EV)/(6*KV)

    ## P-wave modulus, M: Voigt
    MV = KV + (4*GV/3.0)

    #####  Reuss Method
    mc_new = np.matrix(cnew)
    mci_new = mc_new.I

    for i in range(0,6):
        for j in range(0,6):
            snew[i][j]=mci_new[i,j]

    ## bulk: Reuss
    KR = (snew[0][0] + snew[1][1] + snew[2][2]) + 2*(snew[0][1] + snew[1][2] + snew[2][0])
    KR = 1.0/KR

    ## Shear: Reuss
    GR = 4*(snew[0][0] + snew[1][1] + snew[2][2]) -4*(snew[0][1] + snew[1][2] + snew[2][0]) + 3*(snew[3][3] + snew[4][4] + snew[5][5])
    GR = 15.0/GR

    ## Young's: Reuss
    ER = (9*KR*GR)/(3*KR + GR)

    ## Poisson's ratio: Reuss
    Nu_R = (3*KR - ER)/(6*KR)

    ## P-wave modulus, M: Reuss
    MR = KR + (4*GR/3.0)


    ### Voigt-Reuss-Hill Approximation: average of both methods
    ## VRH
    Kvrh = (KV + KR)/2.0
    Gvrh = (GV + GR)/2.0
    Evrh = (EV + ER)/2.0
    Nu_vrh = (Nu_V + Nu_R)/2.0
    Mvrh = (MV + MR)/2.0
    KG_ratio_V = KV/GV
    KG_ratio_R = KR/GR
    KG_ratio_vrh = Kvrh/Gvrh    
    ## Elastic Anisotropy
    ## Zener anisotropy for cubic crystals only
    Az = 2*cnew[3][3]/(cnew[0][0] - cnew[0][1])

    ##Ranganathan and Ostoja-Starzewski method: Phys. Rev. Lett. 101, 055504 (2008).
    ## for any crystalline symmetry: Universal anisotropy index
    AU = (KV/KR) + 5*(GV/GR) - 6.0

    #""'Note: AU is a relative measure of anisotropy with respect to a limiting value. For example, AU does not prove that a crystal having AU = 3 has double the anisotropy of another crystal with AU = 1.5. I""'

    ## log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016)
    #""'AL  CV , CR   is based on the distance between the averaged stiffnesses CV and
    #CR , which is more appropriate. Clearly, AL  CV , CR   is zero when the crystallite is isotropic.

    AL = np.sqrt(5)*2.303*np.log(1 + (AU/5))
    

    print("\n \n                         Voigt     Reuss    Average")
    print("-------------------------------------------------------")
    print("Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (KV, KR, Kvrh))
    print("Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (GV, GR, Gvrh))
    print("Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (EV, ER, Evrh))
    print("Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh))
    print("P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh))
    print("Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) " % (KG_ratio_V, KG_ratio_R, KG_ratio_vrh,  ductile_test(KG_ratio_vrh) ))
    print("-------------------------------------------------------")


    print(" \n \n \t \t Elastic Anisotropy \n ")
    print("Zener anisotropy (true for cubic crystals only) Az = %10.5f" %Az)
    print("Universal anisotropy index (Ranganathan and Ostoja-Starzewski method; PRL 101, 055504 (2008)) Au = %10.5f" %AU)
    print("Log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016) AL = %10.5f " %AL)


    ## Calculation of elastic wave velocities and Debye temperature using the values obtained from Voigt-Reuss-Hill Approximation
    G = Gvrh*1.0e+9   ## converting from GPa to Pascal units (kg/ms^2)
    K = Kvrh*1.0e+9

    ## transverse velocity: Navier's equation
    vt = np.sqrt((G/density))

    ## longitudinal velocity: Navier's equation
    vl = np.sqrt(((3*K + 4*G)/(3.0*density)))

    ## average 
    vm=1.0/(np.cbrt((1./3.)*(2./(vt*vt*vt)+1./(vl*vl*vl))))


    ## Debye temperature calculated using  Orson Anderson's proposal [Ref- J. Phys. Chem. Solids (1963) 24, 909-917]
    q = np.sum(natoms)	
    theta = (h_Planck/kB)*vm*np.cbrt((3*q*N_avogadro*density)/(4*(np.pi)*totalmass))

    ## Melting temperature estimated using empirical relation from Ref: Johnston I, Keeler G, Rollins R and Spicklemire S 1996 
    ## Solid State Physics Simulations, The Consortium for Upper-Level Physics Software (New York: Wiley)
    Tm = 607 + 9.3*Kvrh

    print("\n \t \t  Elastic wave velocities calculated using Navier's equation  (in m/s units) \n")
    print("----------------------------------------------- ")
    print("Longitudinal wave velocity (vl) : %10.5f " %vl)
    print("Transverse wave velocity (vt) : %10.5f " %vt)
    print("Average wave velocity (vm) : %10.5f " %vm)
    print("Debye temperature  (in K) : %10.5f " %theta)
    print("WARNING: Debye model for the atomic displacement is based on a monoatomic crystal, here we consider an average mass if your crystal has several species")
#    print "Atomic displacement at 150, 300 and 450 K  (in A^2) : %10.5f %10.5f %10.5f" %(u2FromDebye(mass,natoms,theta,150.),u2FromDebye(mass,natoms,theta,300.), u2FromDebye(mass,natoms,theta,450.))
    print("\nMelting temperature calculated from empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K)")
    print("Tm (in K)=  %10.5f (plus-minus 555 K) " % Tm)
    print("----------------------------------------------- ")

    
    
    positive_evals(cnew)

    crystal = np.array(['cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'])

    if args.crystal is not None:
            crystaltype = args.crystal
            print("\n \t \t Mechanical Stability Test \n")
            stability_test(cnew, crystaltype)
        
    else:
            print("\n")
            print("WARNING: crystal symmetry class  was not provided by user, it will be taken from the OUTCAR")
            print("One of the following was expected as the second argument: \n 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'")
            crystaltype=''
            spg=int(spglib.get_spacegroup(cell, symprec=1e-5).split()[1][1:-1])
            if (spg>=1 and spg<=2):
                   crystaltype='triclinic'
            if (spg>=3 and spg<=15):
                  crystaltype='monoclinic'
            if (spg>=16 and spg<=74):
                  crystaltype='orthorhombic'
            if (spg>=75 and spg<=142):
                  crystaltype='tetragonal'
            if (spg>=143 and spg<=167):
                if (spg == 155 or spg == 160 or spg == 166):
                      crystaltype='rhombohedral-1'
                else:
                      crystaltype='rhombohedral-2'
            if (spg>=168 and spg<=194):
                  crystaltype='hexagonal'
            if (spg>=195):
                  crystaltype='cubic'
          
            print('From OUTCAR the crystal type is = ', crystaltype)
            stability_test(cnew, crystaltype)        
        
    
    return



def ductile_test(ratio):
    if(ratio > 1.75):
        return "ductile"
    else:
        return "brittle"


    
def stability_test(matrix, crystaltype):
    c = np.copy(matrix)

    if(crystaltype =="cubic"):
        print("Cubic crystal system \n")
        print("Born stability criteria for the stability of cubic system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n ")
        
        ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if(c[0][0] - c[0][1] > 0.0):
            print("Condition (i) satified.")
        else:
            print("Condition (i) NOT satisfied.")
            
        if(c[0][0] + 2*c[0][1] > 0.0):
            print("Condition (ii) satified.")
        else:
            print("Condition (ii) NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iii) satified.")
        else:
            print("Condition (iii) NOT satisfied.")


    if(crystaltype =="hexagonal"):
        print("Hexagonal crystal system \n")
        print("Born stability criteria for the stability of hexagonal system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n ")

        ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if(c[0][0] - c[0][1] > 0.0):
            print("Condition (i) satified.")
        else:
            print("Condition (i) NOT satisfied.")

        if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
            print("Condition (ii) satified.")
        else:
            print("Condition (ii) NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iii) satified.")
        else:
            print("Condition (iii) NOT satisfied.")


    if(crystaltype =="tetragonal"):
        print("Tetragonal crystal system \n")
        print("Born stability criteria for the stability of Tetragonal system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0;   (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n ")

        ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if(c[0][0] - c[0][1] > 0.0):
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if(c[5][5] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if(2*c[0][5]*c[0][5] < c[5][5]*(c[0][0] - c[0][1])):
            print("Condition (v) is satified.")
        else:
            print("Condition (v) is NOT satisfied.")

            
    if(crystaltype =="rhombohedral-1"):
        print("Rhombohedral (class-1): i.e. structures with point group: 3m, -3m and 32 \n")
        print("Born stability criteria for the stability of Rhombohedral-1 class system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n ")

        if(c[0][0] - c[0][1] > 0.0):
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if((c[0][2]*c[0][2]) < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if(c[0][3]*c[0][3] < 0.5*c[3][3]*(c[0][0] - c[0][1])):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")



    if(crystaltype =="rhombohedral-2"):
        print("Rhombohedral (class-2): i.e structures with point group: 3, -3 \n")
        print("Born stability criteria for the stability of Rhombohedral-1 class system are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 + C15^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0;  Note: C15 is added.. \n ")

        if(c[0][0] - c[0][1] > 0.0):
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if(c[0][2]*c[0][2] < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if(c[0][3]*c[0][3] + c[0][4]*c[0][4]  < 0.5*c[3][3]*(c[0][0] - c[0][1])):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")


    if(crystaltype =="orthorhombic"):
        print("Orthorhombic crystal system.... \n")
        print("Born stability criteria for the stability of Orthorhombic systems are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n")
        print("(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0;   (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n")

        ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
        if(c[0][0]  > 0.0):
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satisfied.")

        if(c[0][0]*c[1][1] < c[0][1]*(c[0][1])):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satisfied.")

        if(c[0][0]*c[1][1]*c[2][2] + 2*c[0][1]*c[0][2]*c[1][2] - c[0][0]*c[1][2]*c[1][2] - c[1][1]*c[0][2]*c[0][2] - c[2][2]*c[0][1]*c[0][1] > 0 ):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satisfied.")

        if(c[3][3] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if(c[4][4] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

        if(c[5][5] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satisfied.")

            
    if(crystaltype =="monoclinic"):
        print("Monoclinic crystal system.... \n")
        print("Born stability criteria for the stability of monoclinic systems are: [Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n")
        print("(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  ")
        print(" (vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 ")
        print(" (xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  ")
        print(" (xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  ")
        print("          where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] ")


	if(c[0][0] > 0.0):
            print("Condition (i) is satified.")
        else:
            print("Condition (i) is NOT satified.")


        if(c[1][1] > 0.0):
            print("Condition (ii) is satified.")
        else:
            print("Condition (ii) is NOT satified.")

        if(c[2][2] > 0.0):
            print("Condition (iii) is satified.")
        else:
            print("Condition (iii) is NOT satified.")

        if(c[3][3] > 0.0):
            print("Condition (iv) is satified.")
        else:
            print("Condition (iv) is NOT satified.")

        if(c[4][4] > 0.0):
            print("Condition (v) is satified.")
        else:
            print("Condition (v) is NOT satified.")

        if(c[5][5] > 0.0):
            print("Condition (vi) is satified.")
        else:
            print("Condition (vi) is NOT satified.")


        #for i in range(0, 6):
        #    if(c[i][i]  > 0.0):
        #        print("Condition (%2d) is satified." % (i+1))
        #    else:
        #        print("Condition (%2d) is NOT satified." % (i+1))


        if(c[0][0] + c[1][1] + c[2][2] + 2*(c[0][1] + c[0][2] + c[1][2]) > 0 ):
            print("Condition (vii) is satified.")
        else:
            print("Condition (vii) is NOT satisfied.")

        if (c[2][2]*c[4][4] - c[2][4]*c[2][4] > 0):
            print("Condition (viii) is satified.")
        else:
            print("Condition (viii) is NOT satisfied.")

        if(c[3][3]*c[5][5] - c[3][5]*c[3][5] > 0.0):
            print("Condition (ix) is satified.")
        else:
            print("Condition (ix) is NOT satisfied.")

        if(c[1][1] + c[2][2] - 2*c[1][2] > 0.0):
            print("Condition (x) is satified.")
        else:
            print("Condition (x) is NOT satisfied.")


	if(c[1][1]*(c[2][2]*c[4][4] - c[2][4]*c[2][4]) + 2*c[1][2]*c[1][4]*c[2][4] - c[1][2]*c[1][2]*c[4][4] - c[1][4]*c[1][4]*c[2][2] > 0.0):
            print("Condition (xi) is satified.")
        else:
            print("Condition (xi) is NOT satisfied.")

        g = (c[0][0]*c[1][1]*c[2][2]) - (c[0][0]*c[1][2]*c[1][2]) - (c[1][1]*c[0][2]*c[0][2]) - (c[2][2]*c[0][1]*c[0][1]) + 2.0*(c[0][1]*c[0][2]*c[1][2])
        
        h1 =  2*(c[0][4]*c[1][4]*(c[2][2]*c[0][1] - c[0][2]*c[1][2]) + c[0][4]*c[2][4]*(c[1][1]*c[0][2] - c[0][1]*c[1][2]) + c[1][4]*c[2][4]*(c[0][0]*c[1][2] - c[0][1]*c[0][2])) 
        
        h2 = (c[0][4]*c[0][4]*(c[1][1]*c[2][2] - c[1][2]*c[1][2]) + c[1][4]*c[1][4]*(c[0][0]*c[2][2] - c[0][2]*c[0][2]) + c[2][4]*c[2][4]*(c[0][0]*c[1][1] - c[0][1]*c[0][1]))
        
        x = h1 - h2 + c[4][4]*g

#  print 'x = %10.4f' % x
#  print 'g = %10.4f' % g

        if(x > 0.0):
            print("Condition (xii) is satified.")
        else:
            print("Condition (xii) is NOT satisfied.")

    return 




def elastic_const_2D(cnew):

    """Convert the units from GPa or N/m2 to N/m for two-dimensional systems
        1 GPa = 10^9 N/m2
        
        Here, we use second Piola-Kirchhoff stress method to express the 2D forces per unit length in N/m units.
        Ref: [Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013). ]  
        [Singh et al., Phys. Rev. B 95, 165444 (2017)]
        
        We multiply elastic tensor by the thickness of the simulation cell to consider the vacuum thickness.
        In 2D:  Cij = bulk_Cij * C_latticevector (final units N/m)

        For example: if bulk_Cij = 15 GPa and out-of-plane cell parameter c = 10 Angs.
                  Then  2D_Cij = [15*10^9 N/m2] * [10*10^(-10) m] ; i.e 15*(0.1*c) N/m """
     
    c2d=np.zeros((6,6))
    for i in range(0,6):
        for j in range(0,6):
            c2d[i][j]=cnew[i][j]*0.1*float(lattconst[2])

    print("\n \n Elastic tensor for two-dimensional system in N/m units \n")
    np.set_printoptions(precision=3, suppress=True)
    printMatrix(c2d)

    
    ## Elastic properties
    ## Layer modulus: represents the resistance of a 2D sheet to stretching; Lm = (1/4)*[c11 + c22 + 2*c12]  [Ref: Andrew et al.;  Phys. Rev. B 85, 125428 (2012)]
    Lm = 0.25*(c2d[0][0] + c2d[1][1] + 2*c2d[0][1] )

    ## 2D Young's modulus or in-plane stiffness: Y[10] = [c11c22 - c12^2]/[c22], Y[01] = [c11c22 - c12^2]/[c11]
    Y10 = (c2d[0][0]*c2d[1][1] - c2d[0][1]*c2d[0][1])/(c2d[1][1])
    Y01 = (c2d[0][0]*c2d[1][1] - c2d[0][1]*c2d[0][1])/(c2d[0][0])

    ## 2D Poisson's ratio;  nu10 = c12/c22, nu01 = c12/c11
    nu10 = c2d[0][1]/c2d[1][1]
    nu01 = c2d[0][1]/c2d[0][0]

    ## 2D shear modulus; G2d = C66
    G2d = c2d[5][5]

    print("\n \n Elastic properties in two-dimensions \n")
    print("[Useful refs. Andrew et al.;  Phys. Rev. B 85, 125428 (2012), Peng et al., Acta Mechanica 223 (2012), 2591-2596; Comput. Mater. Sci. 68, 320 (2013);  Mech. Mater. 64, 135 (2013)  ]")
    print("                   ")
    print("-------------------------------------------------------")
    print("2D layer modulus (N/m)         :   %10.3f "  % Lm)
    print("2D Young's modulus Y[10] (N/m) :   %10.3f " % Y10)
    print("2D Young's modulus Y[01] (N/m) :   %10.3f " % Y01)
    print("2D Shear modulus G (N/m)       :   %10.3f " % G2d)
    print("2D Poisson ratio v[10]         :   %10.3f " % nu10)
    print("2D Poisson ratio v[01]         :   %10.3f " % nu01)
    print("-------------------------------------------------------")
    print("Note:  The elastic stabilty test for 2D systems is not yet implemented. ")

    print_warning_2D()
    
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



## calculate elastic properties
def main():

    if (args.dim == '2D'):
        elastic_const_2D(cnew)

    elif (args.dim == '3D'):
        elastic_const_bulk(cnew,snew)

    else:
        elastic_const_bulk(cnew,snew)


if __name__ == '__main__':
    main()
    print("\n  Thanks! See you later. ")



