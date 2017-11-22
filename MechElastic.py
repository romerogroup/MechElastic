#!/usr/bin/env python

''' This script reads elastic tensor from OUTCAR, and calculates elastic and mechanical properties... 
 Authors: Sobhit Singh and Aldo Romero
	  West Virginia University, Morgantown, USA

 To RUN:
 >> python MechElastic.py OUTCAR cubic '''

import numpy as np
import sys
from numpy import linalg as LA

if len(sys.argv) < 2:
  print "OUTCAR file is expected as the first argument \n"
else:
  print "reading OUTCAR ....."


## Defining useful constants
## Avogadro number
N_avogadro = 6.022140857e+23
## Planck's constant
h_Planck = 6.62607004e-34
## Boltzmann constant
kB = 1.38064852e-23

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

s=np.zeros((6,6))
c=np.zeros((6,6))

mass=[]
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
        print i.split()
    if "ions per type =" in i:
        natoms=map(int,i.split()[4:])
    if "LATTYP" in i:
        lattyp=i
    #   break

mass=np.array(mass)
print "Mass of atoms (in g/mol units): "
print(np.array(mass))
print "Number of atoms: %d" %np.sum(natoms)

totalmass=0.0
for i in range(len(natoms)):
        totalmass = totalmass + natoms[i]*mass[i]
print "Total mass (in g/mol): %10.4f " %totalmass
print "Volume of the cell (in Ang^3 units): %10.4f " %vol

## converting the units
vol = vol*1.0e-30     ## from Angstrom to meters
totalmass = totalmass*1.0e-3 ## from gram to kg
density = totalmass/(vol*N_avogadro)

print '\n Density (in kg/m^3 units ): %10.5f' % density


for i in range(0,6):
    l=lines[ii+3+i]
    word = l.split()
    s[i][:]=word[1:7]
    c[i]=map(float,s[i][:])

#print c

mc=np.matrix(c)
mci=mc.I

for i in range(0,6):
    for j in range(0,6):
        s[i][j]=mci[i,j]

print "Printing Cij matrix as read from OUTCAR\n"
print(np.matrix(c))

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

print "\n \n \n printing CNEW: Modified matrix in correct order (in GPa units)....\n"
#np.set_printoptions(precision=3)
np.set_printoptions(precision=3, suppress=True)
print(np.matrix(cnew))

def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)

print('\n Checking if the modified matrix CNEW is symmetric: i.e. Cij = Cji:  %10s' % check_symmetric(cnew))


print "\nEigen Values of the matrix CNEW:"
evals = LA.eigvals(cnew)
print(evals)


check = 0
for i in range(len(evals)):
	if evals[i] > 0.0:
		pass
	else:
		check = 1

#if np.all(evals) > 0.0:
#	print "All eigen values are positive indicating elastic stability."
#else:
#	print "ATTENTION: One or more eigen values are negative indicating elastic instability."


if check == 1:
	print "ATTENTION: One or more eigen values are negative indicating elastic instability."

if check == 0:
	print "All eigen values are positive indicating elastic stability."


## calculate elastic properties

def main():
	elastic_const(cnew,snew)
	
	crystal = np.array(['cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'])
	if len(sys.argv) == 3:
        	if (np.sum(crystal == sys.argv[2]) > 0):
	          print "Given crystal system:  %10s"  % sys.argv[2]
        	  crystaltype = sys.argv[2]
	          print "\n \t \t Mechanical Stability Test \n"
       	  	  stability_test(cnew, crystaltype)
	        else:
        	  print "\t \t WARNING: crystal symmetry class  was not provided by user, it will be taken from the OUTCAR"
		  print "\t \t One of the following was expected as the second argument: 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'"
                  crystaltype=''
                  for x in crystal:
                      if x in lattyp:
                          crystaltype=x
                  if crystaltype=='':
                      crystaltype='monoclinic'
                  print 'From OUTCAR the crystal type is = ',crystaltype
       	  	  stability_test(cnew, crystaltype)
	          #print "\t \t To perform the mechanical stability test you need to enter the crystal symmetry class: one of the following is expected as the second argument."
        	  #print "\t \t 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic' "
	else:
       	  print "\t \t WARNING: crystal symmetry class  was not provided by user, it will be taken from the OUTCAR"
          print "\t \t One of the following was expected as the second argument: 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic'"
          crystaltype=''
          for x in crystal:
                if x in lattyp:
                  crystaltype=x
                if crystaltype=='':
                  crystaltype='monoclinic'
          print 'From OUTCAR the crystal type is = ',crystaltype
       	  stability_test(cnew, crystaltype)
	  #print "\n \t \t WARNING: the mechanical stability test won't be performed. Crystal type is expected as the second argument in order to perform the stability test. \n"
	  #print "\t \t Choose one of the following: 'cubic', 'hexagonal', 'tetragonal', 'rhombohedral-1', 'rhombohedral-2', 'orthorhombic', 'monoclinic' \n"


def elastic_const(cnew, snew):

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
	

	print "\n \n                         Voigt     Reuss    Average"
	print "-------------------------------------------------------"
	print "Bulk modulus   (GPa)  %9.3f %9.3f %9.3f " % (KV, KR, Kvrh)
	print "Shear modulus  (GPa)  %9.3f %9.3f %9.3f " % (GV, GR, Gvrh)
	print "Young modulus  (GPa)  %9.3f %9.3f %9.3f " % (EV, ER, Evrh)
	print "Poisson ratio         %9.3f %9.3f %9.3f " % (Nu_V, Nu_R, Nu_vrh)
	print "P-wave modulus  (GPa) %9.3f %9.3f %9.3f " % (MV, MR, Mvrh)
	print "Bulk/Shear ratio      %9.3f %9.3f %9.3f (%s) " % (KG_ratio_V, KG_ratio_R, KG_ratio_vrh,  ductile_test(KG_ratio_vrh) )
	print "-------------------------------------------------------"


	print " \n \n \t \t Elastic Anisotropy \n "
	print "Zener anisotropy (true for cubic crystals only) Az = %10.5f" %Az
	print "Universal anisotropy index (Ranganathan and Ostoja-Starzewski method; PRL 101, 055504 (2008)) Au = %10.5f" %AU
	print "Log-Euclidean anisotropy parameter by Christopher M. Kube, AIP Advances 6, 095209 (2016) AL = %10.5f " %AL


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

	## melting temperature using empirical relation from Ref: Johnston I, Keeler G, Rollins R and Spicklemire S 1996 
	##							  Solid State Physics Simulations, The Consortium for Upper-Level Physics Software (New York: Wiley)
	Tm = 607 + 9.3*Kvrh

        print "\n \t \t  Elastic wave velocities calculated using Navier's equation  (in m/s units) \n"
        print "----------------------------------------------- "
        print "Longitudinal wave velocity (vl) : %10.5f " %vl
        print "Transverse wave velocity (vt) : %10.5f " %vt
        print "Average wave velocity (vm) : %10.5f " %vm
        print "Debye temperature  (in K) : %10.5f " %theta
        print "\nMelting temperature calculated from empirical relation: Tm = 607 + 9.3*Kvrh \pm 555 (in K)"
	print "Tm (in K)=  %10.5f (plus-minus 555 K) " % Tm
        print "----------------------------------------------- "


	

def ductile_test(ratio):
	if(ratio > 1.75):
		return "ductile"
	else:
		return "brittle"


def stability_test(matrix, crystaltype):
 c = np.copy(matrix)

 if(crystaltype =="cubic"):
   print "Cubic crystal system \n"
   print "Born stability criteria for the stability of cubic system are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 - C12 > 0;    (ii) C11 + 2C12 > 0;   (iii) C44 > 0 \n "

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print "Condition (i) satified."
   else:
      print "Condition (i) NOT satisfied."

   if(c[0][0] + 2*c[0][1] > 0.0):
      print "Condition (ii) satified."
   else:
      print "Condition (ii) NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iii) satified."
   else:
      print "Condition (iii) NOT satisfied."


 if(crystaltype =="hexagonal"):
   print "Hexagonal crystal system \n"
   print "Born stability criteria for the stability of hexagonal system are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0 \n "

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print "Condition (i) satified."
   else:
      print "Condition (i) NOT satisfied."

   if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
      print "Condition (ii) satified."
   else:
      print "Condition (ii) NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iii) satified."
   else:
      print "Condition (iii) NOT satisfied."


 if(crystaltype =="tetragonal"):
   print "Tetragonal crystal system \n"
   print "Born stability criteria for the stability of Tetragonal system are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 - C12 > 0;    (ii) 2*C13^2 < C33(C11 + C12);   (iii) C44 > 0;   (iv) C66 > 0;    (v) 2C16^2 < C66*(C11-C12) \n "

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print "Condition (i) is satified."
   else:
      print "Condition (i) is NOT satisfied."

   if(2*(c[0][2]*c[0][2]) < c[2][2]*(c[0][0] + c[0][1])):
      print "Condition (ii) is satified."
   else:
      print "Condition (ii) is NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iii) is satified."
   else:
      print "Condition (iii) is NOT satisfied."

   if(c[5][5] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."

   if(2*c[0][5]*c[0][5] < c[5][5]*(c[0][0] - c[0][1])):
      print "Condition (v) is satified."
   else:
      print "Condition (v) is NOT satisfied."


 if(crystaltype =="rhombohedral-1"):
   print "Rhombohedral (class-1): i.e. structures with point group: 3m, -3m and 32 \n"
   print "Born stability criteria for the stability of Rhombohedral-1 class system are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0; \n "

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print "Condition (i) is satified."
   else:
      print "Condition (i) is NOT satisfied."

   if((c[0][2]*c[0][2]) < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
      print "Condition (ii) is satified."
   else:
      print "Condition (ii) is NOT satisfied."

   if(c[0][3]*c[0][3] < 0.5*c[3][3]*(c[0][0] - c[0][1])):
      print "Condition (iii) is satified."
   else:
      print "Condition (iii) is NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."


 if(crystaltype =="rhombohedral-2"):
   print "Rhombohedral (class-2): i.e structures with point group: 3, -3 \n"
   print "Born stability criteria for the stability of Rhombohedral-1 class system are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 - C12 > 0;    (ii) C13^2 < (1/2)*C33(C11 + C12);   (iii) C14^2 + C15^2 < (1/2)*C44*(C11-C12) = C44*C66;   (iv)  C44 > 0;  Note: C15 is added.. \n "

   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0] - c[0][1] > 0.0):
      print "Condition (i) is satified."
   else:
      print "Condition (i) is NOT satisfied."

   if(c[0][2]*c[0][2] < (0.5)*c[2][2]*(c[0][0] + c[0][1])):
      print "Condition (ii) is satified."
   else:
      print "Condition (ii) is NOT satisfied."

   if(c[0][3]*c[0][3] + c[0][4]*c[0][4]  < 0.5*c[3][3]*(c[0][0] - c[0][1])):
      print "Condition (iii) is satified."
   else:
      print "Condition (iii) is NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."


 if(crystaltype =="orthorhombic"):
   print "Orthorhombic crystal system.... \n"
   print "Born stability criteria for the stability of Orthorhombic systems are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014)]  \n"
   print "(i) C11 > 0;   (ii) C11*C22 > C12^2;   (iii) C11*C22*C33 + 2C12*C13*C23 - C11*C23^2 - C22*C13^2 - C33*C12^2 > 0;   (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0 \n"
   ## check (i)   keep in mind list starts with 0, so c11 is stored as c00
   if(c[0][0]  > 0.0):
      print "Condition (i) is satified."
   else:
      print "Condition (i) is NOT satisfied."

   if(c[0][0]*c[1][1] < c[0][1]*(c[0][1])):
      print "Condition (ii) is satified."
   else:
      print "Condition (ii) is NOT satisfied."

   if(c[0][0]*c[1][1]*c[2][2] + 2*c[0][1]*c[0][2]*c[1][2] - c[0][0]*c[1][2]*c[1][2] - c[1][1]*c[0][2]*c[0][2] - c[2][2]*c[0][1]*c[0][1] > 0 ):
      print "Condition (iii) is satified."
   else:
      print "Condition (iii) is NOT satisfied."

   if(c[3][3] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."

   if(c[4][4] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."

   if(c[5][5] > 0.0):
      print "Condition (iv) is satified."
   else:
      print "Condition (iv) is NOT satisfied."


 if(crystaltype =="monoclinic"):
   print "Monoclinic crystal system.... \n"
   print "Born stability criteria for the stability of monoclinic systems are : [Ref- Mouhat and Coudert, PRB 90, 224104 (2014), and Wu et al. PRB 76, 054115 (2007)]  \n"
   print "(i) C11 > 0;  (ii)  C22 > 0; (iii)  C33 > 0; (iv)  C44 > 0;   (v)  C55 > 0 ;   (vi)  C66 > 0  "
   print " (vii) [C11 + C22 + C33 + 2*(C12 + C13 + C23)] > 0;    (viii)  C33*C55 - C35^2 > 0;   (ix)  C44*C66 - C46^2 > 0;   (x) C22 + C33 - 2*C23  > 0 "
   print " (xi) C22*(C33*C55 - C35^2) + 2*C23*C25*C35 - (C23^2)*C55 - (C25^2)*C33   > 0  "
   print " (xii)  2*[C15*C25*(C33*C12 - C13*C23) + C15*C35*(C22*C13 - C12*C23) + C25*C35*(C11*C23 - C12*C13)] - [C15*C15*(C22*C33 - C23^2) + C25*C25*(C11*C33 - C13^2) + C35*C35*(C11*C22 - C12^2)] + C55*g > 0  "
   print "          where, g = [C11*C22*C33 - C11*C23*C23 - C22*C13*C13 - C33*C12*C12 + 2*C12*C13*C23 ] "

   for i in range(0, 6):
      if(c[i][i]  > 0.0):
        print "Condition (%2d) is satified." % (i+1)
      else:
        print "Condition (%2d) is NOT satified." % (i+1)


   if(c[0][0] + c[1][1] + c[2][2] + 2*(c[0][1] + c[0][2] + c[1][2]) > 0 ):
      print "Condition (vii) is satified."
   else:
      print "Condition (vii) is NOT satisfied."

   if (c[2][2]*c[4][4] - c[2][4]*c[2][4] > 0):
      print "Condition (viii) is satified."
   else:
      print "Condition (viii) is NOT satisfied."


   if(c[3][3]*c[5][5] - c[3][5]*c[3][5] > 0.0):
      print "Condition (ix) is satified."
   else:
      print "Condition (ix) is NOT satisfied."


   if(c[1][1] + c[2][2] - 2*c[1][2] > 0.0):
      print "Condition (x) is satified."
   else:
      print "Condition (x) is NOT satisfied."


   if(c[1][1]*(c[2][2]*c[2][4] - c[2][4]*c[2][4]) + 2*c[1][2]*c[1][4]*c[2][4] - c[1][4]*c[1][4]*c[2][2] > 0.0):
      print "Condition (xi) is satified."
   else:
      print "Condition (xi) is NOT satisfied."

   g = (c[0][0]*c[1][1]*c[2][2]) - (c[0][0]*c[1][2]*c[1][2]) - (c[1][1]*c[0][2]*c[0][2]) - (c[2][2]*c[0][1]*c[0][1]) + 2.0*(c[0][1]*c[0][2]*c[1][2])
   h1 =  2*(c[0][4]*c[1][4]*(c[2][2]*c[0][1] - c[0][2]*c[1][2]) + c[0][4]*c[2][4]*(c[1][1]*c[0][2] - c[0][1]*c[1][2]) + c[1][4]*c[2][4]*(c[0][0]*c[1][2] - c[0][1]*c[0][2])) 
   h2 = (c[0][4]*c[0][4]*(c[1][1]*c[2][2] - c[1][2]*c[1][2]) + c[1][4]*c[1][4]*(c[0][0]*c[2][2] - c[0][2]*c[0][2]) + c[2][4]*c[2][4]*(c[0][0]*c[1][1] - c[0][1]*c[0][1]))
   x = h1 - h2 + c[4][4]*g

#  print 'x = %10.4f' % x
#  print 'g = %10.4f' % g

   if(x > 0.0):
      print "Condition (xii) is satified."
   else:
      print "Condition (xii) is NOT satisfied."



if __name__ == '__main__':
    main()

