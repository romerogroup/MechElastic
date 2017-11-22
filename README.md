# MechElastic
This python scripts can be used to calculate some important physical properties such as elastic moduli, melting temperature, Debye temperature, elastic wave velocities, elastic anisotropy, etc. for all crystalline systems using the VASP output data from an elastic tensor calculation.

It can also be used to test the mechanical stability of any bulk system. 

The script reads the elastic matrix written in the OUTCAR file as input. 
 
NOTE: In order to evaluate accurate elastic properties and mechanical strength, one must well-converge the elastic constants by increasing the size of k-mesh and energy cutoff in the VASP calculation. 

An example input file (INCAR) for elastic constants calculations in VASP is given below: 

system   =  Si-227

PREC      =  High

NELMIN    =  8

NELM      =  100

EDIFF     =  1E-07

ISIF      =  3

LREAL     = .FALSE.

NSW       = 1

ISMEAR    = 0     

SIGMA     =  0.02

IBRION    = 6

POTIM     = 0.015

NFREE     = 2




********* TO RUN *********

Use the below line:

--PyElastic.py OUTCAR cubic

The script expects 'OUTCAR' file as the first argument and 'crystal type' as the second argument. If not provided by user, this script will read the crystal type from the OUTCAR file. Crystal type is needed only to perform the mechanical stability test. 


********* TO DO *********

Generalize the script to deal with two-dimensional systems.
