# MechElastic
This Python scripts can be used to calculate some important physical properties such as elastic moduli, melting temperature, Debye temperature, elastic wave velocities, elastic anisotropy, etc. for all crystalline systems using the VASP output data from an elastic tensor calculation. It can also be used to test the mechanical stability of any bulk system. 

The script reads the elastic matrix written in the OUTCAR file as input. 

Please cite this paper if you use MechElastic script for your research: 

[Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.054108)


NOTE: In order to evaluate accurate elastic properties and mechanical strength, one must well-converge the elastic constants by increasing the size of k-mesh and energy cutoff in the VASP calculation. 

An example input file (INCAR) for elastic constants calculations in VASP is given below: 

```
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
```



## Installation

```
pip install mechelastic
```



## Running

**Note:** '-d' argument is used to provide the dimensionality of the system (2D or 3D), '-i' argument provides the name of the input OUTCAR_file, and '-c' can be used to provide information related to the crystal type. If crystal type is not provided by the user, this script will read the crystal type from the OUTCAR file. Crystal type is needed only to perform the mechanical stability test for bulk systems.  

Some examples are available in the examples directory.

For bulk Si:
```
MechElastic  -d=3D -i OUTCAR-Si_bulk > output_Si_bulk.log
```

For 2D graphene:

```
MechElastic  -d=2D -i OUTCAR-graphene > output_graphene.log
```

For 2D BN:

```
MechElastic  -d=2D -i OUTCAR-BN_mono > output_BN_monolayer.log
```





