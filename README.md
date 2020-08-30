[![PyPI version](https://badge.fury.io/py/MechElastic.svg)](https://badge.fury.io/py/MechElastic)
[![HitCount](http://hits.dwyl.com/uthpalaherath/romerogroup/mechelastic.svg)](http://hits.dwyl.com/uthpalaherath/romerogroup/mechelastic)
![PyPI - Downloads](https://img.shields.io/pypi/dm/mechelastic)

# MechElastic
This Python library can be used to calculate some important physical properties such as elastic moduli, melting temperature, Debye temperature, elastic wave velocities, elastic anisotropy, etc. for all crystalline systems using the VASP output data from an elastic tensor calculation. It can also be used to test the mechanical stability of any bulk system. 

MechElastic reads the elastic matrix written in the OUTCAR file as input. 

Please cite this paper if you use MechElastic for your research: 

[Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.054108)

BibTeX:

```
@article{PhysRevB.97.054108,
  title = {Elastic, mechanical, and thermodynamic properties of Bi-Sb binaries: Effect of spin-orbit coupling},
  author = {Singh, Sobhit and Valencia-Jaime, Irais and Pavlic, Olivia and Romero, Aldo H.},
  journal = {Phys. Rev. B},
  volume = {97},
  issue = {5},
  pages = {054108},
  numpages = {11},
  year = {2018},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.97.054108},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.97.054108}
}
```

**NOTE:** In order to evaluate accurate elastic properties and mechanical strength, one must well-converge the elastic constants by increasing the size of k-mesh and energy cutoff in the VASP calculation. 

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
pip install MechElastic
```

\* Packages in PyPI are not case sensitive. 

Once installed, use the ``-h`` flag to see a list of options for the stand-alone mode.

```
MechElastic -h
```

## Usage

MechElastic has two modes: stand-alone and library. </br>

### Stand-alone Mode

**Note:** '-d' argument is used to provide the dimensionality of the system (2D or 3D), '-i' argument provides the name of the input OUTCAR_file, and '-c' can be used to provide information related to the crystal type. If crystal type is not provided by the user, MechElastic will read the crystal type from the OUTCAR file. Crystal type is needed only to perform the mechanical stability test for bulk systems.  

Some examples are available in the **examples** directory.

For bulk Si:
```
MechElastic.py -d=3D -i OUTCAR-Si_bulk > output_Si_bulk.log
```

For 2D graphene:

```
MechElastic.py -d=2D -i OUTCAR-graphene > output_graphene.log
```

For 2D BN:

```
MechElastic.py -d=2D -i OUTCAR-BN_mono > output_BN_monolayer.log
```

MechElastic currently supports VASP and Abinit.
To run the calculation for Abinit, two inputs are required.

1. The output from a SCF calculation (abinit.out)
2. The output from the anaddb calculation to retrieve elastic tensors (abinit2.out)

The following command will run MechElastic for an Abinit calculation:

```
MechElastic.py –i abinit.out –ddb abinit2.out –co abinit
```

### Library Mode

Similarly, the above examples can be performed with the library mode after importing MechElastic from Python.

For bulk Si:
```
import mechelastic 
mechelastic.calculate_elastic(code="vasp", dim="3D", infile="OUTCAR-Si_bulk")
```

For 2D graphene:

```
import mechelastic 
mechelastic.calculate_elastic(code="vasp", dim="2D", infile="OUTCAR-graphene")
```

For 2D BN:

```
import mechelastic 
mechelastic.calculate_elastic(code="vasp", dim="2D", infile="OUTCAR-BN_mono")
```

To provide a crystal type (required for the stability test) manualy, the ``crystal`` flag can be used. If not provided, MechElatic will determind the crystal symmetry using spglib.
The stability test is currently only required for 3D systems.

```
import mechelastic 
mechelastic.calculate_elastic(code="vasp", dim="3D", infile="OUTCAR-Si_bulk", crystal="cubic")
```

To run for Abinit:

```
mechelastic.calculate_elastic(code="abinit", infile="abinit.out", anaddbfile="abinit2.out")
```


``mechelastic.calculate_elastic()`` calculates the complete set of elastic properties. However, if one wishes to call particular methods that can be also done through the library mode. For example, given a matrix and a crystaltype, the stability can be determined:

```
import mechelastic

parserclass = mechelastic.parsers.VaspOutcar()
elastic_tensor = parserclass.elastic_tensor
crystaltype = "cubic"
mechelastic.tests.stability.stability_test(elastic_tensor, crystaltype)
```

To determine the crystal symmetry:

```
import mechelastic

parserclass = mechelastic.parsers.VaspOutcar()
elastic_tensor = parserclass.elastic_tensor
cell = parserclass.structure.spglib_cell
mechelastic.utils.crystalutils.crystal_select(elastic_tensor, cell)
```







