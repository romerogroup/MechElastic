[![PyPI version](https://badge.fury.io/py/MechElastic.svg)](https://badge.fury.io/py/MechElastic)
[![HitCount](http://hits.dwyl.com/uthpalaherath/romerogroup/mechelastic.svg)](http://hits.dwyl.com/uthpalaherath/romerogroup/mechelastic)
![PyPI - Downloads](https://img.shields.io/pypi/dm/mechelastic)

# MechElastic

This Python library can be used to calculate some important physical properties such as elastic moduli, melting temperature, Debye temperature, elastic wave velocities, elastic anisotropy, etc. for all crystalline systems using output data from an elastic tensor calculation. It can also be used to test the mechanical stability of any bulk system. MechElastic reads the elastic matrix written in the output files from DFT codes. Additionally, MechElastic allows performing Equation of State (EOS) analysis provided inputs of Volume and Energy or Volume and Pressure.  

Currently supports:

- VASP
- Abinit 
- Quantum Espresso


## Documentation

https://romerogroup.github.io/MechElastic/


Developers
------------
Sobhit Singh <br />
Logan Lang <br />
Viviana Dovale-Farelo <br />
Uthpala Herath <br />
Pedram Tavadze <br />
François-Xavier Coudert <br />
Aldo H. Romero <br />

Mailing list
-------------
Please post your questions on our forum.

https://groups.google.com/d/forum/mechelastic

## Installation

```
pip install mechelastic
```

## Usage

Once installed, use the ``-h`` flag to see a list of options for the stand-alone mode.

```
MechElastic.py -h
```

For example: 

```
MechElastic.py -d=3D -i OUTCAR-Si_bulk 
```

For more information please refer to the documentation. 

Citations
------------------

Please cite the following articles if you use MechElastic for your research: 

- [Sobhit Singh, Irais Valencia-Jaime, Olivia Pavlic, and Aldo H. Romero; Phys. Rev. B 97, 054108 (2018).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.054108)

- [Sobhit Singh, Logan Lang, Viviana Dovale-Farelo, Uthpala Herath, Pedram Tavadze, François-Xavier Coudert, and Aldo H. Romero; Computer Physics Communications 267, 108068 (2021).](https://doi.org/10.1016/j.cpc.2021.108068)

BibTeX:

```
@article{singh2018mechelastic,
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

@article{singh2021mechelastic,
	Author = {Sobhit Singh and Logan Lang and Viviana Dovale-Farelo and Uthpala Herath and Pedram Tavadze and Fran{\c c}ois-Xavier Coudert and Aldo H. Romero},
	Doi = {https://doi.org/10.1016/j.cpc.2021.108068},
	Issn = {0010-4655},
	Journal = {Computer Physics Communications},
	Keywords = {Elastic properties, Mechanical stability, Elastic wave velocities, Elastic anisotropy, 2D materials, High-throughput, DFT, Equation of state},
	Pages = {108068},
	Title = {MechElastic: A Python library for analysis of mechanical and elastic properties of bulk and 2D materials},
	Url = {https://www.sciencedirect.com/science/article/pii/S0010465521001806},
	Volume = {267},
	Year = {2021},
	Bdsk-Url-1 = {https://www.sciencedirect.com/science/article/pii/S0010465521001806},
	Bdsk-Url-2 = {https://doi.org/10.1016/j.cpc.2021.108068}}
}

```









