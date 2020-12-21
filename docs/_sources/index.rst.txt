.. MechElastic documentation master file, created by
   sphinx-quickstart on Tue Aug 18 10:05:06 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MechElastic's documentation!
=======================================

The MechElastic Python package evaluates the mechanical and elastic properties of bulk and 2D materials using the elastic coefficient matrix (:math:`C_{ij}`) obtained from any ab-initio density-functional theory (DFT) code. The current version of this package reads the output of VASP, ABINIT and Quantum Espresso codes (but it can be easily generalized to any other DFT code) and performs the appropriate post-processing of elastic constants as per the requirement of the user. This program can also detect the input structure's crystal symmetry and test the mechanical stability of all crystal classes using the Born-Huang criteria. Various useful material-specific properties such as elastic moduli, longitudinal and transverse elastic wave velocities, Debye temperature, elastic anisotropy, 2D layer modulus, hardness, Pugh's ratio, Cauchy's pressure, Kleinman's parameter, and Lame's coefficients, can be estimated using this program. Another existing feature of this program is to employ the `ELATE package <https://iopscience.iop.org/article/10.1088/0953-8984/28/27/275201/meta>`_ and plot the spatial variation of several elastic properties such as Poisson's ratio, linear compressibility, shear modulus, and Young's modulus in three dimensions. Further, the MechElastic package can plot the equation of state (EOS) curves for energy and pressure for a variety of EOS models such as Murnaghan, Birch, Birch-Murnaghan, and Vinet, by reading the inputted energy/pressure vs. volume data obtained via numerical calculations or experiments. This package is particularly useful for the high-throughput analysis of elastic and mechanical properties of materials.

Currently supports:

1. VASP
2. Abinit 
3. Quantum Espresso


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   contributors
   cite
   dftprep
   tutorials

   modules

	     

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
