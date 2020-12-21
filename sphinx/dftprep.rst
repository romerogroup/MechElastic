.. _labeldftprep:

DFT Preparation
================

This section discusses steps to perform DFT calculations to obtain the elastic tensor reqired to run MechElastic to obtain elastic properties. Examples of these are available in the ``examples`` directory of the github repository. 

The flag ``code`` is to be set in MechElastic functions to select the DFT code.

E.g.

In stand-alone mode::

    MechElastic.py -d=3D -i OUTCAR-Si_bulk -co vasp

In library mode:: 

    import mechelastic 
    mechelastic.calculate_elastic(code="vasp", dim="3D", infile="OUTCAR-Si_bulk")


========
1. VASP
========

- Required files : OUTCAR 
- flag           : code='vasp' (default)


In order to evaluate accurate elastic properties and mechanical strength, one must well-converge the elastic constants by increasing the size of k-mesh and energy cutoff in the VASP calculation. 

An example input file (INCAR) for elastic constants calculations in VASP is given below:: 

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


=========
2. Abinit
=========

- Required files : abinit.out, abinit2.out
- flag           : code='abinit'

abinit.out - The output from a SCF calculation
abinit2.out - The output from the anaddb calculation to retrieve elastic tensors


===================
3. Quantum Espresso
===================

- Required files : scf.out
- flag           : code='qe'

Quantum Espresso v6.5+ is supported. 




