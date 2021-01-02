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



1 . **QE ElaStic** implementation

- Required files : scf.in , ElaStic_2nd.out

- flag           : code='qe_ElaStic'

scf.in - The input file for a SCF calculation from the ElaStic code.

ElaStic_2nd.out - The output file from the ElaStic code.

See \MechElastic\examples\QE_ElaStic for an example of its use.



2 . **QE_thermo_pw** implementation

- Required files : si.elastic.in , si.elastic.out
- flag           : code='qe_thermo_pw'

si.elastic.in - The input file for a SCF calculation for the qe_thermo_pw code.

si.elastic.out - The output file for a SCF calculation for the qe_thermo_pw code.

See \MechElastic\examples\QE_thermo_pw for an example of its use.

Quantum Espresso v6.5+ is supported. 




