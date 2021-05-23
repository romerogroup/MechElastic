Stand-alone Mode
================

Please note that the stand-alone mode may be outdated compared to the library mode.

Arguments:

* '-d' - to provide the dimensionality of the system (2D or 3D)
* '-i' - to provide the name of the input file (output of DFT calculation)
* '-anaddb' - Output of the Abinit anaddb calculation.
* '-c' - used to provide information related to the crystal type. If crystal type is not provided by the user, MechElastic will read the crystal type from the output file. Crystal type is needed only to perform the mechanical stability test for bulk systems.
* '-co' - DFT code
* '-ap' - Flag to adjust pressure in Elastic Tensor (VASP). Default: 1 (True)
* '-v' - Enable verbosity.

Examples:

1. For bulk Si::

    MechElastic.py -d=3D -i OUTCAR-Si_bulk > output_Si_bulk.log

2. For 2D graphene::

    MechElastic.py -d=2D -i OUTCAR-graphene > output_graphene.log

3. For 2D BN::

    MechElastic.py -d=2D -i OUTCAR-BN_mono > output_BN_monolayer.log

4. Bulk Si Abinit calculation::

    MechElastic.py –i abinit.out –anaddb abinit2.out –co abinit
