Elastic constants and base functions 
====================================

1. For bulk Si::

    import mechelastic 
    mechelastic.calculate_elastic(code="vasp", dim="3D", infile="OUTCAR-Si_bulk")

2. For 2D graphene::

    import mechelastic 
    mechelastic.calculate_elastic(code="vasp", dim="2D", infile="OUTCAR-graphene")

3. For 2D BN::
   
    import mechelastic 
    mechelastic.calculate_elastic(code="vasp", dim="2D", infile="OUTCAR-BN_mono")

4. To provide a crystal type (required for the stability test) manualy, the ``crystal`` flag can be used. If not provided, MechElatic will determind the crystal symmetry using spglib. The stability test is currently only required for 3D systems::

    import mechelastic 
    mechelastic.calculate_elastic(code="vasp", dim="3D", infile="OUTCAR-Si_bulk", crystal="cubic")

5. To run elastic constants calculation for Abinit::
   
    import mechelastic 
    mechelastic.calculate_elastic(code="abinit", infile="abinit.out", anaddbfile="abinit2.out")


6. ``mechelastic.calculate_elastic()`` calculates the complete set of elastic properties. However, if one wishes to call particular methods that can be also done through the library mode. For example, given a matrix and a crystaltype, the stability can be determined::

    import mechelastic

    parserclass = mechelastic.parsers.VaspOutcar()
    elastic_tensor = parserclass.elastic_tensor
    crystaltype = "cubic"
    mechelastic.tests.stability.stability_test(elastic_tensor, crystaltype)

7. To determine the crystal symmetry::

    import mechelastic

    parserclass = mechelastic.parsers.VaspOutcar()
    elastic_tensor = parserclass.elastic_tensor
    cell = parserclass.structure.spglib_cell
    mechelastic.utils.crystalutils.crystal_select(elastic_tensor, cell)