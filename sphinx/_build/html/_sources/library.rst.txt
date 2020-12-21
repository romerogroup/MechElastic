Library Mode
============

Similar to the Stand-alone mode, MechElastic can be run with the library mode after importing MechElastic from Python.

Elastic constants and base functions 
------------------------------------

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

Equation of State (EOS) Analysis 
--------------------------------

Provided pressure/energy and volume data, MechElastic can perform EOS analysis using the following models:

1. Vinet
2. Birch
3. Murnaghan
4. Birch-Murnaghan

**INPUT:**
Text file with two columns containing volume and energy (``eostype='energy'``) or volume and pressure (``eostype='pressure'``), respectively, should be passed into the MechElastic package. 


With this, an initial second order parabolic polyfit is performed to obtain the initial fitting parameters, :math:`E_0` or :math:`P_0`, :math:`B_0`, :math:`B_p` and :math:`V_0`. Afterwards, a more accurate fitting for the energy or pressure is performed using the least square method against each EOS model. A central difference scheme is used to obtain the pressure from this fitted energy. If pressure is to be obtained from energy, an integrating scheme is used instead. Finally, a plot of `Energy vs Volume' and `Pressure vs Volume' is obtained comparing the values obtained from the different EOS models with the originally provided raw data. With these plots, an analysis of the phase boundaries could be performed. 

In addition to the fitting, MechElastic also performs a regression analysis for fitting using the Mean-Squared Error (MSE) of the residuals. This allows users to decide the best EOS model to use for their dataset. 

Usage::

    from mechelastic import EOS
    eos_object = EOS('energyvsvolume.dat', eostype='energy')
    eos_object.plot_eos()


