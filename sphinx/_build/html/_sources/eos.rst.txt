Equation of State (EOS) Analysis 
================================

Provided pressure/energy and volume data, MechElastic can perform EOS analysis using the following models:

1. Vinet
2. Birch
3. Murnaghan
4. Birch-Murnaghan

Please find examples for EOS analysis in the ``examples/EOS/`` directory.

plot_eos()
----------

**INPUT:**
Text file with two columns containing volume (:math:`Å^3`) and energy (eV) (``eostype='energy'``) or volume and pressure (GPa) (``eostype='pressure'``), respectively, should be passed into the MechElastic package. 


With this, an initial second order parabolic polyfit is performed to obtain the initial fitting parameters, :math:`E_0` or :math:`P_0`, :math:`B_0`, :math:`B_p` and :math:`V_0`. Afterwards, a more accurate fitting for the energy or pressure is performed using the least square method against each EOS model. A central difference scheme is used to obtain the pressure from this fitted energy. If pressure is to be obtained from energy, an integrating scheme is used instead. Finally, a plot of `Energy vs Volume' and `Pressure vs Volume' is obtained comparing the values obtained from the different EOS models with the originally provided raw data. With these plots, an analysis of the phase boundaries could be performed. The energy or pressure is normalized with the number of atoms in the unit-cell set by ``natoms``. If the output units are to be in atomic units, set ``au=True``.  

In addition to the fitting, MechElastic also performs a regression analysis for fitting using the Mean-Squared Error (MSE) of the residuals. This allows users to decide the best EOS model to use for their dataset. 

Usage::

    from mechelastic import EOS
    eos_object = EOS()
    eos_object.plot_eos(inputfile='EvsV.dat', eostype='energy', natoms=1, au=False)

To set a initial and final value to the volume range use the flag ``vlim=[v_initial, v_final]``. If not set, the minimum and maximum values of the provided dataset would be used.  

plot_enthalpy_curves()
----------------------

This function plots Enthalpy vs. Pressure curves for multiple datasets of different phases and calculates the phase transition pressures. The pressure is calculated through a central difference scheme similar to ``plot_eos()``. 

**INPUT:**
Multiple input files containing two columns of volume (:math:`Å^3`) and energy (eV) for each phase.

Usage::

    from mechelastic import EOS

    eos_object = EOS()

    infiles = ["Rm-3m", "Pm-3m", "Im-3m", "I4_mcm", "C12:m1"]
    natoms = [6, 1, 2, 8, 4]
    eos_object.plot_enthalpy_curves(infiles, natoms, au=True)

The flags ``vlim`` and ``au`` work similar to the ``plot_eos()`` section. 


