Quick ELATE example 
====================

With the output of a given DFT code or direct input of an elastic tensor MechElastic can perform ELATE's anisotropic analysis and produce the same 2d and 3d plots of elastic properties. The original website for ELATE can be found here, <http://progs.coudert.name/elate>.

Examples of the implementation can be found in \MechElastic\examples\ELATE.py or \MechElastic\examples\mechelastic_w_mpDatabase.py



Direct output from a DFT code usage::

    import mechelastic
    # Produce the 2D plots for the Shear modulus 
    # elastic_calc - ['SHEAR','POISSON','YOUNG','LC'] 
    # npoints - number of points to use to make the plot
    mechelastic.calculate_elastic_anisotropy(outfile = 'si.elastic.out',infile = 'si.elastic.in',                                                                            code="qe_thermo_pw",                                                                                                             plot="2D", 
                                             elastic_calc= 'SHEAR'
                                             npoints = 100)
    # Produce the 3D plots for the Shear modulus
    # elastic_calc - ['SHEAR','POISSON','YOUNG','LC'] 
    # npoints - number of points to use to make the plot
    mechelastic.calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', 
                                             code="qe_thermo_pw",                                                    
                                             plot="3D", 
                                             elastic_calc= 'SHEAR'
                                             npoints = 100)
    # Just output the ELATE analysis
    mechelastic.calculate_elastic_anisotropy(outfile = 'si.elastic.out',infile = 'si.elastic.in', code="qe_thermo_pw")



Direct input of a elastic tensor usage::

    from mechelastic.core import ELATE
    from mechelastic.parsers import QE_thermo_pw_Parser
    
    output = QE_thermo_pw_Parser(outfile = 'si.elastic.out', 
                                 infile =   'si.elastic.in' )
    # This is a 6X6 matrix. Input your tensor here.
    elastic_tensor = output.elastic_tensor
    
    row = elastic_tensor.shape[0]
    col = elastic_tensor.shape[1]
    rowsList = []
    for i in range(row):
        columnsList = []
        for j in range(col):
            columnsList.append(round(elastic_tensor[i, j],3))
        rowsList.append(columnsList)
    elastic_tensor = ELATE.ELATE(rowsList)
    
    voigt_shear  = elastic_tensor.voigtShear
    
    elastic_tensor.plot3D(elastic_calc="LC", npoints = 100)
    elastic_tensor.plot2D(elastic_calc="LC", npoints = 100)


.. toctree::
   :maxdepth: 4 

   elate_example
   plotting_2D