# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 01:10:57 2020

@author: lllan
"""

from mechelastic.parsers import QE_thermo_pw_Parser
from mechelastic import calculate_elastic_anisotropy
from mechelastic.core import ELATE


"""This example is using qe_thermo_pw, but can easily be generailized to work with the other codes"""

"""Generate 2d plots with (plot = "2D") and select which elactic property you want (elastic_calc = 'SHEAR' , 'POISSON' , 'YOUNG , 'LC') """

#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'SHEAR')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'YOUNG')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'LC')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'POISSON')



"""Generate 3d plots with (plot = "3D") and select which elactic property you want (elastic_calc = 'SHEAR' , 'POISSON' , 'YOUNG , 'LC') 

Must have pyvista installed
"""

#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="3D", elastic_calc= 'SHEAR')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'YOUNG')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'LC')
#calculate_elastic_anisotropy(outfile = 'si.elastic.out', infile = 'si.elastic.in', code="qe_thermo_pw", plot="2D", elastic_calc= 'POISSON')

"""Call ELATE as an object and get various properties

The ELATE input is a list of row list. You can get them into the correct format with the following code
"""

output = QE_thermo_pw_Parser(outfile = 'si.elastic.out', infile = 'si.elastic.in' )
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