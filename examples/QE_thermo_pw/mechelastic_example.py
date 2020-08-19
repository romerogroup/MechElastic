# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 01:10:57 2020

@author: lllan
"""

from mechelastic.parsers import QE_thermo_pw_Parser
from mechelastic import calculate_elastic
from mechelastic.core import ElasticProperties

calculate_elastic(outfile = 'si.elastic.out' , infile = 'si.elastic.in', crystal = 'cubic', code = "qe_thermo_pw")


parser = QE_thermo_pw_Parser(outfile = 'si.elastic.out' , infile = 'si.elastic.in')

elastic_tensor = parser.elastic_tensor
structure = parser.structure
lattice_constant = parser.lattice_constant
crystal_type = 'cubic'

elastic_properties = ElasticProperties(elastic_tensor, structure, crystal_type)

bulk_modulus_voigt = elastic_properties.bulk_modulus_voigt