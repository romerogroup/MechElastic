# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 01:12:27 2020

@author: lllan
"""

from mechelastic.parsers import QE_ElaStic_Parser
from mechelastic import calculate_elastic
from mechelastic.core import ElasticProperties

calculate_elastic(outfile = "ElaStic_2nd.out" ,infile = "scf.in", crystal = 'cubic', code = "qe_ElaStic")


parser = QE_ElaStic_Parser(outfile = "ElaStic_2nd.out" , infile = "scf.in")


elastic_tensor = parser.elastic_tensor
structure = parser.structure
lattice_constant = parser.lattice_constant
crystal_type = 'cubic'

elastic_properties = ElasticProperties(elastic_tensor, structure, crystal_type)

bulk_modulus_voigt = elastic_properties.bulk_modulus_voigt