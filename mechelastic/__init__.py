# Copyright (C) 2020 Sobhit Singh
#
# This file is part of MechElastic.
#
# MechElastic is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MechElastic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MechElastic.  If not, see <http://www.gnu.org/licenses/>.

from .version import version as __version__
from .version import author as __author__
from .version import copyright as __copyright__
from .version import email as __email__
from .version import status as __status__
from .version import date as __date__


import re
import sys

# import matplotlib.pyplot as plt
import numpy as np

# import pyvista
import prettytable
import spglib

from mechelastic.parsers import VaspOutcar
from mechelastic.parsers import AbinitParser
from mechelastic.comms import printer
from mechelastic.core import Structure, ElasticProperties, ElasticProperties2D, ELATE
from mechelastic.tests import ductile, eigenvals, stability, symmetry
from mechelastic.utils import constants, elements, crystalutils

from .calculate_elastic_anisotropy import calculate_elastic_anisotropy
from .calculate_elastic import calculate_elastic
