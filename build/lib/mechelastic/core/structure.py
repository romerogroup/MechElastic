# -*- coding: utf-8 -*-

import spglib
import numpy as np
from ..utils import elements
from ..utils.constants import N_avogadro


class Structure:
    def __init__(
        self, atoms, fractional_coordinates, lattice,
    ):
        self.fractional_coordinates = np.array(fractional_coordinates)
        self.atoms = np.array(atoms)
        self.lattice = np.array(lattice)

    @property
    def volume(self):
        return abs(np.linalg.det(self.lattice)) * 1e-30

    @property
    def masses(self):
        return [elements.atomic_mass(x) * 1.0e-3 for x in self.atoms]

    @property
    def density(self):
        return np.sum(self.masses) / (self.volume * N_avogadro)

    @property
    def species(self):
        return np.unique(self.atoms).tolist()

    @property
    def nspecies(self):
        return len(self.species)

    @property
    def natoms(self):
        return len(self.atoms)

    @property
    def atomic_numbers(self):
        return [elements.atomic_number(x) for x in self.atoms]

    @property
    def spglib_cell(self):
        return (self.lattice, self.fractional_coordinates, self.atomic_numbers)

    def get_space_group_number(self, symprec=1e-5):
        return spglib.get_symmetry_dataset(self.spglib_cell, symprec)["number"]

    def get_space_group_international(self, symprec=1e-5):
        return spglib.get_symmetry_dataset(self.spglib_cell, symprec)["international"]

    def get_wyckoff_positions(self, symprec=1e-5):
        return spglib.get_symmetry_dataset(self.spglib_cell, symprec)["wyckoffs"]
    
    def to_dict(self, symprec=1e-5):
        return {'atomic_numbers':self.atomic_numbers,
                'atoms':self.atoms.tolist(),
                'density':self.density,
                'fractional_coordinates':self.fractional_coordinates.tolist(),
                'space_group_international':self.get_space_group_international(symprec),
                'space_group_number':self.get_space_group_number(symprec),
                'wyckoff_positions':self.get_wyckoff_positions(symprec),
                'lattice':self.lattice.tolist(),
                'masses':self.masses,
                'natoms':self.natoms,
                'nspecies':self.nspecies,
                'species':self.species,
                'volume':self.volume,
            }