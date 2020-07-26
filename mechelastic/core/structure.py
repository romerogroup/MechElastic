# -*- coding: utf-8 -*-

import spglib
import numpy as np
from ..utils import elements
from ..utils.constants import N_avogadro


class Structure:
    def __init__(
        self, atoms, fractional_coordinates, lattice,
    ):
        self.fractional_coordinates = fractional_coordinates
        self.atoms = atoms
        self.lattice = lattice

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
        return np.unique(self.atoms)

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
