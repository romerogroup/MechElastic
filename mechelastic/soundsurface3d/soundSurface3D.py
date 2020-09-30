# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:02:51 2020

@author: lllan
"""

import numpy as np
import scipy.interpolate as interpolate
from ..core import Isosurface, BrillouinZone


class SoundSurface3D(Isosurface):
    def __init__(self,
                 qpoints=None,
                 soundSpeedArray=None,
                 soundSpeed = None,
                 reciprocal_lattice=None,
                 supercell=[1, 1, 1]):
        """
        

        Parameters
        ----------
        kpoints : (n,3) float
            A list of kpoints used in the DFT calculation, this list
            has to be (n,3), n being number of kpoints and 3 being the
            3 different cartesian coordinates.
        
        soundSpeedArray : (n,) float
            A list of energies of ith band cooresponding to the
            kpoints.

        soundSpeed : float
            Value of the sound speed or any sound speed that one wants to
            find the isosurface with.

        reciprocal_lattice : (3,3) float
            Reciprocal lattice of the structure.


        color : TYPE, optional
            DESCRIPTION. The default is None.

        projection_accuracy : TYPE, optional
            DESCRIPTION. The default is 'Normal'.



        """

        self.qpoints = qpoints
        self.soundSpeedArray = soundSpeedArray
        self.reciprocal_lattice = reciprocal_lattice
        self.supercell = supercell
        self.soundSpeed = soundSpeed

        self.brillouin_zone = self._get_brilloin_zone(supercell)

        Isosurface.__init__(self,
                            XYZ=self.qpoints,
                            V=self.soundSpeedArray,
                            isovalue=self.soundSpeed,
                            algorithm='lewiner',
                            interpolation_factor=1,
                            padding= self.supercell,
                            transform_matrix=self.reciprocal_lattice,
                            boundaries=self.brillouin_zone)


    def _get_brilloin_zone(self, supercell):
        return BrillouinZone(self.reciprocal_lattice, supercell)
