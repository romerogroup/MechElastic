# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:11:30 2020

@author: lllan
"""

import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter
from matplotlib import colors as mpcolors
from matplotlib import cm
from .comms import printer
from .soundsurface3d import SoundSurface3D

from .core import BrillouinZone
from pyvistaqt import BackgroundPlotter
import yaml


def soundSurface3D(
        infile='mesh.yaml',
        soundSpeed = None,
        smoothness_interations = 500,
        supercell=[1, 1, 1],
        color='blue',
        background_color='white',
):
    """
    
    Parameters
    ----------
    infile : str, optional (default ``'mesh.yaml'``)
        Path to the mesh.yaml file of the simulation

        e.g. ``procar='~/mesh.yaml'``

    soundSpeed : float, optional (default ``None``)
        Sound speed at which the sound speed surface is created. In other
        words soundSpeed is the isovalue at which the sound speed surface is
        created.

        e.g. ``sound speed=-5.49``
        


    supercell : list int, optional (default ``[1, 1, 1]``)
        If one wants plot more than the 1st brillouin zone, this
        parameter can be used.

        e.g. ``supercell=[2, 2, 2]``

    color : list str, optional
        Color of isosurface. This argument does not work when
        a 3d file is saved. T
        
        e.g. ``colors=['red', 'blue', 'green']``

    background_color : str, optional (default ``white``)
        Defines the background color.

        e.g. ``background_color='gray'``


    Returns
    -------
    s : pyprocar surface object
        The whole fermi surface added bands

    surfaces : list pyprocar surface objects
        list of fermi surface of each band

    """
    def mag(v):
        return (v[0]**2+v[1]**2+v[2]**2)**0.5

    printer.print_mechelastic()

    
    plotter = BackgroundPlotter()

    rf = open(infile, "r")
    data = yaml.safe_load(rf)
    rf.close()

    lattice = np.array(data["lattice"])
    reciprocal_lattice  = np.array(data["reciprocal_lattice"])
    
    qNum = data["nqpoint"]
    bandNum = len(data["phonon"][0]["band"])
    
    
    qp = np.zeros(shape = [qNum , 3 ])
    sound_velocity = np.zeros(shape= [qNum,bandNum,3])
    band = np.zeros(shape= [qNum,bandNum])
    sound_speed = np.zeros(shape = [qNum,bandNum])
    
    for i in range(qNum):
        qp[i,:] = np.array(data["phonon"][i]['q-position'])
        for j in range(bandNum):
            band[i,j] = data["phonon"][i]['band'][j]["frequency"]
            sound_velocity[i,j,:] = np.array(data["phonon"][i]['band'][j]["group_velocity"])
            sound_speed[i,j] = mag(sound_velocity[i,j,:])

    surfaces = []

    bands = np.arange(len(band[0, :]))
    counter = 0
    
    brillouin_zone = BrillouinZone(reciprocal_lattice, [1,1,1])
    
    plotter.add_mesh(brillouin_zone.pyvista_obj,
           style='wireframe',
           line_width=3.5,
           color='black')
    
    for iband in bands:
    
        surface = SoundSurface3D(
                 qpoints=qp,
                 soundSpeedArray=sound_speed[:,iband],
                 soundSpeed = soundSpeed,
                 reciprocal_lattice=reciprocal_lattice,
                 supercell=[1, 1, 1])
        surfaces.append(surface)
        
    nsurface = len(surfaces)
    for isurface in range(nsurface):
        if surfaces[isurface].pyvista_obj != None:
            smoove = surfaces[isurface].pyvista_obj.smooth(n_iter=smoothness_interations )
            plotter.add_mesh(smoove, color = color)

    plotter.set_background(color=background_color)
    


    plotter.add_axes(xlabel='Kx',
               ylabel='Ky',
               zlabel='Kz',
               line_width=6,
               labels_off=False)
