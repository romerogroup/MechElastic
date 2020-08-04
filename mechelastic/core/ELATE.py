# -*- coding: utf-8 -*-

import json
import math
import os
import platform
import random
import re
import sys
import time

from collections import OrderedDict
from io import StringIO
import requests

import numpy as np
from scipy import optimize


__author__ = "Romain Gaillac and François-Xavier Coudert"
__version__ = "2019.01.09"
__license__ = "MIT"




def make3DPlot(func, legend = '', npoints = 200):
    
      str1 = legend.split("\'")[0]
      str2 = legend.split("\'")[1]
    
      u = np.linspace(0, np.pi, npoints)
      v = np.linspace(0, 2*np.pi, 2*npoints)
      r = np.zeros(len(u)*len(v))
    
      dataX = [[0.0 for i in range(len(v))] for j in range(len(u))]
      dataY = [[0.0 for i in range(len(v))] for j in range(len(u))]
      dataZ = [[0.0 for i in range(len(v))] for j in range(len(u))]
      dataR = [["0.0" for i in range(len(v))] for j in range(len(u))]
    
      count = 0
      for cu in range(len(u)):
        for cv in range(len(v)):
          r_tmp = func(u[cu],v[cv])
          z = r_tmp * np.cos(u[cu])
          x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
          y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
          dataX[cu][cv] = x
          dataY[cu][cv] = y
          dataZ[cu][cv] = z
          dataR[cu][cv] = "'E = "+str(float(int(10*r_tmp))/10.0)+" GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
          count = count+1
      return ((dataX,dataY,dataZ,dataR))


def make3DPlotPosNeg(func, legend = '', npoints = 200):

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, 2*np.pi, 2*npoints)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  for cu in range(len(u)):
    for cv in range(len(v)):
      r_tmp = max(0,func(u[cu],v[cv]))
      z = r_tmp * np.cos(u[cu])
      x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
      y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
      dataX1[cu][cv] = x
      dataY1[cu][cv] = y
      dataZ1[cu][cv] = z
      dataR1[cu][cv] = "'"+"\u03B2 = "+str(float(int(10*r_tmp))/10.0)+" TPa'"+"+'-1'.sup()+"+"', \u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  for cu in range(len(u)):
    for cv in range(len(v)):
      r_tmp = max(0,-func(u[cu],v[cv]))
      z = r_tmp * np.cos(u[cu])
      x = r_tmp * np.sin(u[cu]) * np.cos(v[cv])
      y = r_tmp * np.sin(u[cu]) * np.sin(v[cv])
      dataX2[cu][cv] = x
      dataY2[cu][cv] = y
      dataZ2[cu][cv] = z
      dataR2[cu][cv] = "'"+"\u03B2 = -"+str(float(int(10*r_tmp))/10.0)+" TPa'"+"+'-1'.sup()+"+"', \u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  return ((dataX1,dataY1,dataZ1,dataR1),(dataX2,dataY2,dataZ2,dataR2))

def make3DPlot2(func, legend = '', npoints = 50):

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, np.pi, npoints)
  w = [v[i]+np.pi for i in range(1,len(v))]
  v = np.append(v, w)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  r = [0.0,0.0,np.pi/2.0,np.pi/2.0]
  for cu in range(len(u)):
    for cv in range(len(v)):

      r = func(u[cu],v[cv],r[2],r[3])
      z = np.cos(u[cu])
      x = np.sin(u[cu]) * np.cos(v[cv])
      y = np.sin(u[cu]) * np.sin(v[cv])

      r1_tmp = r[0]
      z1 = r1_tmp * z
      x1 = r1_tmp * x
      y1 = r1_tmp * y
      dataX1[cu][cv] = x1
      dataY1[cu][cv] = y1
      dataZ1[cu][cv] = z1
      dataR1[cu][cv] = "'"+"G'"+"+'min'.sub()+"+"' = "+str(float(int(10*r1_tmp))/10.0)+"GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"

      r2_tmp = r[1]
      z2 = r2_tmp * z
      x2 = r2_tmp * x
      y2 = r2_tmp * y
      dataX2[cu][cv] = x2
      dataY2[cu][cv] = y2
      dataZ2[cu][cv] = z2
      dataR2[cu][cv] = "'"+"G'"+"+'max'.sub()+"+"' = "+str(float(int(10*r1_tmp))/10.0)+"GPa, "+"\u03B8 = "+str(float(int(10*u[cu]*180/np.pi))/10.0)+"\u00B0, "+"\u03c6 = "+str(float(int(10*v[cv]*180/np.pi))/10.0)+"\u00B0'"
      count = count+1

  i = random.randint(0, 100000)
  return ((dataX1,dataY1,dataZ1,dataR1),(dataX2,dataY2,dataZ2,dataR2))


def make3DPlot3(func, legend = '', width = 600, height = 600, npoints = 50):

  str1 = legend.split("\'")[0]
  str2 = legend.split("\'")[1]

  u = np.linspace(0, np.pi, npoints)
  v = np.linspace(0, np.pi, npoints)
  w = [v[i]+np.pi for i in range(1,len(v))]
  v = np.append(v, w)

  dataX1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ1 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR1 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ2 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR2 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  dataX3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataY3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataZ3 = [[0.0 for i in range(len(v))] for j in range(len(u))]
  dataR3 = [["0.0" for i in range(len(v))] for j in range(len(u))]

  count = 0
  r = [0.0, 0.0, 0.0, np.pi/2.0, np.pi/2.0]
  ruv = [[r for i in range(len(u))] for j in range(len(v))]
  for cu in range(len(u)):
    for cv in range(len(v)):
       ruv[cv][cu] = func(u[cu],v[cv],r[3],r[4])

  for cu in range(len(u)):
    for cv in range(len(v)):

      z = np.cos(u[cu])
      x = np.sin(u[cu]) * np.cos(v[cv])
      y = np.sin(u[cu]) * np.sin(v[cv])

      r = ruv[cv][cu]
      r1_tmp = r[0]
      dataX1[cu][cv] = r1_tmp * x
      dataY1[cu][cv] = r1_tmp * y
      dataZ1[cu][cv] = r1_tmp * z
      dataR1[cu][cv] = "'"+"\u03BD'"+"+'min'.sub()+"+"' = "+str(float(int(100*r1_tmp))/100.0)+", "+"\u03B8 = "+str(float(int(100*u[cu]*180/np.pi))/100.0)+"\u00B0, "+"\u03c6 = "+str(float(int(100*v[cv]*180/np.pi))/100.0)+"\u00B0'"

      r2_tmp = r[1]
      dataX2[cu][cv] = r2_tmp * x
      dataY2[cu][cv] = r2_tmp * y
      dataZ2[cu][cv] = r2_tmp * z
      dataR2[cu][cv] = float(int(100*r2_tmp))/100.0 

      r3_tmp = r[2]
      dataX3[cu][cv] = r3_tmp * x
      dataY3[cu][cv] = r3_tmp * y
      dataZ3[cu][cv] = r3_tmp * z
      dataR3[cu][cv] = "'"+"\u03BD'"+"+'max'.sub()+"+"' = "+str(float(int(100*r3_tmp))/100.0)+", "+"\u03B8 = "+str(float(int(100*u[cu]*180/np.pi))/100.0)+"\u00B0, "+"\u03c6 = "+str(float(int(100*v[cv]*180/np.pi))/100.0)+"\u00B0'"
      count = count+1

  return ((dataX1,dataY1,dataZ1,dataR1),(dataX2,dataY2,dataZ2,dataR2),(dataX3,dataY3,dataZ3,dataR3))




# Polar plot functions
################################################################################################


def makePolarPlot(func, legend = '', p="xy", npoints = 90):
  u = np.linspace(0, 2*np.pi, npoints)
  r = list(map(func, u))
  if (p=="xy"):
    x = r * np.cos(u)
    y = r * np.sin(u)
  else:
    y = r * np.cos(u)
    x = r * np.sin(u)

  return ((x,y))

def makePolarPlotPosNeg(func, legend = '', p="xy", npoints = 90):
  
  u = np.linspace(0, 2*np.pi, npoints)
  r = list(map(lambda x: max(0, func(x)), u))
  if (p=="xy"):
    x1 = r * np.cos(u)
    y1 = r * np.sin(u)
  else:
    y1 = r * np.cos(u)
    x1 = r * np.sin(u)
  r = list(map(lambda x: max(0, -func(x)), u))
  if (p=="xy"):
    x2 = r * np.cos(u)
    y2 = r * np.sin(u)
  else:
    y2 = r * np.cos(u)
    x2 = r * np.sin(u)

  return ((x1,y1),(x2,y2))

def makePolarPlot2(func, legend = '', p="xy", npoints = 61):


  u = np.linspace(0, 2*np.pi, npoints)
  r = list(map(func, u))

  if (p=="xy"):
    x1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
    y1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
    x2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
    y2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
  else:
    y1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
    x1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
    y2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
    x2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])

  return ((x1,y1),(x2,y2))

def makePolarPlot3(func, legend = '', p="xy", npoints = 61):


  u = np.linspace(0, 2*np.pi, npoints)
  r = list(map(func, u))

  if (p=="xy"):
    x1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
    y1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
    x2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
    y2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
    x3 = np.array([ ir[2] * np.cos(iu) for ir, iu in zip(r,u) ])
    y3 = np.array([ ir[2] * np.sin(iu) for ir, iu in zip(r,u) ])
  else:
    y1 = np.array([ ir[0] * np.cos(iu) for ir, iu in zip(r,u) ])
    x1 = np.array([ ir[0] * np.sin(iu) for ir, iu in zip(r,u) ])
    y2 = np.array([ ir[1] * np.cos(iu) for ir, iu in zip(r,u) ])
    x2 = np.array([ ir[1] * np.sin(iu) for ir, iu in zip(r,u) ])
    y3 = np.array([ ir[2] * np.cos(iu) for ir, iu in zip(r,u) ])
    x3 = np.array([ ir[2] * np.sin(iu) for ir, iu in zip(r,u) ])

  return ((x1,y1),(x2,y2),(x3,y3))


########################################################################################################

def dirVec(theta, phi):
  return [ math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta) ]

def dirVec1(theta, phi, chi):
  return [ math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta) ]

def dirVec2(theta, phi, chi):
  return [ math.cos(theta)*math.cos(phi)*math.cos(chi) - math.sin(phi)*math.sin(chi),
          math.cos(theta)*math.sin(phi)*math.cos(chi) + math.cos(phi)*math.sin(chi),
          - math.sin(theta)*math.cos(chi) ]


# Functions to minimize/maximize
def minimize(func, dim):
  if dim == 2:
    r = ((0, np.pi), (0, np.pi))
    n = 25
  elif dim == 3:
    r = ((0, np.pi), (0, np.pi), (0, np.pi))
    n = 10

  # TODO -- try basin hopping or annealing
  return optimize.brute(func, r, Ns = n, full_output = True, finish = optimize.fmin)[0:2]

def maximize(func, dim):
  res = minimize(lambda x: -func(x), dim)
  return (res[0], -res[1])


class ELATE:
    def __init__(self, s):
        
        self.elas = Elastic(s)
        self.elasList = s
        
        
        
        minE = minimize(self.elas.Young, 2)
        maxE = maximize(self.elas.Young, 2)
        minLC = minimize(self.elas.LC, 2)
        maxLC = maximize(self.elas.LC, 2)
        minG = minimize(self.elas.shear, 3)
        maxG = maximize(self.elas.shear, 3)
        minNu = minimize(self.elas.Poisson, 3)
        maxNu = maximize(self.elas.Poisson, 3)
        
        self.voigtE = self.elas.averages()[0][0]
        self.reussE = self.elas.averages()[1][0]
        self.hillE  = self.elas.averages()[2][0]
        self.max_E = maxE[1]
        self.min_E = minE[1]
        self.min_axis_E = tuple(dirVec(*minE[0]))
        self.max_axis_E = tuple(dirVec(*maxE[0]))
        self.anis_E = maxE[1]/minE[1]
        
        self.voigtLC = self.elas.averages()[0][1]
        self.reussLC = self.elas.averages()[1][1]
        self.hillLC  = self.elas.averages()[2][1]
        self.max_LC = maxLC[1]
        self.min_LC = minLC[1]
        self.min_axis_LC = tuple(dirVec(*minLC[0]))
        self.max_axis_LC = tuple(dirVec(*maxLC[0]))
        if minLC[1] > 0:
            self.anis_LC = maxLC[1]/minLC[1]
        else:
            self.anis_LC = "&infin;"
            
        self.voigtShear = self.elas.averages()[0][2]
        self.reussShear = self.elas.averages()[1][2]
        self.hillShear  = self.elas.averages()[2][2]
        self.max_Shear = maxG[1]
        self.min_Shear = minG[1]
        self.min_axis_Shear = tuple(dirVec1(*minG[0]))
        self.max_axis_Shear = tuple(dirVec1(*maxG[0]))
        self.mix_2nd_axis_Shear = tuple(dirVec2(*minG[0]))
        self.max_2nd_axis_Shear = tuple(dirVec2(*maxG[0]))
        self.anis_Shear = maxG[1]/minG[1]
        
        self.voigtPoisson = self.elas.averages()[0][3]
        self.reussPoisson = self.elas.averages()[1][3]
        self.hillPoisson  = self.elas.averages()[2][3]
        self.max_Poisson = maxNu[1]
        self.min_Poisson = minNu[1]
        self.min_axis_Poisson = tuple(dirVec1(*minNu[0]))
        self.max_axis_Poisson = tuple(dirVec1(*maxNu[0]))
        self.min_2nd_axis_Poisson = tuple(dirVec2(*minNu[0]))
        self.max_2nd_axis_Poisson = tuple(dirVec2(*maxNu[0]))
        if minNu[1]*maxNu[1] > 0:
            self.anis_Poisson = maxNu[1]/minNu[1]
        else:
            self.anix_Poisson = "&infin;"
            
            
        
    def YOUNG2D(self):
        data1 = makePolarPlot(lambda x: self.elas.Young([np.pi / 2, x]), "Young's modulus in (xy) plane", "xy")
        data2 = makePolarPlot(lambda x: self.elas.Young([x, 0]), "Young's modulus in (xz) plane", "xz")
        data3 = makePolarPlot(lambda x: self.elas.Young([x, np.pi / 2]), "Young's modulus in (yz) plane", "yz")
        return (data1,data2,data3)
    
    def YOUNG3D(self):
    
          if self.elas.isOrthorhombic():
            self.elas = ElasticOrtho(self.elas)
        
        
          data = make3DPlot(lambda x, y: self.elas.Young_2(x, y), "Young's modulus")
        
          return data
    
    def LC2D(self):
        data1 = makePolarPlotPosNeg(lambda x: self.elas.LC([np.pi / 2, x]), "linear compressibility in (xy) plane", "xy")
        data2 = makePolarPlotPosNeg(lambda x: self.elas.LC([x, 0]), "linear compressibility in (xz) plane", "xz")
        data3 = makePolarPlotPosNeg(lambda x: self.elas.LC([x, np.pi / 2]),  "linear compressibility in (yz) plane", "yz")
        return (data1,data2,data3)
    
    def LC3D(self):
    
          if self.elas.isOrthorhombic():
            self.elas = ElasticOrtho(self.elas)
        
          data = make3DPlotPosNeg(lambda x, y: self.elas.LC_2(x, y), "Linear compressiblity")
        
          return data
    
    def SHEAR2D(self):
        data1 = makePolarPlot2(lambda x: self.elas.shear2D([np.pi / 2, x]), "Shear modulus in (xy) plane", "xy")
        data2 = makePolarPlot2(lambda x: self.elas.shear2D([x, 0]), "Shear modulus in (xz) plane", "xz")
        data3 = makePolarPlot2(lambda x: self.elas.shear2D([x, np.pi / 2]), "Shear modulus in (yz) plane", "yz")
        return (data1,data2,data3)
    
    def SHEAR3D(self):
    
          if self.elas.isOrthorhombic():
            self.elas = ElasticOrtho(self.elas)
        
          data = make3DPlot2(lambda x, y, g1, g2: self.elas.shear3D(x, y, g1, g2), "Shear modulus")
        
        
          return data
    
    def POISSON2D(self):
        data1 = makePolarPlot3(lambda x: self.elas.Poisson2D([np.pi / 2, x]), "Poisson's ratio in (xy) plane", "xy")
        data2 = makePolarPlot3(lambda x: self.elas.Poisson2D([x, 0]), "Poisson's ratio in (xz) plane", "xz")
        data3 = makePolarPlot3(lambda x: self.elas.Poisson2D([x, np.pi / 2]), "Poisson's ratio in (yz) plane", "yz")
        return (data1,data2,data3)
    
    def POISSON3D(self):
    
          if self.elas.isOrthorhombic():
            self.elas = ElasticOrtho(self.elas)
        
          data = make3DPlot3(lambda x, y, g1, g2: self.elas.poisson3D(x, y, g1, g2), "Poisson's ratio")
        
          return data
    
    
    ##############################################################################
    # Plotting functions
    #############################################################################
    def plot_3D(self,elastic_calc = ''):
        import pyvista as pv 
        from pyvistaqt import BackgroundPlotter
        
        plotter = BackgroundPlotter()
        
        x = None
        y = None
        z = None
        r = None
        
        if elastic_calc == 'POISSON':
            func = self.POISSON3D()
            colors = ['red','green','blue']
            for ix,icolor in zip(range(len(func)),colors):
                x = np.array(func[ix][0])
                y = np.array(func[ix][1])
                z = np.array(func[ix][2])
                r = np.array(func[ix][3])
                if np.all((func[ix][0] == 0)):
                    continue
                else:
                    grid = pv.StructuredGrid(x, y, z)
                    if (ix == 2):
                        plotter.add_mesh(grid,opacity = 0.25, color = icolor )
                    else:
                        plotter.add_mesh(grid, opacity = 0.50, color = icolor )
                    
        elif elastic_calc == 'SHEAR':
            func = self.SHEAR3D()
            colors = ['green','blue']
            for ix,icolor in zip(range(len(func)),colors):
                x = np.array(func[ix][0])
                y = np.array(func[ix][1])
                z = np.array(func[ix][2])
                r = np.array(func[ix][3])
                if np.all((func[ix][0] == 0)):
                    continue
                else:
                    grid = pv.StructuredGrid(x, y, z)
                    if (ix == 2):
                        plotter.add_mesh(grid,opacity = 0.25, color = icolor )
                    else:
                        plotter.add_mesh(grid, opacity = 0.50, color = icolor )
                        
        elif elastic_calc == 'LC':
            func = self.LC3D()
            colors = ['green','red']
            for ix,icolor in zip(range(len(func)),colors):
                x = np.array(func[ix][0])
                y = np.array(func[ix][1])
                z = np.array(func[ix][2])
                r = np.array(func[ix][3])
                if np.all((func[ix][0] == 0)):
                    continue
                else:
                    grid = pv.StructuredGrid(x, y, z)
                    if (ix == 2):
                        plotter.add_mesh(grid,opacity = 0.25, color = icolor )
                    else:
                        plotter.add_mesh(grid, opacity = 0.50, color = icolor )
    
        elif elastic_calc == 'YOUNG':
            func = self.YOUNG3D()
            colors = ['green']
            for ix,icolor in zip(range(len(func)),colors):
                x = np.array(func[ix][0])
                y = np.array(func[ix][1])
                z = np.array(func[ix][2])
                r = np.array(func[ix][3])
                if np.all((func[ix][0] == 0)):
                    continue
                else:
                    grid = pv.StructuredGrid(x, y, z)
                    if (ix == 2):
                        plotter.add_mesh(grid,opacity = 0.25, color = icolor )
                    else:
                        plotter.add_mesh(grid, opacity = 0.50, color = icolor )
                        
        plotter.add_axes()
        plotter.show_grid()
        
    def plot_2D(self,elastic_calc =''):
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        subTitles = ['XY','XZ','YZ']
        
        if elastic_calc == 'POISSON':
            func = self.POISSON2D()
            colors = ['red','green','blue']
            labels = ['Poisson - Neg', 'Poisson - Pos','Poisson - Max']
            fig.suptitle('Poisson\'s ratio')
            for iplane,title in zip(range(len(func)),subTitles):
                
                ax = fig.add_subplot(1,3,iplane+1)
                ax.set_title(title)
                ax.get_yaxis().set_visible(False)
                for iplot,color in zip(range(len(func[0])),colors):
                    plt.plot(func[iplane][iplot][0],func[iplane][iplot][1], color = color)
                    
        elif elastic_calc == 'SHEAR':
            func = self.SHEAR2D()
            colors = ['green','blue']
            labels = ['Shear - Max', 'Shear -']
            fig.suptitle('Shear modulus')
            for iplane,title in zip(range(len(func)),subTitles):
                
                ax = fig.add_subplot(1,3,iplane+1)
                ax.set_title(title)
                ax.get_yaxis().set_visible(False)
                for iplot,color in zip(range(len(func[0])),colors):
                    plt.plot(func[iplane][iplot][0],func[iplane][iplot][1], color = color)
        
        elif elastic_calc == 'LC':
            func = self.LC2D()
            colors = ['green','red']
            labels = ['LC-positive','LC-negative']
            fig.suptitle('Linear Compression')
            for iplane,title in zip(range(len(func)),subTitles):
            
                ax = fig.add_subplot(1,3,iplane+1)
                ax.set_title(title)
                ax.get_yaxis().set_visible(False)
                for iplot,color,label in zip(range(len(func[0])),colors,labels):
                    plt.plot(func[iplane][iplot][0],func[iplane][iplot][1], color = color,label = label)
    
                    
        elif elastic_calc == 'YOUNG':
            func = self.YOUNG2D()
            color = 'green'
            label = 'YOUNG'
            fig.suptitle('Young\'s Modulus')
            for iplane, title in zip(range(len(func)),subTitles):
                ax = fig.add_subplot(1,3,iplane+1)
                ax.set_title(title)
                ax.get_yaxis().set_visible(False)
                plt.plot(func[iplane][0],func[iplane][1], color = color, label = label)
       
    
    def print_properties(self):
#        print("\n \n Input: compliance matrix(coefficients in GPa) of\n")
#        for x in self.elasList:
#            print(str(x))
            
        print("\n")
        print("\n")
        print("Average properties")
        
        print("\n \n                         Voigt     Reuss    Hill")
        print("-------------------------------------------------------")
        print(
            "Bulk modulus   (GPa)  %9.3f %9.3f %9.3f "
            % (self.elas.averages()[0][0], self.elas.averages()[1][0], self.elas.averages()[2][0])
        )
        print(
            "Shear modulus  (GPa)  %9.3f %9.3f %9.3f "
            % (self.elas.averages()[0][1], self.elas.averages()[1][1], self.elas.averages()[2][1])
        )
        print(
            "Young modulus  (GPa)  %9.3f %9.3f %9.3f "
            % (self.elas.averages()[0][2], self.elas.averages()[1][2], self.elas.averages()[2][2])
        )
        print(
            "Poisson ratio         %9.3f %9.3f %9.3f "
            % (self.elas.averages()[0][3], self.elas.averages()[1][3], self.elas.averages()[2][3]))
        
        print("-------------------------------------------------------")
        #printer.printMatix(elastic_rowsList)
        print("\n")
        print("Eigenvalues of compliance matrix")
        print("\n \n   lamda_1  lamda_2  lamda_3  lamda_4  lamda_5  lamda_6")    
        print("---------------------------------------------------------------")
        eigenval = sorted(np.linalg.eig(self.elas.CVoigt)[0])
        print("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f" % tuple(eigenval))
        print("--------------------------------------------------------------")
        if eigenval[0] <= 0:
            print('<div class="error">Stiffness matrix is not definite positive, crystal is mechanically unstable<br/>')
        print("\n")
        print("Variations of the elastic moduli")
#        print("\n \n   Young's Modulus                Linear Compressibility                   Shear Modulus              Poisson's Ratio")
#        print("    Emin         Emax             betaMin            betaMax               Gmin        Gmax           Nu_min      Nu_max")
#        print("----------------------------------------------------------------------------------------------")
        minE = minimize(self.elas.Young, 2)
        maxE = maximize(self.elas.Young, 2)
        minLC = minimize(self.elas.LC, 2)
        maxLC = maximize(self.elas.LC, 2)
        minG = minimize(self.elas.shear, 3)
        maxG = maximize(self.elas.shear, 3)
        minNu = minimize(self.elas.Poisson, 3)
        maxNu = maximize(self.elas.Poisson, 3)
        #print(type((minE,maxE,minLC,maxLC,minG,maxG,minNu,maxNu)))
        #print("%9.3f GPa  %9.3f GPa  %9.3f TPa^-1   %9.3f TPa^-1    %9.3f GPa  %9.3f GPa  %9.3f  %9.3f" % (minE[1],maxE[1],minLC[1],maxLC[1],minG[1],maxG[1],minNu[1],maxNu[1]))
        
        anisE = '%8.4f' % (maxE[1]/minE[1])
        anisLC = ('%8.4f' % (maxLC[1]/minLC[1])) if minLC[1] > 0 else "&infin;"
        anisG = '%8.4f' % (maxG[1]/minG[1])
        anisNu = ('%8.4f' % (maxNu[1]/minNu[1])) if minNu[1]*maxNu[1] > 0 else "&infin;"
        
        #print("       %s                     %s                                %s                         %s" % (anisE,anisLC,anisG,anisNu))
        #print(dirVec1(*minG[0]))
        
        print("\n \n                                  Min       Max   ||   Anisotropy")
        print("----------------------------------------------------------------------------------------")
        print(
            "Young's Modulus     (GPa)      %9.3f %9.3f || %s "
            % (minE[1], maxE[1], anisE)
        )
        print(
            "         Min Axis:    %s      \n         Max Axis:    %s "
            % (str(tuple(dirVec(*minE[0]))), str(tuple(dirVec(*maxE[0]))))
        )
        print("----------------------------------------------------------------------------------------")   
        print(
            "Linear Compression   (TPa^-1)  %9.3f %9.3f || %s "
            % (minE[1], maxE[1], anisLC)
        )
        print(
            "         Min Axis:    %s      \n         Max Axis:    %s "
            % (str(tuple(dirVec(*minLC[0]))), str(tuple(dirVec(*maxLC[0]))))
        )
        print("----------------------------------------------------------------------------------------")   
        print(
            "Shear Modulus   (GPa)          %9.3f %9.3f || %s "
            % (minE[1], maxE[1], anisG)
        )
        print(
            "         Min Axis:    %s      \n         Max Axis:    %s "
            % (str(tuple(dirVec1(*minG[0]))), str(tuple(dirVec1(*maxG[0]))))
        )
        print(
            "         Second Min Axis:    %s      \n         Second Max Axis:    %s "
            % (str(tuple(dirVec2(*minG[0]))), str(tuple(dirVec2(*maxG[0]))))
        )
        print("----------------------------------------------------------------------------------------")   
        print(
            "Poisson Ratio                 %9.3f %9.3f || %s "
            % (minE[1], maxE[1], anisNu)
        )
        print(
            "         Min Axis:    %s      \n         Max Axis:    %s "
            % (str(tuple(dirVec1(*minNu[0]))), str(tuple(dirVec1(*maxNu[0]))))
        )
        print(
            "         Second Min Axis:    %s      \n         Second Max Axis:    %s "
            % (str(tuple(dirVec2(*minNu[0]))), str(tuple(dirVec2(*maxNu[0]))))
        )
#        print(
#            "Shear modulus  (GPa)  %9.3f %9.3f %9.3f "
#            % (self.elas.averages()[0][1], self.elas.averages()[1][1], self.elas.averages()[2][1])
#        )
#        print(
#            "Young modulus  (GPa)  %9.3f %9.3f %9.3f "
#            % (self.elas.averages()[0][2], self.elas.averages()[1][2], self.elas.averages()[2][2])
#        )
#        print(
#            "Poisson ratio         %9.3f %9.3f %9.3f "
#            % (self.elas.averages()[0][3], self.elas.averages()[2][3], self.elas.averages()[2][3]))
#        
#        print("-------------------------------------------------------")
#        print("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f" % tuple(dirVec(*minE[0])) + tuple(dirVec(*maxE[0])) + tuple(dirVec(*minLC[0])) + tuple(dirVec(*maxLC[0])) + tuple(dirVec1(*minG[0])) +  tuple(dirVec1(*maxG[0])) + tuple(dirVec1(*minNu[0])) + tuple(dirVec1(*maxNu[0])))
#        print("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f" % tuple(dirVec2(*minG[0])) + tuple(dirVec2(*maxG[0])) + tuple(dirVec2(*minNu[0])) + tuple(dirVec2(*maxNu[0])))



class Elastic:
  """An elastic tensor, along with methods to access it"""

  def __init__(self, s):
    """Initialize the elastic tensor from a string"""
    
    if not s:
      raise ValueError("no matrix was provided")

    # Argument can be a 6-line string, a list of list, or a string representation of the list of list
    try:
      if type(json.loads(s)) == list: s = json.loads(s)
    except:
      pass

    if type(s) == str:
      # Remove braces and pipes
      s = s.replace("|", " ").replace("(", " ").replace(")", " ")

      # Remove empty lines
      lines = [line for line in s.split('\n') if line.strip()]
      if len(lines) != 6:
        raise ValueError("should have six rows")

      # Convert to float
      try:
        mat = [list(map(float, line.split())) for line in lines]
      except:
        raise ValueError("not all entries are numbers")
    elif type(s) == list:
      # If we already have a list, simply use it
      mat = s
    else:
      raise ValueError("invalid argument as matrix")

    # Make it into a square matrix
    mat = np.array(mat)
    if mat.shape != (6,6):
      # Is it upper triangular?
      if list(map(len, mat)) == [6,5,4,3,2,1]:
        mat = [ [0]*i + mat[i] for i in range(6) ]
        mat = np.array(mat)

      # Is it lower triangular?
      if list(map(len, mat)) == [1,2,3,4,5,6]:
        mat = [ mat[i] + [0]*(5-i) for i in range(6) ]
        mat = np.array(mat)

    if mat.shape != (6,6):
      raise ValueError("should be a square matrix")

    # Check that is is symmetric, or make it symmetric
    if np.linalg.norm(np.tril(mat, -1)) == 0:
      mat = mat + np.triu(mat, 1).transpose()
    if np.linalg.norm(np.triu(mat, 1)) == 0:
      mat = mat + np.tril(mat, -1).transpose()
    if np.linalg.norm(mat - mat.transpose()) > 1e-3:
      raise ValueError("should be symmetric, or triangular")
    elif np.linalg.norm(mat - mat.transpose()) > 0:
      mat = 0.5 * (mat + mat.transpose())

    # Store it
    self.CVoigt = mat

    # Put it in a more useful representation
    try:
      self.SVoigt = np.linalg.inv(self.CVoigt)
    except:
      raise ValueError("matrix is singular")

    VoigtMat = [[0, 5, 4], [5, 1, 3], [4, 3, 2]]
    def SVoigtCoeff(p,q): return 1. / ((1+p//3)*(1+q//3))

    self.Smat = [[[[ SVoigtCoeff(VoigtMat[i][j], VoigtMat[k][l]) * self.SVoigt[VoigtMat[i][j]][VoigtMat[k][l]]
                     for i in range(3) ] for j in range(3) ] for k in range(3) ] for l in range(3) ]
    return

  def isOrthorhombic(self):
    def iszero(x): return (abs(x) < 1.e-3)
    return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
            and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
            and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
            and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5]))

  def isCubic(self):
      def iszero(x): return (abs(x) < 1.e-3)
      return (iszero(self.CVoigt[0][3]) and iszero(self.CVoigt[0][4]) and iszero(self.CVoigt[0][5])
                and iszero(self.CVoigt[1][3]) and iszero(self.CVoigt[1][4]) and iszero(self.CVoigt[1][5])
                and iszero(self.CVoigt[2][3]) and iszero(self.CVoigt[2][4]) and iszero(self.CVoigt[2][5])
                and iszero(self.CVoigt[3][4]) and iszero(self.CVoigt[3][5]) and iszero(self.CVoigt[4][5])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[0][0] - self.CVoigt[1][1]) and iszero(self.CVoigt[0][0] - self.CVoigt[2][2])
                and iszero(self.CVoigt[3][3] - self.CVoigt[4][4]) and iszero(self.CVoigt[3][3] - self.CVoigt[5][5])
                and iszero(self.CVoigt[0][1] - self.CVoigt[0][2]) and iszero(self.CVoigt[0][1] - self.CVoigt[1][2]))

  def Young(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

  def Young_2(self,x,y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/r

  def LC(self, x):
        a = dirVec(x[0], x[1])
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

  def LC_2(self, x, y):
        a = dirVec(x, y)
        r = sum([ a[i]*a[j] * self.Smat[i][j][k][k]
                  for i in range(3) for j in range(3) for k in range(3) ])
        return 1000 * r

  def shear(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r = sum([ a[i]*b[j]*a[k]*b[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return 1/(4*r)

  def Poisson(self, x):
        a = dirVec(x[0], x[1])
        b = dirVec2(x[0], x[1], x[2])
        r1 = sum([ a[i]*a[j]*b[k]*b[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        r2 = sum([ a[i]*a[j]*a[k]*a[l] * self.Smat[i][j][k][l]
                  for i in range(3) for j in range(3) for k in range(3) for l in range(3) ])
        return -r1/r2

  def averages(self):
        A = (self.CVoigt[0][0] + self.CVoigt[1][1] + self.CVoigt[2][2]) / 3
        B = (self.CVoigt[1][2] + self.CVoigt[0][2] + self.CVoigt[0][1]) / 3
        C = (self.CVoigt[3][3] + self.CVoigt[4][4] + self.CVoigt[5][5]) / 3
        a = (self.SVoigt[0][0] + self.SVoigt[1][1] + self.SVoigt[2][2]) / 3
        b = (self.SVoigt[1][2] + self.SVoigt[0][2] + self.SVoigt[0][1]) / 3
        c = (self.SVoigt[3][3] + self.SVoigt[4][4] + self.SVoigt[5][5]) / 3
    
        KV = (A + 2*B) / 3
        GV = (A - B + 3*C) / 5
    
        KR = 1 / (3*a + 6*b)
        GR = 5 / (4*a - 4*b + 3*c)
    
        KH = (KV + KR) / 2
        GH = (GV + GR) / 2
    
        return [ [KV, 1/(1/(3*GV) + 1/(9*KV)), GV, (1 - 3*GV/(3*KV+GV))/2],
                 [KR, 1/(1/(3*GR) + 1/(9*KR)), GR, (1 - 3*GR/(3*KR+GR))/2],
                 [KH, 1/(1/(3*GH) + 1/(9*KH)), GH, (1 - 3*GH/(3*KH+GH))/2] ]

  def shear2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.shear([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.shear([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun))

  def shear3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.005
        def func1(z): return self.shear([x, y, z])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.shear([x, y, z])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        return (float(r1.fun), -float(r2.fun), float(r1.x), float(r2.x))

  def Poisson2D(self, x):
        ftol = 0.001
        xtol = 0.01
        def func1(z): return self.Poisson([x[0], x[1], z])
        r1 = optimize.minimize(func1, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.Poisson([x[0], x[1], z])
        r2 = optimize.minimize(func2, np.pi/2.0, args=(), method = 'Powell', options={"xtol":xtol, "ftol":ftol})#, bounds=[(0.0,np.pi)])
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun))

  def poisson3D(self, x, y, guess1 = np.pi/2.0, guess2 = np.pi/2.0):
        tol = 0.005
        def func1(z): return self.Poisson([x, y, z])
        r1 = optimize.minimize(func1, guess1, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        def func2(z): return -self.Poisson([x, y, z])
        r2 = optimize.minimize(func2, guess2, args=(), method = 'COBYLA', options={"tol":tol})#, bounds=[(0.0,np.pi)])
        return (min(0,float(r1.fun)), max(0,float(r1.fun)), -float(r2.fun), float(r1.x), float(r2.x))


class ElasticOrtho(Elastic):
  """An elastic tensor, for the specific case of an orthorhombic system"""

  def __init__(self, arg):
        """Initialize from a matrix, or from an Elastic object"""
        if type(arg) == str:
          Elastic.__init__(self, arg)
        elif isinstance(arg, Elastic):
          self.CVoigt = arg.CVoigt
          self.SVoigt = arg.SVoigt
          self.Smat = arg.Smat
        else:
          raise TypeError("ElasticOrtho constructor argument should be string or Elastic object")

  def Young(self, x):
        ct2 = math.cos(x[0])**2
        st2 = 1 - ct2
        cf2 = math.cos(x[1])**2
        sf2 = 1 - cf2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1/(ct2**2*s33 + 2*cf2*ct2*s13*st2 + cf2*ct2*s55*st2 + 2*ct2*s23*sf2*st2 + ct2*s44*sf2*st2 + cf2**2*s11*st2**2 + 2*cf2*s12*sf2*st2**2 + cf2*s66*sf2*st2**2 + s22*sf2**2*st2**2)

  def LC(self, x):
        ct2 = math.cos(x[0])**2
        cf2 = math.cos(x[1])**2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        return 1000 * (ct2 * (s13 + s23 + s33) + (cf2 * (s11 + s12 + s13) + (s12 + s22 + s23) * (1 - cf2)) * (1 - ct2))

  def shear(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        sf2 = sf*sf
        cx = math.cos(x[2])
        cx2 = cx*cx
        sx = math.sin(x[2])
        sx2 = 1 - cx2
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
        r = (
              ct2*ct2*cx2*s44*sf2 + cx2*s44*sf2*st2*st2 + 4*cf**3*ct*cx*(-2*s11 + 2*s12 + s66)*sf*st2*sx
              + 2*cf*ct*cx*sf*(ct2*(s44 - s55) + (4*s13 - 4*s23 - s44 + s55 - 4*s12*sf2 + 4*s22*sf2 - 2*s66*sf2)*st2)*sx
              + s66*sf2*sf2*st2*sx2 + cf**4*st2*(4*ct2*cx2*s11 + s66*sx2)
              + ct2*(2*cx2*(2*s33 + sf2*(-4*s23 - s44 + 2*s22*sf2))*st2 + s55*sf2*sx2)
              + cf**2*(ct2*ct2*cx2*s55 + ct2*(-2*cx2*(4*s13 + s55 - 2*(2*s12 + s66)*sf2)*st2 + s44*sx2)
                       + st2*(cx2*s55*st2 + 2*(2*s11 - 4*s12 + 2*s22 - s66)*sf2*sx2))
            )
        return 1/r

  def Poisson(self, x):
        ct = math.cos(x[0])
        ct2 = ct*ct
        st2 = 1 - ct2
        cf = math.cos(x[1])
        sf = math.sin(x[1])
        cx = math.cos(x[2])
        sx = math.sin(x[2])
        s11 = self.Smat[0][0][0][0]
        s22 = self.Smat[1][1][1][1]
        s33 = self.Smat[2][2][2][2]
        s44 = 4 * self.Smat[1][2][1][2]
        s55 = 4 * self.Smat[0][2][0][2]
        s66 = 4 * self.Smat[0][1][0][1]
        s12 = self.Smat[0][0][1][1]
        s13 = self.Smat[0][0][2][2]
        s23 = self.Smat[1][1][2][2]
    
        return (
      (-(ct**2*cx**2*s33*st2) - cf**2*cx**2*s13*st2*st2 - cx**2*s23*sf**2*st2*st2 + ct*cx*s44*sf*st2*(ct*cx*sf + cf*sx) -
              ct**2*s23*(ct*cx*sf + cf*sx)**2 - cf**2*s12*st2*(ct*cx*sf + cf*sx)**2 - s22*sf**2*st2*(ct*cx*sf + cf*sx)**2 +
              cf*ct*cx*s55*st2*(cf*ct*cx - sf*sx) - cf*s66*sf*st2*(ct*cx*sf + cf*sx)*(cf*ct*cx - sf*sx) -
              ct**2*s13*(cf*ct*cx - sf*sx)**2 - cf**2*s11*st2*(cf*ct*cx - sf*sx)**2 - s12*sf**2*st2*(cf*ct*cx - sf*sx)**2)/
            (ct**4*s33 + 2*cf**2*ct**2*s13*st2 + cf**2*ct**2*s55*st2 + 2*ct**2*s23*sf**2*st2 + ct**2*s44*sf**2*st2 +
              cf**4*s11*st2*st2 + 2*cf**2*s12*sf**2*st2*st2 + cf**2*s66*sf**2*st2*st2 + s22*sf**4*st2*st2)
        )
#########################################################################################
# Data genertaation
################################################################################################
# ELATE : basic usage of the tool, only 2D plots
# YOUNG3D : visualize Young's modulus in 3D
# LC3D : visualize Linear compressiblity in 3D
# SHEAR3D : visualize Shear modulus in 3D
# POISSON3D : visualize Poisson ratio in 3D
################################################################################################




      