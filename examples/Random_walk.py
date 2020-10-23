# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 12:52:52 2020

@author: lllan
"""

import numpy as np
import matplotlib.pyplot as plt

import random
def rwalk3D(size):
    x      = 0.0
    y      = 0.0
    z      = 0.0
    outArr = []
    xArr   = np.array([x])
    yArr   = np.array([y])
    zArr   = np.array([z])
    
    for i in range(0, size+1):
        x += (random.random() - 0.5)*2
        y += (random.random() - 0.5)*2
        z += (random.random() - 0.5)*2

        xArr = np.append(xArr, x)
        yArr = np.append(yArr, y)
        zArr = np.append(zArr, z)
    
    
    outArr = np.array(np.append(outArr,xArr))    
    outArr = np.append(outArr,yArr)
    outArr = np.append(outArr,zArr)

    outArr = np.reshape(outArr,(3,len(xArr)))

    return outArr

np.random.seed(41)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*rwalk3D(100))
ax.plot(*rwalk3D(100))
ax.plot(*rwalk3D(100))
ax.plot(*rwalk3D(100))
plt.show()


def rwalk3Dsum(size, N_trials):
    outArr = []
    def calcSum(size):
        x      = 0.0
        y      = 0.0
        z      = 0.0
        outArr = np.array([0,0,0])
    
        for i in range(0, size+1):
            x += (random.random() - 0.5)*2
            y += (random.random() - 0.5)*2
            z += (random.random() - 0.5)*2
        
        outArr = [x,y,z]
        return outArr
    
    for i in range(0,N_trials):
        outArr.append(calcSum(size))
    print(outArr)
    return outArr


plt.hist(rwalk3Dsum(100, 4000), bins='auto');