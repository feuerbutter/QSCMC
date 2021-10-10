"""
Utility functions for BLR 1 qubit example
Created on 29-06-2021
By CqtLW
"""

import numpy as np

cgs = 0.1327
cg = np.sqrt(cgs)
sgs = 1 - cgs

def coordToProb(coord_xy,N_total):
    # convert coordinate (x,y) in the bloch circle
    # not true, it is actually rhos
    p = np.zeros((N_total,3))
    p[:,0] = 1/2 * sgs * (1+coord_xy[:,0]) 
    p[:,1] = 1/4 * (1+cgs-sgs*coord_xy[:,0]+2*cg*coord_xy[:,1])
    p[:,2] = 1/4 * (1+cgs-sgs*coord_xy[:,0]-2*cg*coord_xy[:,1])
    return p

def probToCoord(p,N_total):
    coord_xy = np.zeros((N_total,2))
    coord_xy[:,0] = 2 / sgs * p[:,0] - 1
    coord_xy[:,1] = (p[:,1]-p[:,2]) / cg 
    return coord_xy