# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 12:06:07 2018
\n pdf for Wischart \n
@author: cqt
"""

import numpy as np
import numpy.linalg as LA
import math

def UnnorCvmDen(rho,Ncol,d,CVMin):   
    lg = np.real((Ncol-d)*np.log(LA.det(rho))-(Ncol*d)*
                 np.log(np.trace(rho@CVMin,offset=0,axis1=1,axis2=2)))
    return lg
    
def GetLogNorm(Ncol,d,CVMin):
    gl = np.zeros([d,1])
    for i in np.arange(d):
        gl[i] = math.lgamma(Ncol-i)
        
    lnor = math.lgamma(Ncol*d) - d*(d-1)/2*np.log(np.pi) - np.sum(gl) \
    + Ncol * np.log(np.real(LA.det(CVMin)))
    return lnor
    