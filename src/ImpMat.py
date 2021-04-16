# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 17:36:45 2018

@author: cqt
"""

import numpy as np
import scipy.io as sio
import h5py

def v73(filename,varname):
    f = h5py.File(filename,'r')
    var = f.get(varname)
    var = np.array(var)
    return var



def vo(filename,varname,ifstkm):
    
    f = sio.loadmat(filename)
    var = f[varname]
    if ifstkm == 1:
        var = np.moveaxis(var,2,0)
    return var