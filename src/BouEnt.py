# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 14:09:57 2019

@author: cqt

Criteria for Bound Entanglement ppn >0 and ccn > 1
"""

import numpy as np
from numpy.linalg import eigvalsh,svd

def Realign(Rho,d1,d2):
    Nt = np.size(Rho,axis=0)
    Rrho = np.zeros([Nt,d1**2,d2**2],dtype=np.complex_)
    
    for i in np.arange(0,d1):
        for j in np.arange(0,d1):
            Srho = Rho[:,i*d2:(i+1)*d2,j*d2:(j+1)*d2]
            Rrho[:,j*d1+i,:] = np.reshape(Srho,(Nt,d2*d2))
    return Rrho

def ParTran(Rho,d1,d2):
    Nt = np.size(Rho,axis=0)
    Trho = np.zeros([Nt,d1*d2,d2*d1],dtype=np.complex_)
    Srho = np.zeros([Nt,d2,d2],dtype=np.complex_)
    
    for i in np.arange(0,d1):
        for j in np.arange(0,d1):
            Srho = Rho[:,i*d2:(i+1)*d2,j*d2:(j+1)*d2]
            Trho[:,i*d2:(i+1)*d2,j*d2:(j+1)*d2] = np.transpose(Srho,(0,2,1))
    return Trho

def GetPnCn(rho,d1,d2):
    ppn = np.min(eigvalsh(ParTran(rho,d1,d2)),axis=1)
    ccn = np.sum(np.real(svd(Realign(rho,d1,d2),compute_uv=False)),axis=1)
    return ppn,ccn