# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 17:59:11 2018

@author: cqt
"""

import numpy.linalg as LA
import numpy as np
from concurrent.futures import ThreadPoolExecutor
# from numba import jit

def getpomc(pom):
    d = np.size(pom,1)
    nprob = np.size(pom,0)
    pomc = np.transpose(np.reshape(np.transpose(pom,(0,2,1)),(d**2,nprob)))
    
    return pomc

def GetRho(corp,pom):
    pomc = getpomc(pom)
    pomci = LA.inv(pomc) # inverse of pomc
    cs = np.shape(corp)
    Nt = int(cs[0])
    d = int(np.sqrt(cs[1]))
    rhor = corp@pomci
    rho = np.reshape(rhor,(Nt,d,d))
    return rho

def GetP(rho,pom):
    pomc = getpomc(pom)
    d = np.size(rho,axis=1)
    Nt = np.size(rho,0)
    rhor = np.reshape(rho,(Nt,d**2))
    corp = np.real(rhor@pomc)
    return corp

def AtoR(Ab,Nt):
    Abh = np.conj(np.transpose(Ab,(0,2,1)))
    rhob = np.matmul(Ab,Abh)
    trrb = np.reshape((np.trace(rhob,offset=0,axis1=1,axis2=2)),(Nt,1,1))
    rhosb = rhob / trrb
    return rhosb

def CBtoR(Nt,d,B,C):
    Bb = np.reshape(B,(Nt,d,d))
    Cb = np.reshape(C,(Nt,d,d))
    Ab = Bb+Cb*1j
    Abh = np.conj(np.transpose(Ab,(0,2,1)))
    rhob = np.matmul(Ab,Abh)
    trrb = np.reshape((np.trace(rhob,offset=0,axis1=1,axis2=2)),(Nt,1,1))
    rhosb = rhob / trrb
    return rhosb

    
def IsPhy(var,meig,pomci,vartype):
    
    if vartype == 'p':
        #cs = np.shape(var)
        #Nt = int(cs[0])
        #d = int(np.sqrt(cs[1]))
        #rhor = var@LA.inv(pomc)
        #rho = np.reshape(rhor,(Nt,d,d))
        rho = GetRho(var,pomci)
    elif vartype == 'r':
        rho = var
        
    ei = LA.eigvalsh(rho)
    # Nt = np.size(rho,axis=0)
    # ei = np.zeros((Nt,1))
    # for i in np.arange(0,Nt):
    #     ei[i] = np.min(LA.eigvalsh(rho[i,:,:]))
    minei = np.min(ei,axis=1)
    isph = minei > meig
    # isph = ei > meig
    
    return isph

# @jit(nopython=False,parallel=True)
def IsPhy_C(rho):
    Nt = np.size(rho,axis=0)
    isph = np.zeros((Nt,))
    for i in np.arange(0,Nt):
        try:
            LA.cholesky(rho[i,:,:])
            isph[i] = 1
        except LA.LinAlgError:
            isph[i] = 0
    return isph

def GetLik(corp,pop):
    logL = np.real((np.sum(np.log(np.abs(corp))*pop,axis=1)))
    return logL
# def IsPhy_C2(rho):
#     try:
#         LA.cholesky(rho)
#         return True
#     except LA.LinAlgError:
#         return False

# def IsPhy_CP(rho):
#     with ThreadPoolExecutor(max_workers=8) as executor:
#         future1 = executor.submit(IsPhy_C2,rho)
#     return future1.result()

