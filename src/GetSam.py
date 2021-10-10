# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 17:48:32 2018

@author: cqtlw
"""
import numpy as np
import numpy.linalg as la
import numpy.random as npr

def GetWisSam(Ncolw,Ncols,Nw,Ns,Nt,d,rhop,rhops,rhopso,pomc):
    #each row of corp is a sample
    Nb = Nt-Nw-Ns
    delrho = rhops-rhopso

    CVMin = (la.inv(rhop)+np.eye(d)* (d**2) /(Ncolw-d))
    CVM = np.eye(d)/(la.inv(rhop)+np.eye(d)* d**2 /(Ncolw-d))
#    kra = la.matrix_power(CVM,1/2);      

    CVMins = (la.inv(rhopso)+np.eye(d)* d**2 /(Ncols-d))
    CVMs = np.eye(d)/(la.inv(rhopso)+np.eye(d)* (d**2) /(Ncols-d))
#    kras = la.matrix_power(CVMs,1/2);

    
    Cw = npr.randn(Nw,d,Ncolw)*1j
    Bw = npr.randn(Nw,d,Ncolw)
    Aw = Bw+Cw
    Awh = np.conj(np.transpose(Aw,(0,2,1)))
    rhow = np.matmul(Aw,Awh)
    
    rhosw = np.matmul(rhow,CVM)
    trrw = np.reshape((np.trace(rhosw,offset=0,axis1=1,axis2=2)),(Nw,1,1))
    rhosw = rhosw / trrw
    
    Cb = npr.randn(Nb,d,d)*1j
    Bb = npr.randn(Nb,d,d)
    Ab = Bb+Cb
    Abh = np.conj(np.transpose(Ab,(0,2,1)))
    rhob = np.matmul(Ab,Abh)
    trrb = np.reshape((np.trace(rhob,offset=0,axis1=1,axis2=2)),(Nb,1,1))
    rhosb = rhob / trrb
    
    Cs = npr.randn(Ns,d,Ncols)*1j
    Bs = npr.randn(Ns,d,Ncols)
    As = Bs+Cs
    Ash = np.conj(np.transpose(As,(0,2,1)))
    rhos = np.matmul(As,Ash)
    
    rhoss = np.matmul(rhos,CVMs)
    trrs = np.reshape((np.trace(rhoss,offset=0,axis1=1,axis2=2)),(Ns,1,1))
    rhoss = rhoss / trrs + delrho
    
    rho = np.concatenate((rhosw,rhosb,rhoss),axis=0)
    
    rhor = np.reshape(rho,(Nt,d**2))
    corp = np.real(rhor@pomc)

    return rho,corp,CVMin,CVMins
    #return rho,corp

def UniR(Nt,d):
    Ab = npr.randn(Nt,d,d+1)
    Abh = np.transpose(Ab,(0,2,1))
    rhob = np.matmul(Ab,Abh)
    trrb = np.reshape((np.trace(rhob,offset=0,axis1=1,axis2=2)),(Nt,1,1))
    rhosb = rhob / trrb
    return rhosb

def GetUniSam(Nt,d):
    Cb = npr.randn(Nt,d,d)
    Bb = npr.randn(Nt,d,d)
    Ab = Bb+Cb*1j
    Abh = np.conj(np.transpose(Ab,(0,2,1)))
    rhob = np.matmul(Ab,Abh)
    trrb = np.reshape((np.trace(rhob,offset=0,axis1=1,axis2=2)),(Nt,1,1))
    rhosb = rhob / trrb
    # B = np.reshape(Bb,(Nt,d**2))
    # C = np.reshape(Cb,(Nt,d**2))
    return rhosb

def GetDiriSam(pop,Nt):
    # params = [a1, a2, ..., ak]
    # sample = [random.gammavariate(a,1) for a in params]
    # sample = [v/sum(sample) for v in sample]
    nprob = np.size(pop)
    npop = np.reshape(pop,(nprob,))
    sample = np.array([npr.gamma(a,1.,(Nt,)) for a in npop])
    sample = sample.T
    s = np.sum(sample,axis=1)
    sample = sample / s[:,None]
    return sample