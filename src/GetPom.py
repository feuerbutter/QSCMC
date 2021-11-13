# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 16:51:03 2018

@author: cqtlw
"""



import numpy as np
from numpy import sqrt
from numpy.random import randn
from numpy.linalg import eigvalsh
from numpy.linalg import inv
from numpy.random import rand
from scipy.linalg import sqrtm

def Trine():
    Tr = np.zeros((3,2,2),dtype=complex)
    Tr[0,:,:] = np.sqrt(2/3) * [[1.,0.],[0.,0.]]
    Tr[1,:,:] = 1/4 * [[np.sqrt(2/3),-np.sqrt(2)],[-np.sqrt(2),np.sqrt(6)]]
    Tr[2,:,:] = 1/4 * [[np.sqrt(2/3),np.sqrt(2)],[np.sqrt(2),np.sqrt(6)]]
    return Tr

def DisTri(phi):
    phi2 = np.pi - phi
    phi3 = np.pi + phi
    Tr = np.zeros((3,2,2),dtype=complex)
    Tr[0,:,:] = [[1.,1.],[1.,1.]]
    Tr[1,:,:] = [[1,np.exp(-1j*phi2)],[np.exp(1j*phi2),1]]
    Tr[2,:,:] = [[1,np.exp(-1j*phi3)],[np.exp(1j*phi3),1]]
    return Tr

def Pauli():
    S = np.zeros((4,2,2),dtype=complex)
    S[0,:,:] = np.eye(2)
    S[1,:,:] = [[0.,1.],[1.,0.]]
    S[2,:,:] = [[0.,-1j],[1j,0.]]
    S[3,:,:] = [[1.,0.],[0.,-1.]]
    return S

S = Pauli()

def Tetra():
#    tvec1 = np.array([1.,1./sqrt(3.),-1./sqrt(3.),-1./sqrt(3.)])/4
#    tvec1 = np.array([1.,-1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)])/4
#    tvec1 = np.array([1.,-1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)])/4
#    tvec1 = np.array([1.,1./sqrt(3.),-1./sqrt(3.),-1./sqrt(3.)])/4
    
    tvec = np.array([[1.,1./sqrt(3.),-1./sqrt(3.),-1./sqrt(3.)],
                      [1.,-1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)],
                      [1.,-1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)],
                      [1.,1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)]])/4
#    S = GetPom.Pauli()
    T = np.zeros((4,2,2),dtype=complex)
    for i in np.arange(4):
        tempv = np.reshape(tvec[i,:],(4,1,1))
        T[i,:,:] = np.sum(tempv*S,0)
        
    return T

T = Tetra()       

def nTetra(n):
    pom = T
    for i in np.arange(n-1):
        tempn = 2 ** (2*(i+1))
        d = 2 ** (i+2)
        tempom = np.zeros((4*tempn,d,d),dtype=complex)
        for j in np.arange(4):
            for k in np.arange(tempn):              
                tempom[tempn*(j)+k,:,:] = np.kron(pom[k,:,:],T[j,:,:])
                
        pom = tempom
        
    return pom
        
def RandPom(d,nprob,scale):
    ex = -1
    pom = np.zeros((nprob,d,d),dtype=np.complex_)
    
    while ex < 0:
        AR = randn(nprob-1,d,d)
        AI = randn(nprob-1,d,d)
        A = AR + AI*1j
        Ah = np.conj(np.transpose(A,(0,2,1)))
        pom[0:-1,:,:] = scale * np.matmul(A,Ah)
        pom[-1,:,:] = np.eye(d) - np.sum(pom[0:-1,:,:],axis=0)

        ex = np.min(eigvalsh(pom[-1,:,:]))
    
    return pom

def QtsPom():
    pom = np.zeros((9,3,3),dtype=np.complex_)
    ga = np.exp(1j*2*np.pi/3)
    ga2 = ga ** 2

    psi = np.array([[0., 1., 1.], [0., ga, ga2], [0., ga2, ga], [1., 0., 1.], [ga2, 0., ga], [ga, 0., ga2], [1., 1., 0.], [ga, ga2, 0.], [ga2, ga, 0.]])
    psi = psi.T

    for i in np.arange(0,9):
        ps = np.reshape(psi[:,i],(3,1))
        pom[i,:,:] = np.kron(ps,ps.conj().T) / 6
        

    return pom

def SqrtPom(d):
    nprob = d**2
    a = rand(nprob,d,1)
    ah = np.conj(np.transpose(a,(0,2,1)))
    ao = np.zeros((nprob,d,d),dtype=complex)
    for i in np.arange(0,nprob):
        ao[i,:,:] = np.kron(a[i,:,:],ah[i,:,:])

    A = np.sum(ao,axis=0)
    Asi = sqrtm(inv(A))
    pom = Asi @ ao @ Asi

    return pom