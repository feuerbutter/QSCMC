"""
The Run File for S(C)MC
Created on 01-08-2019
By CqtLW
"""

import numpy as np
import numpy.linalg as LA
from GetSam import Uni
import ptor
import GetPom as gp
import numpy.random as npr
import scipy.io as sio


import time

t1 = time.perf_counter()  #start the clock

#%% User Supplied Quantities

#---Ref sample
Nt = np.int(1e4) # sample size
nq = 1 # no. of qubit
d = 2 ** nq # dimension
nprob = d ** 2 # no. of probability d^2

rhoref = Uni(Nt,d) # A uniform reference sample
# to use the mle as rhoref, we write 

pom = gp.nTetra(nq) # tedrahedron POM
pomc = ptor.getpomc(pom,d,nprob) # put each POM as a column
pomci = LA.inv(pomc) # inverse of pomc
corpref = ptor.GetP(rhoref,pomc,Nt) # matrix of dimension Nt*nprob, storing the p values of the sample (as rows)

# corpref = np.copy(corppy) # to use the last generated sample as the new reference

#---Ref PDF
# It could be either a function of rho and/or probabilities
# Recommend to use log to avoid under/overflow

def refPDF(rho,corp):   
    return 0 # return 0 if using uniform distribution
# if not uniform, we asign values to the corresponding column entries

#---Target PDF
    
# pop = 10*np.ones([1,nprob])
# pop = np.array([[10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6]])
#pop = np.array([[36,13,64,71,14,16,7,15,60,10,84,63,64,9,55,71,8,12,10,16,16,48,67,62,9,64,75,63,10,74,60,73,65,14,62,66,9,57,76,53,82,78,128,22,61,44,25,27,56,12,52,66,14,76,56,78,45,47,22,27,66,68,25,102]])

pop = np.array([[100,200,250,450]]) # sigeye = 2e-2, scale = 1e5

ppop = pop/np.sum(pop) # this is the peak probabilities if the peak is physical
loglpeak = ptor.GetLik(ppop,pop) # log of the likelihood of the peak, using the maximum likelihood estimator as the peak
# if the peak is not physical, use the matlab qse_apg.m to get the state of the MLE, then get its probability coordinates

def tarPDF(rho,corp):
     tarfx = np.real((np.sum(np.log(np.abs(corp))*pop,axis=1))) # the log of the likelihood function for the target posterior
     return tarfx
#    return 0 # this is put 0 because we are aiming to get uniform distribution for a certain region/ states with certain properties

#---Constr Func
scale = 1e5 # the steepness of the step function. the larger it is, the more samples satisfies the constrain
# scale should not be too large that the code breaks down as the points is stuck
# if the no. of samples satisfying the constraint is samll and the code breaks down for larger scale, we need to increase the betak (the no. of distribution steps)



def Constr(rho,corp,beta):
    # for a constraint that is gradually imposed / if a tightening constraint is needed to acheive (VAL>0)    
    #--- constrx1 is used for sampling within a given BLR
#    VAL = ptor.GetLik(corp,pop)-(-1)-loglpeak # for the BLR, VAL = ptor.GetLik(corp,pop)-(log of lambda BLR)-loglpeak
#    constrx1 = np.log((1+np.tanh(scale*beta*VAL))/2) # the Heaviside step function, beta changes for each distribution step

    # if a strict constraint is required
    constrx2 = np.log(ptor.IsPhy_C(rho)) # the physicality constraint
    
    # total constraint
#    constrx = constrx1 + constrx2
    return constrx2

#%% Adjustable Para

#MC Kernel, Stepsize

# Prog() defines how one walks, paranow is the parameters (in a row) for the current point
# sigeye is used to adjust the acceptance rate to be between 0.20 to 0.30, best is 0.234
sigeye = 2e-3 #Gaussian Kernel with Cov. Mat. propto I
def Prog(paranow): #Progate in the Prob. Space using eye
     paranew = np.zeros(np.shape(paranow)) # initialization of the size of the new parameter
     Step = sigeye*np.eye(nprob-1) # the Cov.Mat. used in the Gaussian kernel
     # in prob. space, we do not walk the last coordinate as it is fixed by the other (0:-1) coordinates
     paranew[:,0:-1] = npr.multivariate_normal(np.zeros([1,nprob-1]).flatten(), Step, Nt) # diff. walk use diff. Step
     paranew = paranow + paranew
     paranew[:,-1] = 1 - np.sum(paranew[:,0:-1],axis = 1)
     return paranew

#sigcov = 3.3e-2 #Gaussian Kernel with Cov. Mat propto Cov. Mat of the Current Sample
#def Prog(paranow): #Progate in the Prob. Space using Cov
 #   paranew = np.zeros(np.shape(paranow))
  #  Sig = np.cov(paranow[:,:-1], None, False) # the Cov. Mat. of the current sample
   # Step = sigcov**2 / (nprob-1)**(1/3) * Sig # give the Cov.Mat. for the Gaussian kernal
    #paranew[:,0:-1] = npr.multivariate_normal(np.zeros([1,nprob-1]).flatten(), Step, Nt)
    #paranew = paranow + paranew
    #paranew[:,-1] = 1 - np.sum(paranew[:,0:-1],axis = 1)
    #return paranew


#Distribution Step
betak = 40 # no. of distribution steps
nstep = 15 # no. of MC iteration per distribution step given in stat. papers as 15
nstepl = 20 # no. of MC iteration for the last distribution step. should not be less than 15 to reduced the correlation

#Resampling Threshold
rthres = 4/5 # 1/2 is often used, the larger it is the more often the resampling is done. 
# here we want it larger to reduce correlation. if rthres is too large, it takes longer CPU time.


#%% Running of Algorithm

# define the intermediate distributions, in the log scale
def intPDF(rho,corp,beta):
    intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + Constr(rho,corp,beta)
    return intx

Nthres = Nt * rthres # the resampling threshold
betai = 0 # initial distribution is refPDF
delbeta = 1 / betak
rhonow = np.copy(rhoref)
corpnow = np.copy(corpref)
# use np.copy() instead of assigning the pointers using "=" directly to keep the original rhoref/corpref
ess = Nt # initiate the effective sample size Nt
w = np.ones([Nt,]) / Nt # initiate the weight to be 1/Nt for each sample point

for betakk in np.arange(1,betak+1,1): # in matlab this would be for betakk=2:betak
    beta = betai + betakk * delbeta
    print(['beta is', beta])
    intgx = intPDF(rhonow,corpnow,beta) # the intermediate distribution
    if betakk == betak: # the last distribution step
        nstep = nstepl
    # the MC runs for each distribution step
    for nstepk in np.arange(0,nstep):
        if nstepk == 0:
            ogx = intPDF(rhonow,corpnow,beta-delbeta) # the current reference
            w = w * np.exp(intgx-ogx) # update weight
            w = w / np.sum(w) 
            ess = 1 / np.sum(w**2) # if this is very close to the sample size, the sample is good and no need resampling
            print(['ess is',ess])
            
            if (ess < Nthres) + (betakk == betak): # at the last step, if ess<Nthres we do reseampling
                print('resam \n')
                rsind = npr.choice(Nt,Nt,p=w) # resampling index, this line automaticly does the resampling
                rhonow = rhonow[rsind,:,:]
                corpnow = corpnow[rsind,:]
                intgx = intgx[rsind] # the distribution probability aftet he resampling, without reevaluating the values
                # reassign the weight
                w = 1 / Nt *np.ones([Nt,])
                ess = 1 / np.sum(w**2)
                
        corpnew = Prog(corpnow) # MC propagating in probability space
        rhonew = ptor.GetRho(corpnew,pomci)
        intfx = intPDF(rhonew,corpnew,beta) # the desity after MC iteraction

        #accept-reject
        al = intfx - intgx
        draw = np.log(npr.rand(np.size(al)))
        acc = draw < al # this is the standard accept and rejection
        # if we only accept with certain condition, we modify the draw accordingly
        rhonow[acc,:,:] = rhonew[acc,:,:]
        corpnow[acc,:] = corpnew[acc,:]
        intgx[acc] = intfx[acc]
        accr = np.sum(acc) / Nt
    print(accr)
    print(['time is', time.perf_counter() - t1]) # desplay elapsed time

# the variables to be stored and its file name
phnow =   (Constr(rhonow,corpnow,1)==0)
rhopy = rhonow[phnow,:,:]
corppy = corpnow[phnow,:]
sio.savemat('pycorp1q1e4S.mat', mdict={'corppy':corppy,'rhopy':rhopy})
