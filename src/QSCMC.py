"""
Major functions for QSCMC
"""

import numpy as np
import numpy.linalg as LA
from GetSam import GetUniSam
import ptor
import GetPom as gp
import numpy.random as npr
import scipy.io as sio
import BouEnt as be

import BLR1q as BLR

def refPDF(rho,corp):   
    return 0 # return 0 if using uniform distribution
    # if not uniform, we asign values to the corresponding column entries

def Prog(paranow,step_size,kern='C'): #Propagate in the Prob. Space using Cov
    n_para = np.size(paranow,axis=1)
    N_total = np.size(paranow,axis=0)

    if kern == 'I':
        Step = step_size * np.eye(n_para-1)
    elif kern == 'C':
        Sig = np.cov(paranow[:,:-1], None, False)
        Step = step_size ** 2 / (n_para-1)**(1/3) * Sig

    paranew = np.zeros(np.shape(paranow))
    paranew[:,0:-1] = npr.multivariate_normal(np.zeros([1,n_para-1]).flatten(), Step, N_total)
    paranew = paranow + paranew
    paranew[:,-1] = 1 - np.sum(paranew[:,0:-1],axis=1)
    return paranew

def dstrProg(betak,n_step,N_total,rhonow,corpnow,pom,step_size,intPDF,kern='C'):
    """This does SMC Propagation in p space, tarfunc(beta,corp)"""

    n_step_last = n_step + 5

    #Resampling Threshold
    rthres = 4/5 # 1/2 is often used, the larger it is the more often the resampling is done. 
    # here we want it larger to reduce correlation. if rthres is too large, it takes longer CPU time.
    Nthres = N_total * rthres # the resampling threshold

    delbeta = 1 / betak
    ess = N_total # initiate the effective sample size N_total
    w = np.ones([N_total,]) / N_total # initiate the weight to be 1/N_total for each sample point

    for beta in np.linspace(delbeta,1,betak):
        intgx = intPDF(rhonow,corpnow,beta) # the intermediate distribution
        print('beta=',beta)
        if beta == 1: # the last distribution step
            n_step = n_step_last
        # the MC runs for each distribution step
        for n_step_idx in np.arange(0,n_step):
            if n_step_idx == 0:
                ogx = intPDF(rhonow,corpnow,beta-delbeta) # the current reference
                w = w * np.exp(intgx-ogx) # update weight
                w = w / np.sum(w)
                ess = 1 / np.sum(w**2) # if this is very close to the sample size, the sample is good and no need resampling
                print(['ess is',np.round(ess)])

                if (ess < Nthres) + (beta == 1): # at the last step, if ess<Nthres we do reseampling
                    print('resam \n')
                    rsind = npr.choice(N_total,N_total,p=w) # resampling index, this line automaticly does the resampling
                    rhonow = rhonow[rsind,:,:]
                    corpnow = corpnow[rsind,:]
                    intgx = intgx[rsind] # the distribution probability aftet he resampling, without reevaluating the values
                    # reassign the weight
                    w = 1 / N_total *np.ones([N_total,])
                    ess = 1 / np.sum(w**2)

            corpnew = Prog(corpnow,step_size,kern) # MC propagating in probability space
            rhonew = ptor.GetRho(corpnew,pom)
            intfx = intPDF(rhonew,corpnew,beta) # the desity after MC iteraction

            #accept-reject
            al = intfx - intgx
            draw = np.log(npr.rand(np.size(al)))
            acc = draw < al # this is the standard accept and rejection
            # if we only accept with certain condition, we modify the draw accordingly
            rhonow[acc,:,:] = rhonew[acc,:,:]
            corpnow[acc,:] = corpnew[acc,:]
            intgx[acc] = intfx[acc]
            accr = np.sum(acc) / N_total
        print('betaaccr=',accr)
    return rhonow,corpnow
    

def genTarSamPos(N_total,n_qubit,pop,betak,n_step,step_size):
    n_step_last = n_step + 5

    def tarPDF(rho,corp):
        tarfx = np.real((np.sum(np.log(np.abs(corp))*pop,axis=1))) # the log of the likelihood function for the target posterior
        return tarfx
    def constr(rho):
        constr_val = np.log(ptor.IsPhy_C(rho)) 
        return constr_val
    def intPDF(rho,corp,beta):
        intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constr(rho)
        return intx
    d = 2 ** n_qubit # dimension
    # n_prob = d ** 2 # no. of probability d^2

    rhoref = GetUniSam(N_total,d) # A uniform reference sample
    # to use the mle as rhoref, we would need to get the rhomle first and make N_total copies
    pom = gp.nTetra(n_qubit) # tedrahedron POM
    corpref = ptor.GetP(rhoref,pom) # matrix of dimension N_total*n_prob, storing the p values of the sample (as rows)

    # initialize the propgation
    rhonow = np.copy(rhoref)
    corpnow = np.copy(corpref)

    # run betak distribution propagation steps
    rhonow,corpnow = dstrProg(betak,n_step,N_total,rhonow,corpnow,pom,step_size,intPDF)

    return rhonow,corpnow

def genTarSamBE(rho_ref,N_total,d1,d2,betak,n_step,step_size,scale_ppn,scale_cnn,kern='I'):
    n_step_last = n_step + 5

    def tarPDF(rho,corp):
        return 0
    def constr(rho,beta):
        ppn,ccn = be.GetPnCn(rho,d1,d2)
        # an approximation of the Heaviside function using tanh
        constr_val = np.log((1+np.tanh(scale_ppn*beta*ppn))/2) + np.log((1+np.tanh(scale_cnn*beta*(ccn-1)))/2)
        constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) # The BE states should of course be physical 

        # for 2x4 subsequent MC in either the ppn or ccn direction
        # constr_val = np.log((1+np.tanh(scale_cnn*beta*(ccn-1)))/2)
        # constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) + np.log(ppn>0) 
        # constr_val = ppn
        # constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) + np.log(ccn>1) # The BE states should of course be physical 

        return constr_val
    def intPDF(rho,corp,beta):
        intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constr(rho,beta)
        return intx
    
    d = d1 * d2
    n_prob = d ** 2

    # create POM for the composite system, the default is using square root pom, user can change it when suitable
    if (d1 & (d1-1) ==0) and d1 !=0:
        n_qubit_1 = int(np.log2(d1))
        pom1 = gp.nTetra(n_qubit_1)

    elif d1 == 3:
        pom1 = gp.QtsPom()
    else:
        if isinstance(d1,int):
            print("Using the sqrtPom in the GetPom module.")
            pom1 = gp.SqrtPom(d1)

    if (d2 & (d2-1) ==0) and d2 !=0:
        n_qubit_1 = int(np.log2(d2))
        pom2 = gp.nTetra(n_qubit_1)

    elif d2 == 3:
        pom2 = gp.QtsPom()
    else:
        if isinstance(d2,int):
            print("Using the sqrtPom in the GetPom module.")
            pom2 = gp.SqrtPom(d2)


    pom = np.zeros((n_prob,d,d),dtype=np.complex_)

    for idx_1 in np.arange(0,d1**2):
        for idx_2 in np.arange(0,d2**2):
            pom[idx_1*d2**2+idx_2,:,:] = np.kron(pom1[idx_1,:,:],pom2[idx_2,:,:])

    print('POM obtained')

    # rhoref = GetUniSam(N_total,d) # A uniform reference sample
    # # to use the mle as rhoref, we would need to get the rhomle first and make N_total copies
    corpref = ptor.GetP(rho_ref,pom) # matrix of dimension N_total*n_prob, storing the p values of the sample (as rows)

    ppn,ccn = be.GetPnCn(rho_ref,d1,d2)
    beno = (ppn>0) & (ccn>1)
    print('original no. of ppn is ', np.sum(ppn>0), '\noriginal no. of ccn is ', np.sum(ccn>1), '\noriginal no. of BE is ', np.sum(beno))

    rhonow = np.copy(rho_ref)
    corpnow = np.copy(corpref)
    rhonow,corpnow = dstrProg(betak,n_step,N_total,rhonow,corpnow,pom,step_size,intPDF,kern)
    
    return rhonow,corpnow

def genTarSamBLR(N_total,n_qubit,pop,betak,n_step,step_size,rho_mle,scale_BLR,log_lambda,kern='I'):
    n_step_last = n_step + 5

    pom = gp.nTetra(n_qubit) # tedrahedron POM
    corp_mle = ptor.GetP(rho_mle,pom)
    logL_peak = ptor.GetLik(corp_mle,pop)

    def tarPDF(rho,corp):
        return 0
    # the constraint function for imposing bounded likelihood
    def constr(rho,corp,beta):
        VAL = ptor.GetLik(corp,pop)-(-1)-logL_peak-log_lambda # for the BLR, VAL = ptor.GetLik(corp,pop)-(log of lambda BLR)-loglpeak
        # an approximation of the Heaviside function using tanh
        constr_val = np.log((1+np.tanh(scale_BLR*beta*VAL))/2) 
        constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) # make sure the states remain physical
        return constr_val
    def intPDF(rho,corp,beta):
        intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constr(rho,corp,beta)
        return intx
    d = 2 ** n_qubit # dimension
    # n_prob = d ** 2 # no. of probability d^2

    rhoref = GetUniSam(N_total,d) # A uniform reference sample
    # to use the mle as rhoref, we would need to get the rhomle first and make N_total copies
    corpref = ptor.GetP(rhoref,pom) # matrix of dimension N_total*n_prob, storing the p values of the sample (as rows)

    # initialize the propgation
    rhonow = np.copy(rhoref)
    corpnow = np.copy(corpref)

    # run betak distribution propagation steps
    rhonow,corpnow = dstrProg(betak,n_step,N_total,rhonow,corpnow,pom,step_size,intPDF,kern)

    return rhonow,corpnow

def dstrProg_1q(betak,n_step,N_total,coordnow,corpnow,step_size,intPDF,kern='C'):
    """This does SMC Propagation in p space, tarfunc(beta,corp)"""

    n_step_last = n_step + 5

    #Resampling Threshold
    rthres = 4/5 # 1/2 is often used, the larger it is the more often the resampling is done. 
    # here we want it larger to reduce correlation. if rthres is too large, it takes longer CPU time.
    Nthres = N_total * rthres # the resampling threshold

    delbeta = 1 / betak
    ess = N_total # initiate the effective sample size N_total
    w = np.ones([N_total,]) / N_total # initiate the weight to be 1/N_total for each sample point

    for beta in np.linspace(delbeta,1,betak):
        intgx = intPDF(coordnow,corpnow,beta) # the intermediate distribution
        print('beta=',beta)
        if beta == 1: # the last distribution step
            n_step = n_step_last
        # the MC runs for each distribution step
        for n_step_idx in np.arange(0,n_step):
            if n_step_idx == 0:
                ogx = intPDF(coordnow,corpnow,beta-delbeta) # the current reference
                w = w * np.exp(intgx-ogx) # update weight
                w = w / np.sum(w)
                ess = 1 / np.sum(w**2) # if this is very close to the sample size, the sample is good and no need resampling
                print(['ess is',np.round(ess)])

                if (ess < Nthres) + (beta == 1): # at the last step, if ess<Nthres we do reseampling
                    print('resam \n')
                    rsind = npr.choice(N_total,N_total,p=w) # resampling index, this line automaticly does the resampling
                    coordnow = coordnow[rsind,:]
                    corpnow = corpnow[rsind,:]
                    intgx = intgx[rsind] # the distribution probability aftet he resampling, without reevaluating the values
                    # reassign the weight
                    w = 1 / N_total *np.ones([N_total,])
                    ess = 1 / np.sum(w**2)

            corpnew = Prog(corpnow,step_size,kern) # MC propagating in probability space
            coordnew = BLR.probToCoord(corpnew,N_total)
            intfx = intPDF(coordnew,corpnew,beta) # the desity after MC iteraction

            #accept-reject
            al = intfx - intgx
            draw = np.log(npr.rand(np.size(al)))
            acc = draw < al # this is the standard accept and rejection
            # if we only accept with certain condition, we modify the draw accordingly
            coordnow[acc,:] = coordnew[acc,:]
            corpnow[acc,:] = corpnew[acc,:]
            intgx[acc] = intfx[acc]
            accr = np.sum(acc) / N_total
        print('betaaccr=',accr)
    return coordnow,corpnow



def genTarSamBLR_1q(coordref,corpref,N_total,n_qubit,pop,betak,n_step,step_size,rho_mle,scale_BLR,logL_peak,log_lambda,kern='I'):
    n_step_last = n_step + 5

    def tarPDF(coord,corp):
        return 0
    # the constraint function for imposing bounded likelihood
    def constr(coord,corp,beta):
        VAL = ptor.GetLik(corp,pop)-logL_peak-log_lambda # for the BLR, VAL = ptor.GetLik(corp,pop)-(log of lambda BLR)-loglpeak
        # an approximation of the Heaviside function using tanh
        constr_val = np.log((1+np.tanh(scale_BLR*beta*VAL))/2) 
        phys = (coord[:,0]**2+coord[:,1]**2) <= 1
        constr_val = constr_val + np.log(phys) # make sure the states remain physical
        return constr_val
    def intPDF(coord,corp,beta):
        intx = tarPDF(coord,corp) * (beta) + refPDF(coord,corp) * (1-beta) + constr(coord,corp,beta)
        return intx
    d = 2 ** n_qubit # dimension
    # n_prob = d ** 2 # no. of probability d^2

 
    # initialize the propgation
    coordnow = np.copy(coordref)
    corpnow = np.copy(corpref)

    # run betak distribution propagation steps
    coordnow,corpnow = dstrProg_1q(betak,n_step,N_total,coordnow,corpnow,step_size,intPDF,kern)

    return coordnow,corpnow
# def genTarSam(N_total,n_qubit,sam_type,constr_type,pop,betak,n_step,step_size):

#     n_step_last = n_step + 5
#     def refPDF(rho,corp):   
#         return 0 # return 0 if using uniform distribution
#         # if not uniform, we asign values to the corresponding column entries
#     def tarPDF(rho,corp):
#         if sam_type == 'posterior':
#             tarfx = np.real((np.sum(np.log(np.abs(corp))*pop,axis=1))) # the log of the likelihood function for the target posterior
#         else:
#             tarfx = 0
#         return tarfx


    # # the constraint function for imposing bounded likelihood
    # def constrBLR(rho,corp,beta,scale_BLR,logL_peak):
    #     VAL = ptor.GetLik(corp,pop)-(-1)-logL_peak # for the BLR, VAL = ptor.GetLik(corp,pop)-(log of lambda BLR)-loglpeak
    #     constr_val = np.log((1+np.tanh(scale_BLR*beta*VAL))/2) 
    #     return constr_val

#     def constrBE(rho,corp,beta,d1,d2,scale_ppn,scale_ccn):
#         ppn,ccn = be.GetPnCn(rho,d1,d2)
#         constr_val = np.log((1+np.tanh(scalep*beta*ppn))/2) + np.log((1+np.tanh(scale*beta*(ccn-1)))/2)
#         constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) # The BE states should of course be physical 
#         return constr_val

#     def constrPhys(rho):
#         constr_val = np.log(ptor.IsPhy_C(rho)) 
#         return constr_val

#     # def Constr(rho,corp,beta,scale_BLR=0,scalep=0):
#     #     # for a constraint that is gradually imposed / if a tightening constraint is needed to acheive (VAL>0)    
#     #     if constr_type == 'BLR':
#     #         #--- constrx1 is used for sampling within a given BLR
#     #         VAL = ptor.GetLik(corp,pop)-(-1)-loglpeak # for the BLR, VAL = ptor.GetLik(corp,pop)-(log of lambda BLR)-loglpeak
#     #         constr_val = np.log((1+np.tanh(scale*beta*VAL))/2) 
#     #     elif constr_type == 'BE':
#     #         ppn,ccn = be.GetPnCn(rho,d1,d2)
#     #         constr_val = np.log((1+np.tanh(scalep*beta*ppn))/2) + np.log((1+np.tanh(scale*beta*(ccn-1)))/2)
#     #         constr_val = constr_val + np.log(ptor.IsPhy_C(rho)) # The BE states should of course be physical 
#     # #    constrx1 = np.log((1+np.tanh(scale*beta*VAL))/2) # the Heaviside step function, beta changes for each distribution step
#     #     # the physicality constraint
#     #     elif constr_type == 'physical':
#     #         # if a strict constraint is required
#     #         constr_val = np.log(ptor.IsPhy_C(rho)) 
        
#     #     # total constraint
#     # #    constrx = constrx1 + constrx2
#     #     return constr_val

#     def intPDF(rho,corp,beta):
#         if constr_type == 'BLR':
#             intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constrBLR(rho,corp,beta,scale_BLR,logL_peak)
#         elif constr_type == 'BE':
#             intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constrBE(rho,corp,beta,d1,d2,scale_ppn,scale_ccn)
#         elif constr_type = 'phys':
#             intx = tarPDF(rho,corp) * (beta) + refPDF(rho,corp) * (1-beta) + constrPhys(rho)
#         return intx


#     def Prog(paranow): #Propagate in the Prob. Space using Cov
#         paranew = np.zeros(np.shape(paranow))
#         Sig = np.cov(paranow[:,:-1], None, False)
#         Step = step_size ** 2 / (n_prob-1)**(1/3) * Sig
#         paranew[:,0:-1] = npr.multivariate_normal(np.zeros([1,n_prob-1]).flatten(), Step, N_total)
#         paranew = paranow + paranew
#         paranew[:,-1] = 1 - np.sum(paranew[:,0:-1],axis = 1)
#         return paranew

#     d = 2 ** n_qubit # dimension
#     n_prob = d ** 2 # no. of probability d^2

#     rhoref = GetUniSam(N_total,d) # A uniform reference sample
#     # to use the mle as rhoref, we would need to get the rhomle first and make N_total copies

#     pom = gp.nTetra(n_qubit) # tedrahedron POM
#     corpref = ptor.GetP(rhoref,pom) # matrix of dimension N_total*n_prob, storing the p values of the sample (as rows)
    
#     #Resampling Threshold
#     rthres = 4/5 # 1/2 is often used, the larger it is the more often the resampling is done. 
#     # here we want it larger to reduce correlation. if rthres is too large, it takes longer CPU time.

#     Nthres = N_total * rthres # the resampling threshold
#     betai = 0 # initial distribution is refPDF
#     delbeta = 1 / betak
#     rhonow = np.copy(rhoref)
#     corpnow = np.copy(corpref)
#     # use np.copy() instead of assigning the pointers using "=" directly to keep the original rhoref/corpref
#     ess = N_total # initiate the effective sample size N_total
#     w = np.ones([N_total,]) / N_total # initiate the weight to be 1/N_total for each sample point

#     for betakk in np.arange(1,betak+1,1): # in matlab this would be for betakk=2:betak
#         beta = betai + betakk * delbeta
#         print(['beta is', beta])
#         intgx = intPDF(rhonow,corpnow,beta) # the intermediate distribution
#         if betakk == betak: # the last distribution step
#             n_step = n_step_last
#         # the MC runs for each distribution step
#         for n_step_idx in np.arange(0,n_step):
#             if n_step_idx == 0:
#                 ogx = intPDF(rhonow,corpnow,beta-delbeta) # the current reference
#                 w = w * np.exp(intgx-ogx) # update weight
#                 w = w / np.sum(w) 
#                 ess = 1 / np.sum(w**2) # if this is very close to the sample size, the sample is good and no need resampling
#                 print(['ess is',ess])
                
#                 if (ess < Nthres) + (betakk == betak): # at the last step, if ess<Nthres we do reseampling
#                     print('resam \n')
#                     rsind = npr.choice(N_total,N_total,p=w) # resampling index, this line automaticly does the resampling
#                     rhonow = rhonow[rsind,:,:]
#                     corpnow = corpnow[rsind,:]
#                     intgx = intgx[rsind] # the distribution probability aftet he resampling, without reevaluating the values
#                     # reassign the weight
#                     w = 1 / N_total *np.ones([N_total,])
#                     ess = 1 / np.sum(w**2)
                    
#             corpnew = Prog(corpnow) # MC propagating in probability space
#             rhonew = ptor.GetRho(corpnew,pom)
#             intfx = intPDF(rhonew,corpnew,beta) # the desity after MC iteraction

#             #accept-reject
#             al = intfx - intgx
#             draw = np.log(npr.rand(np.size(al)))
#             acc = draw < al # this is the standard accept and rejection
#             # if we only accept with certain condition, we modify the draw accordingly
#             rhonow[acc,:,:] = rhonew[acc,:,:]
#             corpnow[acc,:] = corpnew[acc,:]
#             intgx[acc] = intfx[acc]
#             accr = np.sum(acc) / N_total
#         print(['ar is', accr])
#         # print(['time is', time.perf_counter() - t1]) # desplay elapsed time
#     phnow = (Constr(rhonow,corpnow,1)==0)
#     rhopy = rhonow[phnow,:,:]
#     corppy = corpnow[phnow,:]
#     # sio.savemat('pycorp2q1e5c1e3bk10(newmain).mat', mdict={'corppy':corppy,'rhopy':rhopy})
#     return rhonow,corpnow




