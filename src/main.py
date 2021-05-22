"""
The Run File for QS(C)MC
Created on 19-04-2021
By Weijun Li
"""

import QSCMC as Q
import numpy as np
import time 
import scipy.io as sio

t1 = time.perf_counter()  #start the clock

#----Running Parameters (sam_type independent)
N_total = int(1e5) # no. of points in the sample
betak = 5
step_size = 1 #Gaussian Kernel with Cov. Mat propto Cov. Mat of the Current Sample
n_step = 15


# Specifying the type of target sample, 
#----sam_type specific paras--------------------

#----Target Posterior
# 'posterior' for a posterior distribution with the data counts,
# or set to anything else for Bound-entanglement(BE) or bounded likelihood region
# sam_type = 'posterior' 
# constr_type = 'phys'
n_qubit = 2 # no. of qubit
pop = 10 * np.array([[10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6]]) # data counts for the posterior
# rho_acc,corp_acc = Q.genTarSam(N_total,n_qubit,sam_type,constr_type,pop,betak,n_step,step_size)
rho_acc,corp_acc = Q.genTarSamPos(N_total,n_qubit,pop,betak,n_step,step_size)

#----Bounded likelihood region
sam_type = ''
constr_type = 'BLR'
scale_BLR = 1e5

pop = 10 * np.array([[10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6]]) # data counts for the posterior
logL_peak = # supplied by user, or obtained through superfast qse in MATLAB


#----Bound Entanglement

d1 = 3
d2 = 3

scale_ppn = 1e3
scale_cnn = 1e3





print(['time is', time.perf_counter() - t1])
sio.savemat('pycorp2q1e5c1e3bk5.mat', mdict={'corppy':corp_acc,'rhopy':rho_acc})
