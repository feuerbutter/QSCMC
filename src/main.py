"""
The Run File for QS(C)MC
Created on 19-04-2021
By Weijun Li
"""
#%%
import QSCMC as Q
import numpy as np
import time 
import scipy.io as sio
import BouEnt as be
import ptor as ptor
import GetSam as gs
from os.path import dirname, join as pjoin

t1 = time.perf_counter()  #start the clock

# setting the directory for the data folder
data_dir = dirname('../data/')

#----Running Parameters (sam_type independent)
N_total = int(1e3) # no. of points in the sample
n_step = 15


# Specifying the type of target sample, 
#----sam_type specific paras--------------------

# #----Target Posterior
# # 'posterior' for a posterior distribution with the data counts,
# # or set to anything else for Bound-entanglement(BE) or bounded likelihood region
# # sam_type = 'posterior' 
# # constr_type = 'phys'

#--Sample paras for 2qb pop with 1000 counts
# betak = 10
# step_size = 1.

# n_qubit = 2 # no. of qubit
# pop = 10 * np.array([[10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6]]) # data counts for the posterior
# # rho_acc,corp_acc = Q.genTarSam(N_total,n_qubit,sam_type,constr_type,pop,betak,n_step,step_size)
# rho_acc,corp_acc = Q.genTarSamPos(N_total,n_qubit,pop,betak,n_step,step_size)


#----Bounded Likelihood Region---------------------

# betak = 1000
# step_size = 1e-8

# # user needs to put the rho_mle in the data folder and the following line imports it
# mat_name = 'mle2qb.mat' # name of the mat data file that stores the rho_mle
# var_name = 'rhomle' # the variable name that rho_mle is stored as
# rho_mle = sio.loadmat(pjoin(data_dir,mat_name))[var_name]

# n_qubit = 2 # no. of qubit
# pop = np.array([[10, 4, 6, 4, 7, 6, 5, 6, 5, 6, 10, 6, 5, 6, 8, 6]]) # data counts for the posterior

# scale_BLR = 1e5

# rho_acc,corp_acc = Q.genTarSamBLR(N_total,n_qubit,pop,betak,n_step,step_size,rho_mle,scale_BLR)




#----Bound Entanglement---------------------------------

#--Sample paras for 3x3
betak = 20
step_size = 3e-9 
d1 = 3
d2 = 3
d = d1 * d2
scale_ppn = 1e4
scale_cnn = 3e3


#--Sample paras for 2x4
# betak = int(20)
# step_size = 5e-7
# d1 = 2
# d2 = 4
# d = d1 * d2
# scale_ppn = 2e2
# scale_cnn = 8e2

#——————————————————————————Start from an existing reference sample
#mat_name = 'BE2x4ppnonly.mat' # name of the mat data file that stores the starting reference sample 
#var_name = 'rhopy' # the variable name that rho_mle is stored as
#rho_ref = sio.loadmat(pjoin(data_dir,mat_name))[var_name]
#ppn_ref,ccn_ref = be.GetPnCn(rho_ref,d1,d2)

#—————————for picking out the ccn>1 states
# ccn_ind = (ccn_ref>1)
# n_ccn = np.sum(ccn_ind)
# print(n_ccn)
# n_not_ccn = N_total - n_ccn
# rho_ccn = np.zeros(np.shape(rho_ref),dtype=np.complex_)
# rho_ccn[:n_ccn,:,:] = rho_ref[ccn_ind,:,:]
# rho_ccn[n_ccn:,:,:] = rho_ccn[0:n_not_ccn,:,:]
# rho_ref = rho_ccn

#—————————for picking out the ppn>0 states
# ppn_ind = (ppn_ref>0)
# n_ppn = np.sum(ppn_ind)
# print(n_ppn)
# n_not_ppn = N_total - n_ppn
# rho_ppn = np.zeros(np.shape(rho_ref),dtype=np.complex_)
# rho_ppn[:n_ppn,:,:] = rho_ref[ppn_ind,:,:]
# rho_ppn[n_ppn:,:,:] = rho_ppn[0:n_not_ppn,:,:]
# rho_ref = rho_ppn

#——————————————————————————Start from a uniform sample
rho_ref = gs.GetUniSam(N_total,d) # A uniform reference sample
# to use the mle as rhoref, we would need to get the rhomle first and make N_total copies

#%%
rho_acc,corp_acc = Q.genTarSamBE(rho_ref,N_total,d1,d2,betak,n_step,step_size,scale_ppn,scale_cnn)

ppn_acc,ccn_acc = be.GetPnCn(rho_acc,d1,d2)
be_acc = (ppn_acc>0) & (ccn_acc>1)
n_be = np.sum(be_acc)

print('no. of ppn is ', np.sum(ppn_acc>0), '\nno. of ccn is ', np.sum(ccn_acc>1), '\nno. of BE is ', np.sum(be_acc))

print(['No. of bound entangled state is ', n_be])


print(['time is', time.perf_counter() - t1])
sio.savemat(pjoin(data_dir,'BE3x3_example.mat'), mdict={'corppy':corp_acc,'rhopy':rho_acc,'ppn':ppn_acc,'ccn':ccn_acc})





# %%
