# QSCMC

"""
Manual for the SCMC package
created on 7th Aug 2019
by CqtLWJ
"""
## Introduction
This package utilizes the idea of Sequentially Constraint Monte Carlo to produce quantum samples with different properties. The current codes support:
- Samples with a target posterior
- Uniform samples within a bounded likelihood region
- Samples with bound entanglement for different bipartite systems

The code has a modular structure, so it should be straight-forward to extend its application to other samples. Users are very welcome to extend application to produce different samples.

## Getting Started
This package has two main folders. One is *src*, where the source codes are located, and the other is *data*, the default directory for data storage.

The sampling can be run in the *src* directory with 'python3 main.py'. The default settings would make a 3x3 BE sample in the data folder.

## Basic Use
Users can modidy the *main.py* accordingly to produce different samples.

For all cases, users would need to supply a set of parameters to start the codes running. Examplary parameters are given, but would require modification for different target distributions. 

First, the user can choose the target sample type by setting the variable *sam_type* to be *BE*,*POS* or *BLR*. Currently, only these three are supported. Users are welcome to add in their own. 

After the *sam_type* is set, the user would need to go to the corresponding section in the *sam_type specific paras* to edit the parameters accordingly. The parameters should have been properly commented. If anything is unclear, please feel free to raise an issue or contact us directly. 

In general, the user has to give a properly defined target distribution. For instance, the power of probabilities (*pop* in the code) for a posterior, or the dimensions of the two systems for bound entangled samples (*d1* and *d2*). After that, one would have to calibrate the number of distribution steps, *betak*, MC step size, *step_size* and the corresponding scale parameters of enforcing the constraints, *scale_xxx*. Once the user is satisfied with the yield and speed, the user can set the option *save_output* to be *True* and supply a *filename* to save the final sample. The data would be stored in the *data* folder. 


## Advanced Use

## Utility Functions
This section gives a brief discription of the various utlity functions. Not all are used in the main program and some remain as legacy.

### BouEnt 
It is mainly concerned with Bound Entanglement and includes the following sub-functions:
- *Realign*     produce the realign form of the matrix
- *ParTran*     partial tranpose w.r.t the 2nd system
- *GetPnCn*     returns the minimum eigenvalue of the partial transpose and the sum of the singluar values of the realigned matrix

### GetPom
It provides the construction of several commonly used POM:
- *Trine*
- *DisTri*: Distorted Trine
- *Pauli*
- *Tetra*: Tetrahedron for 1 qubit
- *nTetra*: Tetrahedron for n qubits
- *RandPom*: random pom, "scale" should be chosen small such that I - sum(POM) > 0
- *QtsPom*: Qutrit SIC pom
- *SqrtPom*: Square root pom

### GetSam 
It is mainly concerned with producing reference samples:
- *GetWisSam*    obtain a Wishart reference sample
- *UniR(Nt,d)*   obtain a uniform sample of "Nt" real "dxd" density matrices
- *GenUniSam(Nt,d)*    obtain a uniform sample of "Nt" complex "dxd" density matrices
- *GetDiriSam*   obtain a Dirichelet sample with the powers as (pop-1)

### ImpMat 
This is used to import the variables in Matlab format
- *v73*	for storing large variables (>5Gb)
- *vo* 	for normal use, "ifskm" should be set to 1 if the to-be-imported variable is a stack of matrix, otherwise, 0 e.g. imp.vo('mle4qb.mat','rhomle',0)

### ptor
It is mainly concerned with mapping from (p)robabilities to (r)ho. 
- *getpomc* is to vectorise the pom for calculating probabilities
- *GetRho*  is to Get the corresponding (Rho) from the probabilities
- *GetP*    is to Get the (P)robabilities from rho
- *AtoR*    is to Get (R)ho from A with the Adagger A construction
- *CBtoR*   is similar to AtoR, with B/C as the real/imaginary part of A 
- *IsPhy*   is to determine if the corresponding vartype, 'r'(ho) or 'p'(robabilities), is "physical" (the mininum eigenvalue of the density matrix is > "meig"), this function returns an array of Boolean
- *IsPhy_C*  is to determine if the "rho" is postive semidefinite by using Cholesky decomposition, this is faster than IsPhy, return an integer array of '1' or '0'
- *GetLik* gets the likelihood of a particular probability point for a target posterior

### ws
This is to calculate probability for a Wishart sample:
- *UnnorCvmDen*  returns the unormalized value
- *GetLogNorm*   returns the log of the normalization 

SMCmainVerX 
	is the main run file, where one needs to provide the reference
sample, "rhoref" in the state space or "corpref" in the probability space, 
the reference distribution, "refPDF", the target distribuition, "tarPDF", 
the constraint function, if any, "Constr", and the MC propapagtion, "prog".

The functions, "refPDF", "tarPDF", "Constr" and "prog" should be a function 
of the density matrix, "rho", or the probabilities, "corp" or both. The user
can adjust accordingly.

"rho", as a stack of N dxd matrix, should have the dimension of of Nxdxd.
"corp", as the collection of probabilities for the N points, should ave 
a dimension of NxNprob.

"prog" is a funcation that takes in the current point and output the next MC
point. The defaullt is a Gaussian kernel, with the choice of using the 
(I * sigeye) or (Sig * sigcov) as the covariance matrix for the Gaussian 
kernel, where "Sig" is the covariance matrix of the current points.

One should calibrate the required number of distribution steps, "betak", and 
the corresponding MC stepsize, either "sigeye" or "sigcov".

"betak" should be chosen such that convergence has been reached, which is left 
for the user to check. "sigeye" or "sigcov" should be adjusted such that the 
acceptanace rate should be about 23.4% after it has stablized. 









