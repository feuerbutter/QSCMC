# QSCMC

"""
Manual for the SCMC package
created on 7th Aug 2019
by CqtLWJ
"""
## Introduction
This is the software package of our paper, [Sequentially constrained Monte Carlo sampler for quantum states](https://arxiv.org/abs/2109.14215). It utilizes the idea of Sequentially Constraint Monte Carlo to produce quantum samples with different properties. The current codes support:
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
This section should give enough information such that the user could apply this algorithm for uses other than the already supported three. 

In order to do so, one would need to modify the source code file, *QSCMC.py*. The subfunction *Prog* performs a single MC propagation with a *step_size* and a kernel (the default is the covariance matrix of the current parameters, specified by the value of *kern*, *C* for covariance matrix while *I* for identity), while *dstrProg* uses *Prog* to carry out the Sequentially Constraint MC. These two functions are structually the same for all applications. What varies across applications is the target distribution. Hence, for user-defined application, the user only needs to supply the function to calculate the target distribution for a given point in the probability space.

First, one would create a subfunction in *QSMC.py* like the form of the other *genTarSamXXX*. It should take in *N_total*, the number of points in the sample, *betak*, the number of distribution steps, *n_step*, the number of MC steps for each distribution step, *step_size*, the step size for each MC step, and other necessary parameters used to define the target distribution, such as *n_qubit*,*pop* for the *POS* case.

Inside the function body of *genTarSamXXX*, one would need to define the following subfunctions as a function of the quantum state, *rho*, and/or the probability, *corp*:
- *tarPDF(rho,corp)* returns the logarithm of target probability density function value for a given *rho* and/or *corp*.
- *constr(rho,corp,beta)* calculates the logarithm of the approximated constraint for a given sample point for a given *beta*. For example, the hyperbolic approximation of the Heaviside function.
- *refPDF(rho,corp)* returns the logarithm of reference probability density function value for a given *rho* and/or *corp*. The default is using a uniform reference sample, which is the easiest to construct. If one decides to use a different reference sample, then one has to provide both the way to calculate the PDF and the method to generate such a sample.

If the user is not sure how to construct the any of the functions, please refer to our paper for more details. After defining these functions, the intermedate function, *intPDF*, is automatically construsted according to the geometric combination we used in the paper. The user can also explore other possiblities. 

The algorithm starts by generating the reference sample *rhoref*. As from practice, we found that it is more convenient to propagate in the probability space, we would project the quantum state samples to the probability space, *corpref*, for propogation and project it back to the state space if necessary. Then a true copy of the *rhorref* and *corpref* is created as *rhonow* and *corpnow*. Then *rhonow* and *corpnow* are passed to *dstrProg* along with other necessary parameters. The *dstrProg* would then run the specified number of times and produce the final outputs, stored in *rhonow* and *corpnow*.  

## Some additional comments
### Dimension of variables
*rho*, as a stack of N dxd matrix, should have the dimension of of Nxdxd.
*corp*, as the collection of probabilities for the N points, should have a dimension of NxNprob.

### MC Calibration
*Prog* is a funcation that takes in the current point and output the next MC point. The defaullt is a Gaussian kernel, with the choice of using the *kern* flag, *I* for (I * sigeye) or *C* (Sig * sigcov) as the covariance matrix for the Gaussian kernel, where "Sig" is the covariance matrix of the current points.

One should calibrate the required number of distribution steps, *betak*, and the corresponding MC stepsize, *step_size*. *betak* should be chosen such that convergence has been reached, which is left for the user to check. *step_size* should be adjusted such that the acceptanace rate should be about 23.4% after it has stablized. 

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

