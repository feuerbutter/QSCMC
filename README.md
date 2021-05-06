# QSCMC

ReadMe.txt
Long ago
8 Aug 2019
H
Han Rui edited an item
Text
ReadMe.txt
7 Aug 2019
w
You uploaded an item
Text
ReadMe.txt
"""
Manual for the SCMC package
created on 7th Aug 2019
by CqtLWJ
"""

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

GetPom
	provides tthe construction of several commonly used POM.
	-Trine
	-Distorted Trine
	-Pauli
	-Tetra
	-nTetrahedron(n)
	-random pom, "scale" should be chosen small such that I - sum(POM) > 0
	-Qutrit SIC pom
	-Square root pom
ptor
	is mainly concerned with mapping from (p)robabilities to (r)ho. 
	-getpomc is to vectorise the pom for calculating probabilities
	-GetRho  is to Get the corresponding (Rho) from the probabilities
	-GetP    is to Get the (P)robabilities from rho
	-AtoR    is to Get (R)ho from A with the Adagger A construction
	-CBtoR   is similar to AtoR, with B/C as the real/imaginary part of A 
	-IsPhy   is to determine if the corresponding vartype, 'r'(ho) or 
	 	 'p'(robabilities), is "physical" (the mininum eigenvalue
		  of the density matrix is > "meig"), this function returns 
		  an array of Boolean
	-IsPhy_C  is to determine if the "rho" is postive semidefinite by using 
		  Cholesky decomposition, this is faster than IsPhy, return an 
		  integer array of '1' or '0'

GetSam 
	is mainly concerned with producing reference samples
	-GetWisSam    obtain a Wishart reference sample
	-UniR(Nt,d)   obtain a uniform sample of "Nt" real "dxd" density matrices
	-Uni(Nt,d)    obtain a uniform sample of "Nt" complex "dxd" density matrices
	-GetDiriSam   obtain a Dirichelet sample with the powers as (pop-1)

BouEnt 
	is mainly concerned with Bound Entanglement 
	-Realign     produce the realign form of the matrix
	-ParTran     partial tranpose w.r.t the 2nd system
	-GetPnCn     returns the minimum eigenvalue of the partial transpose and the
		      sum of the singluar values of the realigned matrix
ws
	is to calculate probability for a Wishart Sample
	-UnnorCvmDen  returns the unormalized value
	-GetLogNorm   returns the log of the normalization 

ImpMat 
	is to import the variables in Matlab format
	-v73	for storing large variables (>5Gb)
	-vo 	for normal use, "ifskm" should be set to 1 if the to-be-imported variable
		is a stack of matrix, otherwise, 0 e.g. imp.vo('mle4qb.mat','rhomle',0)