This folder consists of files required to compute the log of marginal likelihood under the Wishart prior. 

1. Wishart_graphical_evidence.m
### Description: This matlab file implements our proposed procedure to compute log-marginal under the Wishart prior, according to the simulation settings mentioned in the paper. 
### Usage:
%%% set seed, q, n, $\alpha$, burnin and nmc (number of saved mc samples).
%%% set the scale matrix V (Scale_matrix). It is stored in the folder Scale_matrix.
%%% True precision matrix is generated from the function call wishrnd(scale_matrix, alpha).
%%% n x q data matrix is generated and stored in the variable "xx" and sample covariance matrix, S = xx'*xx.
%%% xx and S are stored in respective folders, which will be used to compute log-marginal from competing procedures.
%%% As closed form log-marginal exists for Wishart, it is computed. Following this, 25 random permutations of 1:q are generated and 
%%% log-marginal is computed for all permutations, and parallelism is used while computing log-marginal in each permutation using the parallel-for "parfor".
%%% Function calls to "Wishart_Hao_wang.m", "Wishart_last_col_fixed.m" are made to run the required Gibbs samplers for the procedure.
%%% Mean and sd of the log-marginal from HM estimate and our procedure are printed at the end of computation. 

2. Skilling_method.m
### Description: This matlab file estimates the log of marginal likelihood by implementing the Nested sampling approach proposed by Skilling, 2006. 
### Usage:
%%% Requires the same problem dimensions and settings as specified in Wishart_graphical_evidence.m
%%% Reads the n x q data matrix and q x q scale matrix stored after data generation in Wishart_graphical_evidence.m
%%% Function calls are made to "prior_sampling_for_Neal_and_skilling.m" to propose $\Omega$. 
%%% log-marginal is estimated for 25 times in parallel, using the parallel-for "parfor".
%%% Mean and sd of the log-marginal estimates are printed at the end of computation. 


3. AIS_Neal.m
### Description: This matlab file estimates the log of marginal likelihood by implementing the Annealed importance sampling proposed by Neal, 2001. 
### Usage:
%%% Requires the same problem dimensions and settings as specified in Wishart_graphical_evidence.m
%%% Reads the q x q sample covariance matrix and scale matrix stored after data generation in Wishart_graphical_evidence.m
%%% Function calls are made to "wishrnd" to propose $\Omega$. 
%%% log-marginal is estimated for 25 times in parallel, using the parallel-for "parfor".
%%% Mean and sd of the log-marginal estimates are printed at the end of computation. 


