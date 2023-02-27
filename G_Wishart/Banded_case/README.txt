This folder consists of files required to compute the log of marginal likelihood under the G-Wishart prior (when the precision matrix is banded). 

1. G_Wishart_graphical_evidence.m
### Description: This matlab file implements our proposed procedure to compute log-marginal under the G-Wishart prior, according to the simulation settings mentioned in the paper. 
### Usage:
%%% set seed, q, n, $\delta$.
%%% set the adjacency matrix G (G_mat_adj) as a banded matrix.
%%% G = ((g_{ij})) where g_{ij}=1 if |i-j|<=k, 0 otherwise. 'k' is called the banded-parameter; for example if k=1, G is tri-diagonal and so is $\Omega$. G is stored in the folder G_mat.
%%% set the scale matrix V (Scale_matrix). It is stored in the folder Scale_matrix.
%%% True precision matrix is generated from the function call prior_sampling(q,100,100, G_mat_adj, scale_matrix, delta, 0).
%%% n x q data matrix is generated and stored in the variable "xx" and sample covariance matrix, S = xx'*xx.
%%% xx and S are stored in respective folders, which will be used to compute log-marginal from competing procedures.
%%% set burnin and nmc (number of saved mc samples). Following this, 25 random permutations of 1:q are generated and stored in Matrix_of_rand_order.
%%% log-marginal is computed in parallel for all the random permutations, using the parallel-for "parfor".
%%% Function calls to "G_Wishart_Hao_wang.m", "G_Wishart_last_col_fixed.m" are made to run the required Gibbs samplers for the procedure.
%%% Function calls to "logmvgamma.m" are made to compute the log of multivariate gamma function when required. 
%%% Mean and sd of the log-marginal from HM estimate and our procedure are printed at the end of computation. 

2. Skilling_method.m
### Description: This matlab file estimates the log of marginal likelihood by implementing the Nested sampling approach proposed by Skilling, 2006. 
### Usage:
%%% Requires the same problem dimensions and settings as specified in G_Wishart_graphical_evidence.m
%%% Reads the n x q data matrix, adjacency matrix and the scale matrix stored after data generation in G_Wishart_graphical_evidence.m
%%% Function calls are made to "prior_sampling_for_Neal_and_skilling.m" to propose $\Omega$. 
%%% log-marginal is estimated for 25 times in parallel, using the parallel-for "parfor".
%%% Mean and sd of the log-marginal estimates are printed at the end of computation. 


3. AIS_Neal.m
### Description: This matlab file estimates the log of marginal likelihood by implementing the Annealed importance sampling proposed by Neal, 2001. 
### Usage:
%%% Requires the same problem dimensions and settings as specified in G_Wishart_graphical_evidence.m
%%% Reads the n x q data matrix, adjacency matrix and the scale matrix stored after data generation in G_Wishart_graphical_evidence.m
%%% Function calls are made to "prior_sampling_for_Neal_and_skilling.m" to propose $\Omega$. 
%%% log-marginal is estimated for 25 times in parallel, using the parallel-for "parfor".
%%% Mean and sd of the log-marginal estimates are printed at the end of computation. 

4. AKM.R
### Description: This R file computes the log of marginal likelihood using the method proposed by Atay-Kayis and Massam, 2005 and implemented in a R package, "BDgraph" by Mohammadi and Wit, 2015. 
### Usage:
%%% Requires the same problem dimensions and settings as specified in G_Wishart_graphical_evidence.m
%%% Installs the R package "BDgraph" if the R package isn't installed. 
%%% Reads the scale matrix, adjacency matrix, sample covariance matrix and the 25 random permutations stored after data generation in G_Wishart_graphical_evidence.m.
%%% Function calls are made to "gnorm()", to compute the log-marginal for each of the 25 stored permutations. 
%%% Mean and sd of the log-marginal estimates are printed at the end of computation. 

