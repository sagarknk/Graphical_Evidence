This folder consists of files required to compute the log of marginal likelihood under the Bayesian graphical lasso (BGL) prior. 

1. BGL_graphical_evidence.m --> This matlab file implements our proposed procedure according to
the simulation settings mentioned in the paper. To generate the true precision matrix, the code calls "prior_sampling.m". With the returned output, it generates the n x p data matrix y_1,...,y_p and stores it in the folder X_mat. Similarly, the sample covariance matrix is stored in the folder S_mat. Coming to the estimation of the terms I, III and IV, the code internally calls the function 
"BGL_Hao_wang.m" to run the unrestricted MCMC sampler in the evaluation of the term IV_j. Next,
it calls "BGL_last_col_fixed.m" to run the restricted MCMC sampler in the evaluation of the term
IV_j. Terms I_j are evaluated on the fly and the sum of all terms III_j are evalauated the end.

2. Skilling_method.m --> This matlab file estimates the log of marginal likelihood by implementing
the Nested sampling approach proposed by Skilling, 2006. Internally, this calls the function 
"prior_sampling_for_Neal_and_skilling.m" to propose updates of \Omega. 

3. AIS_Neal.m --> This matlab file estimates the log of marginal likelihood by implementing
the Annealed importance sampling proposed by Neal, 2001. Internally, this calls the function 
"prior_sampling_for_Neal_and_skilling.m" to propose updates of \Omega. 


