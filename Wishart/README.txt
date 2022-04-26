This folder consists of files required to compute the log of marginal likelihood under the Wishart prior. 

1. Wishart_graphical_evidence.m --> This matlab file implements our proposed procedure according to
the simulation settings mentioned in the paper. At first, true precision matrix is sampled from the Wishart distribution, given the degree of freedom: alpha and the scale_matrix (stored in the folder Scale_mat) and then it generates the n x p data matrix y_1,...,y_p and stores it in the folder X_mat. Similarly, the sample covariance matrix is stored in the folder S_mat. As the exact form of log-marginal is known in the case of Wishart, it is also computed. The function "logmvgamma.m" is called to get the value of logarithm of multivariate Gamma functions, which we encounter in the normalizing constants of the Wishart. Coming to the estimation of the terms I, III and IV, the code internally calls the function "Wishart_Hao_wang.m" to run the unrestricted MCMC sampler in the evaluation of the term IV_j. Next, it calls "Wishart_last_col_fixed.m" to run the restricted MCMC sampler in the evaluation of the term IV_j. Terms I_j and III_j are evaluated on the fly, at every step of the telescoping sum. 

2. Skilling_method.m --> This matlab file estimates the log of marginal likelihood by implementing
the Nested sampling approach proposed by Skilling, 2006. 

3. AIS_Neal.m --> This matlab file estimates the log of marginal likelihood by implementing
the Annealed importance sampling proposed by Neal, 2001. 



