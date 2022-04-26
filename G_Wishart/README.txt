This folder consists of files required to compute the log of marginal likelihood under the G-Wishart prior. 

1. G_Wishart_graphical_evidence.m --> This matlab file implements our proposed procedure according to the simulation settings mentioned in the paper. Before generating the data, we fix the adjacency matrix "G" and store it in the folder G_mat. We fix the scale matrix "V" to a diagonal matrix with "p" in the diagonals and store it in the folder Scale_matrix. To generate the true precision matrix, the code calls "prior_sampling.m". With the returned output, it generates the n x p data matrix y_1,...,y_p and stores it in the folder X_mat. Similarly, the sample covariance matrix is stored in the folder S_mat. The 100 different permuations of y_1,...,y_p are stored in the folder Matrix_of_random_orders. Data from all these folders is later read by AKM.R, which is explained further in #4 below. Coming to the estimation of the terms I, III and IV, the code internally calls the function "G_Wishart_Hao_wang.m" to run the unrestricted MCMC sampler in the evaluation of the term IV_j. Next, it calls "G_Wishart_last_col_fixed.m" to run the restricted MCMC sampler in the evaluation of the term IV_j. Terms I_j and III_j are evaluated on the fly, at every step of the telescoping sum. 

2. Skilling_method.m --> This matlab file estimates the log of marginal likelihood by implementing
the Nested sampling approach proposed by Skilling, 2006. Internally, this calls the function 
"prior_sampling_for_Neal_and_skilling.m" to propose updates of \Omega. 

3. AIS_Neal.m --> This matlab file estimates the log of marginal likelihood by implementing
the Annealed importance sampling proposed by Neal, 2001. Internally, this calls the function 
"prior_sampling_for_Neal_and_skilling.m" to propose updates of \Omega. 

4. AKM.R --> This R file computes the log of marginal likelihood using the method proposed by Atay-Kayis and Massam, 2005 and implemented in a R package, "BRgraph" by Mohammadi and Wit, 2015. 
