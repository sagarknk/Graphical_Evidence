if(!require(BDgraph)){
  install.packages("BDgraph")
  library(BDgraph)
}

set.seed(123456789)

q = 5
n = 2*q;
delta = 2
dof = 2*delta + 2 ## delta = (dof-2)/2 = 1
banded_param = 1

Scale_mat = read.csv(paste0("./Scale_matrix/Scale_mat_q_",q,"_n_",n,"_delta_",delta,"_banded_param_",banded_param,".csv"),header = FALSE)
Scale_mat = as.matrix(Scale_mat)

G_mat = read.csv(paste0("./G_mat/G_mat_q_",q,"_n_",n,"_delta_",delta,"_banded_param_",banded_param,".csv"),header = FALSE)
G_mat = as.matrix(G_mat)

S_mat = read.csv(paste0("./S_mat/S_mat_q_",q,"_n_",n,"_delta_",delta,"_banded_param_",banded_param,".csv"),header = FALSE)
S_mat = as.matrix(S_mat)

Matrix_of_random_orders = read.csv(paste0("./Matrix_of_rand_order/Random_order_mat_q_",q,"_n_",n,"_delta_",delta,"_banded_param_",banded_param,".csv"),header = FALSE)
Matrix_of_random_orders = as.matrix(Matrix_of_random_orders)

num_orders = nrow(Matrix_of_random_orders)

AKM_final_estimate = rep(0,num_orders)
AKM_estimate_prior = rep(0,num_orders)
AKM_estimate_posterior = rep(0,num_orders)

start_time = rep(0, num_orders)
end_time = rep(0, num_orders)

total_iter = 1.2e4

for(order in 1:num_orders)
{
  print(order)
  start_time[order] = Sys.time()
  G_mat_ordered = G_mat[Matrix_of_random_orders[order,],Matrix_of_random_orders[order,]]
  S_mat_ordered = S_mat[Matrix_of_random_orders[order,],Matrix_of_random_orders[order,]]
  Scale_mat_ordered = Scale_mat[Matrix_of_random_orders[order,],Matrix_of_random_orders[order,]]
  AKM_estimate_prior[order] = gnorm(adj = G_mat_ordered,b = dof,D = Scale_mat_ordered,iter = total_iter)
  
  AKM_estimate_posterior[order] = gnorm(adj = G_mat_ordered,b = dof + n,D = Scale_mat_ordered + S_mat_ordered,iter = total_iter)
  
  AKM_final_estimate[order]  = -(n*q/2)*log(2*pi) + AKM_estimate_posterior[order] - AKM_estimate_prior[order]
  
  end_time[order] = Sys.time()
}

### Removing outliers 
temp_AKM_final_est = AKM_final_estimate[!AKM_final_estimate %in% boxplot.stats(AKM_final_estimate)$out]

print(mean(temp_AKM_final_est))
print(sd(temp_AKM_final_est))
print(mean(end_time - start_time))
