function [omega_save] = ...
    prior_sampling_for_Neal_and_Skilling(p, G_mat_adj, scale_matrix, delta,start_point_omega)

nmc = 1;
burnin = 0;

omega_save = zeros(p,p,nmc);
%%%% ind_noi_all stores the indicices {1,2,...p}\{i} for the i^th column

ind_noi_all = zeros(p-1,p);
for i = 1:p
    if i==1
        ind_noi = [2:p]';
    elseif i==p
        ind_noi = [1:p-1]';
    else
        ind_noi = [1:i-1,i+1:p]';
    end
    
    ind_noi_all(:,i) = ind_noi;
end

% set initial values
Omega = start_point_omega;
Sigma = inv(Omega); %%% added for numerical stability

for iter = 1:(burnin+nmc)
    
    %if(mod(iter,100)==0)
        %fprintf('iter = %d \n',iter);
    %end
    
    %%% Gibb's sampler for Omega with Hao-Wang's decomposition
    
    for i = 1:p
        V_mat_22 = scale_matrix(i,i);
        ind_noi = ind_noi_all(:,i);
        V_mat_12 = scale_matrix(ind_noi,i);
        
        Sigma_11 = Sigma(ind_noi,ind_noi); sigma_12 = Sigma(ind_noi,i);
        sigma_22 = Sigma(i,i);
        
        %%% sample gamma and beta
        
        gamma_param = gamrnd(delta + 1 , 2/(V_mat_22));

        inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;

        inv_C = (V_mat_22)*inv_Omega_11; 

        %%%%%%%%%%%%%%%% finding which elements in beta to sample %%%%%%%%%
        G_mat_current_col = G_mat_adj(ind_noi,i);
        logi_which_ones = (G_mat_current_col ==1);
        
        if sum(logi_which_ones)>=1
            inv_C_required = inv_C(logi_which_ones, logi_which_ones);
            inv_C_chol_required = chol(inv_C_required);
            V_mat_12_required = V_mat_12(logi_which_ones,1);
            mu_i_reduced = -inv_C_required\V_mat_12_required;
            beta_reduced = mu_i_reduced +  inv_C_chol_required\randn(sum(logi_which_ones),1);
            
            beta = zeros(p-1,1);
            beta(logi_which_ones,1) = beta_reduced;
        else
            beta = zeros(p-1,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        omega_12 = beta; omega_22 = gamma_param + beta'*inv_Omega_11*beta;
        
        %%% update Omega and Sigma
        Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
        Omega(i,i) = omega_22;
        
        temp = inv_Omega_11*beta;
        Sigma_11 = inv_Omega_11 + temp*temp'/gamma_param;
        sigma_12 = -temp/gamma_param; sigma_22 = 1/gamma_param;
        Sigma(ind_noi,ind_noi) = Sigma_11; Sigma(i,i) = sigma_22;
        Sigma(i,ind_noi) = sigma_12; Sigma(ind_noi,i) = sigma_12;
    end
    
    if iter > burnin
        omega_save(:,:,iter-burnin) = Omega;
    end
    
end
end

