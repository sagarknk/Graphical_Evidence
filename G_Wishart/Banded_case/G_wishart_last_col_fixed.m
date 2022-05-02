function omega_save = ...
    G_wishart_last_col_fixed_new(S,n,burnin,nmc,delta,fixed_last_col,...
    scale_matrix_this_order_reduced,  G_mat_adj_this_order_reduced, Matrix_2be_added_Gibbs,...
    post_mean_omega)

%%% S: Sample covariance matrix
%%% n: sample size
%%% burnin: burn-in for MCMC
%%% nmc: number of samples to be saved after burn-in
%%% fixed_last_col: \omega_12^*

[p] = size(S,1);
omega_save = zeros(p,p,nmc);

%%% S_reduced is sample covariance matrix for first p-1 columns, as the
%%% lastone remains fixed

S_reduced = S(1:(p-1), 1:(p-1));
p_reduced = p-1;
Matrix_2be_added_Gibbs_reduced = Matrix_2be_added_Gibbs(1:p_reduced, 1:p_reduced);
scale_matrix_required = scale_matrix_this_order_reduced(1:p_reduced, 1:p_reduced);
G_mat_adj_required = G_mat_adj_this_order_reduced(1:p_reduced, 1:p_reduced);

%%%% ind_noi_all stores the indicices {1,2,...p}\{i} for the i^th column

if p_reduced ~=1
    ind_noi_all = zeros(p_reduced-1,p_reduced);
    for i = 1:p_reduced
        if i==1
            ind_noi = [2:p_reduced]';
        elseif i==p_reduced
            ind_noi = [1:p_reduced-1]';
        else
            ind_noi = [1:i-1,i+1:p_reduced]';
        end
        
        ind_noi_all(:,i) = ind_noi;
    end
else
    % do nothing 
end

Omega_reduced = post_mean_omega(1:p_reduced, 1:p_reduced);

omega_pp = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:(burnin + nmc)

%fprintf("%d\n", iter);    
%     if(mod(iter,500)==0)
%         fprintf('iter = %d \n',iter);
%     end
    
    %%% First we update omega_pp which is nothing but sampling \omega_22
    %%% with \omega_12^* held fixed
    
    inv_Omega_11 = inv(Omega_reduced);
    s_22 = S(p,p);
    V_mat_22 = scale_matrix_this_order_reduced(p,p);
    
    %%%%% sample omega_22 %%%%%%%%%%%%%
    gamma = gamrnd(delta + n/2 + 1,2/(s_22+V_mat_22));
    omega_pp =  gamma + ...
        fixed_last_col'*inv_Omega_11*fixed_last_col ; 
    
    if p_reduced ~=1
        Omega_reduced_tilda = Omega_reduced - ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');
        
        temp_matrix_2be_added = (1/omega_pp)*(fixed_last_col*fixed_last_col');
        
        %%%%% sample Omega_11_tilde %%%%%%%%%%%%%
        for i=1:p_reduced
            
            V_mat_22 = scale_matrix_required(i,i);
            ind_noi = ind_noi_all(:,i);
            V_mat_12 = scale_matrix_required(ind_noi,i);
            
            s_21_tilda = S_reduced(ind_noi,i); s_22_tilda = S_reduced(i,i);
            vec_2be_added_21 = -1 * Matrix_2be_added_Gibbs_reduced(ind_noi,i);
            temp_vec_2be_added_21 = -1 * temp_matrix_2be_added(ind_noi,i);
            
            %%% Note the -1 above. This is done to make sure to respect the
            %%% crucial indicator function.
            
            gamma = gamrnd(delta + n/2 + 1, 2/(s_22_tilda+V_mat_22));
            
            tilda_W_11 = Omega_reduced_tilda(ind_noi, ind_noi);
            inv_Omega_11 = inv(tilda_W_11);
            
            inv_C = (s_22_tilda+V_mat_22)*inv_Omega_11;
            %%%%%%%%%%%%%%%% finding which elements in beta to sample %%%%%%%%%
            G_mat_current_col = G_mat_adj_required(ind_noi,i);
            logi_which_ones = (G_mat_current_col ==1);
            logi_which_zeros = (G_mat_current_col ==0);
            
            if sum(logi_which_ones)>=1
                inv_C_required = inv_C(logi_which_ones, logi_which_ones);
                inv_C_chol_required = chol(inv_C_required);
                
                V_mat_12_required = V_mat_12(logi_which_ones,1);
                s_21_tilda_required = s_21_tilda(logi_which_ones,1);
                
                if sum(logi_which_zeros)>=1
                    vec_2be_added_21_required = vec_2be_added_21(logi_which_zeros,1);
                    inv_C_not_required = inv_C(logi_which_zeros, logi_which_ones);
                    
                    vec_2be_added_21_required_modified = ...
                        vec_2be_added_21_required' * inv_C_not_required;
                    
                    vec_2be_added_21_required_modified = ...
                        vec_2be_added_21_required_modified' ;
                    
                    temp_vec_2be_added_21_required = temp_vec_2be_added_21(logi_which_zeros,1);
                    
                    temp_vec_2be_added_21_required_modified = ...
                        temp_vec_2be_added_21_required' * inv_C_not_required;
                    
                    temp_vec_2be_added_21_required_modified = ...
                        temp_vec_2be_added_21_required_modified' ;
                    
                    
                    mu_i_reduced = -inv_C_required\...
                        (V_mat_12_required + s_21_tilda_required + ...
                        vec_2be_added_21_required_modified + ...
                        temp_vec_2be_added_21_required_modified);
                    
                    beta_reduced = ...
                        mu_i_reduced +  inv_C_chol_required\randn(sum(logi_which_ones),1);
                    
                else
                    mu_i_reduced = -inv_C_required\(V_mat_12_required + s_21_tilda_required);
                    beta_reduced = mu_i_reduced +  inv_C_chol_required\randn(sum(logi_which_ones),1);
                end
                
                
                beta = zeros(p_reduced-1,1);
                
                if sum(logi_which_ones)== p_reduced-1
                    beta(logi_which_ones,1) = beta_reduced;
                else
                    beta(logi_which_ones,1) = beta_reduced;
                    beta(logi_which_zeros,1) = vec_2be_added_21_required + ...
                        temp_vec_2be_added_21_required;
                end
                
            else
                beta =vec_2be_added_21 + temp_vec_2be_added_21;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Omega_reduced_tilda(i,ind_noi) = omega_12; Omega_reduced_tilda(ind_noi,i) = omega_12;
            Omega_reduced_tilda(i,i) = omega_22;
        end
        
        Omega_reduced = Omega_reduced_tilda + ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');
        
    else
        s_22 = S_reduced(1,1);
        V_mat_22 = scale_matrix_required(1,1);
        gamma =  gamrnd(delta + n/2 + 1, 2/(s_22+V_mat_22));
        Omega_reduced = gamma + fixed_last_col'*inv(omega_pp)*fixed_last_col;
    end

    if iter > burnin
        omega_save(1:p_reduced,1:p_reduced,iter-burnin) = Omega_reduced;
        omega_save(p, 1:p_reduced, iter-burnin) = fixed_last_col';
        omega_save(1:p_reduced, p , iter-burnin) = fixed_last_col;
        omega_save(p,p,iter-burnin) = omega_pp;
    end
    
end

end

