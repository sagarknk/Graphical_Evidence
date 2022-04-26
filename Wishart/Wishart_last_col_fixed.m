function omega_save = ...
    Wishart_last_col_fixed(S,n,burnin,nmc,dof,fixed_last_col)

%%% S: Sample covariance matrix
%%% n: sample size
%%% burnin: burn-in for MCMC
%%% nmc: number of samples to be saved after burn-in
%%% dof: degree of freedom of wishart (alpha in the main code)
%%% fixed_last_col: \omega_12^*

[p] = size(S,1);
omega_save = zeros(p,p,nmc);

%%% S_reduced is sample covariance matrix for first p-1 columns, as the
%%% lastone remains fixed

S_reduced = S(1:(p-1), 1:(p-1));
p_reduced = p-1;

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

Omega_reduced = eye(p_reduced);
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
    
    %%%%% sample omega_22 %%%%%%%%%%%%%
    gamma = gamrnd((dof + n -p +1)/2,2/(s_22+1));
    omega_pp =  gamma + ...
        fixed_last_col'*inv_Omega_11*fixed_last_col ; 
    
    if p_reduced ~=1
        Omega_reduced_tilda = Omega_reduced - ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');
        
        %%%%% sample Omega_11_tilde %%%%%%%%%%%%%
        for i=1:p_reduced
            
            ind_noi = ind_noi_all(:,i);
            s_21_tilda = S_reduced(ind_noi,i); s_22_tilda = S_reduced(i,i);
            gamma = gamrnd((dof + n -p +1)/2, 2/(s_22_tilda+1));
            
            tilda_W_11 = Omega_reduced_tilda(ind_noi, ind_noi);
            inv_Omega_11 = inv(tilda_W_11);
      
            inv_C = (s_22_tilda+1)*inv_Omega_11; 
            mu_i = -inv_C\s_21_tilda;
            inv_C_chol = chol(inv_C);
            beta = mu_i+ inv_C_chol\randn(p_reduced-1,1);
            
            omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Omega_reduced_tilda(i,ind_noi) = omega_12; Omega_reduced_tilda(ind_noi,i) = omega_12;
            Omega_reduced_tilda(i,i) = omega_22;
        end
        
        Omega_reduced = Omega_reduced_tilda + ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');
        
    else
        s_22 = S_reduced(1,1);
        gamma =  gamrnd((dof + n -p +1)/2, 2/(s_22+1));
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

