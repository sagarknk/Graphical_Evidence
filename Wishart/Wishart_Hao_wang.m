function omega_save = ...
    Wishart_Hao_wang(S,n,burnin, nmc,dof)

%%% S: Sample covariance matrix
%%% n: sample size
%%% burnin: burn-in for MCMC
%%% nmc: number of samples to be saved after burn-in
%%% dof: degree of freedom of wishart (alpha in the main code)

[p] = size(S,1);

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
Omega = eye(p); 

for iter = 1:(burnin+nmc)
    
%     if(mod(iter,500)==0)
%         fprintf('iter = %d \n',iter);
%     end
    
    %%% Gibb's sampler for Omega with Hao-Wang's decomposition
    
    for i = 1:p
        ind_noi = ind_noi_all(:,i);
        s_21 = S(ind_noi,i); s_22 = S(i,i);
        %%% sample gamma and beta
        
        gamma = gamrnd((dof + n -p +1)/2 , 2/(s_22 + 1));   
        % gamma with shape=(dof + n -p +1)/2, rate=(s_22+1)/2 or scale = 2/(s_22 + 1)
     
        inv_Omega_11 = inv(Omega(ind_noi, ind_noi));
        inv_C = (s_22 + 1)*inv_Omega_11; 
        inv_C_chol = chol(inv_C);
        mu_i = -inv_C\s_21;
        beta = mu_i+ inv_C_chol\randn(p-1,1);
        
        omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
        
        %%% update Omega, Sigma, Lambda_sq, Nu
        Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
        Omega(i,i) = omega_22;
    end
    
    if iter > burnin
        omega_save(:,:,iter-burnin) = Omega;
    end
    
end

end

