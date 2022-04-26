function [True_Omega, omega_save,tau_sq_save] = prior_sampling(p, BURNIN, NMC, lambda, seed_i)

rng(123456789 + 100*seed_i)
nmc = NMC;
burnin = BURNIN;

omega_save = zeros(p,p,nmc);
tau_sq_save = zeros(p,p,nmc);

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
Sigma = inv(Omega); %%% added for numerical stability
TAU_sq(1:p,1:p) = 1;
Nu(1:p,1:p) = 1;

for iter = 1:(burnin+nmc)
    
    %if(mod(iter,100)==0)
        %fprintf('iter = %d \n',iter);
    %end
    
    %%% Gibb's sampler for Omega with Hao-Wang's decomposition
    
    for i = 1:p
        ind_noi = ind_noi_all(:,i);
        Sigma_11 = Sigma(ind_noi,ind_noi); sigma_12 = Sigma(ind_noi,i);
        sigma_22 = Sigma(i,i);
        %s_21 = S(ind_noi,i); s_22 = S(i,i);
        tau_sq_12  = TAU_sq(ind_noi,i);
        nu_12 = Nu(ind_noi,i);
        %%% sample gamma and beta
        
        gamma_param = gamrnd(1 , 2*lambda);
        
        inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
        
        inv_C = diag(1./(tau_sq_12*lambda^2)) + (1/lambda)*inv_Omega_11; 
        inv_C_chol = chol(inv_C);
        
        beta = inv_C_chol\randn(p-1,1);
        
        omega_12 = beta; omega_22 = gamma_param + beta'*inv_Omega_11*beta;
        
        %%% sample tau_sq and nu
        rate = omega_12.^2/(2*lambda^2)+1./nu_12;
        tau_sq_12 = 1./gamrnd(1,1./rate);    % random inv gamma with shape=1, rate=rate
        nu_12 = 1./gamrnd(1,1./(1+1./tau_sq_12));    % random inv gamma with shape=1, rate=1+1/lambda_sq_12
        
        %%% update Omega and Sigma
        Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
        Omega(i,i) = omega_22;
        temp = inv_Omega_11*beta;
        Sigma_11 = inv_Omega_11 + temp*temp'/gamma_param;
        sigma_12 = -temp/gamma_param; sigma_22 = 1/gamma_param;
        Sigma(ind_noi,ind_noi) = Sigma_11; Sigma(i,i) = sigma_22;
        Sigma(i,ind_noi) = sigma_12; Sigma(ind_noi,i) = sigma_12;
        TAU_sq(i,ind_noi) = tau_sq_12; TAU_sq(ind_noi,i) = tau_sq_12;
        Nu(i,ind_noi) = nu_12; Nu(ind_noi,i) = nu_12;
        
    end
    
    if iter > burnin
        omega_save(:,:,iter-burnin) = Omega;
        tau_sq_save(:,:,iter-burnin) = TAU_sq;
    end
    
end

True_Omega = omega_save(:,:,nmc);
end
