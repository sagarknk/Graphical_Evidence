function [omega_save,tau_save] = ...
    BGL_Hao_wang(S,n,burnin, nmc,lambda,Matrix_2be_added_Gibbs)

%%% S: Sample covariance matrix
%%% n: sample size
%%% burnin: burn-in for MCMC
%%% nmc: number of samples to be saved after burn-in

[p] = size(S,1);

omega_save = zeros(p,p,nmc);
tau_save = zeros(p,p,nmc);

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
Sigma = inv(Omega);

TAU(1:p,1:p) = 1;

for iter = 1:(burnin+nmc)
    
    %if(mod(iter,1000)==0)
        %fprintf('iter = %d \n',iter);
    %end
    
    %%% Gibb's sampler for Omega with Hao-Wang's decomposition
    for i = 1:p
        ind_noi = ind_noi_all(:,i);
        
        Sigma_11 = Sigma(ind_noi,ind_noi); sigma_12 = Sigma(ind_noi,i);
        sigma_22 = Sigma(i,i);
        
        s_21 = S(ind_noi,i); s_22 = S(i,i);
        vec_2be_added_21 = Matrix_2be_added_Gibbs(ind_noi,i);
        tau_12  = TAU(ind_noi,i);
        %%% sample gamma and beta
        
        gamma_param = gamrnd(n/2 + 1 , 2/(s_22 + lambda));   
        % Gamma with shape=n/2 + 1, rate=(s_22+lambda)/2
        % or scale = 2/(s_22 + lambda)
        
        inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
        
        inv_C = diag(1./tau_12) + (s_22 + lambda)*inv_Omega_11; 
        inv_C_chol = chol(inv_C);
        mu_i = -inv_C\(s_21  + vec_2be_added_21./tau_12);
        beta = mu_i+ inv_C_chol\randn(p-1,1);
        
        omega_12 = beta; omega_22 = gamma_param + beta'*inv_Omega_11*beta;
        
        %%% update Omega
        Omega(i,ind_noi) = omega_12; Omega(ind_noi,i) = omega_12;
        Omega(i,i) = omega_22;
        
        temp = inv_Omega_11*beta;
        Sigma_11 = inv_Omega_11 + temp*temp'/gamma_param;
        sigma_12 = -temp/gamma_param; sigma_22 = 1/gamma_param;
        Sigma(ind_noi,ind_noi) = Sigma_11; Sigma(i,i) = sigma_22;
        Sigma(i,ind_noi) = sigma_12; Sigma(ind_noi,i) = sigma_12;
      
        %%% TAU
        %%%% This is inverse Gaussian with mu' = sqrt(lambda^2/omega_ij^2)
        %%%% and lambda' = lambda^2
        
        mu_prime = sqrt(lambda^2./...
            ((omega_12 + vec_2be_added_21).*(omega_12 + vec_2be_added_21)));
        lambda_prime = lambda^2;
        
        %%% sampler for inverse-Gaussian from wiki
        %%% https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
        
%         rand_nu = randn(p-1,1);
%         rand_y = rand_nu.*rand_nu;
%         rand_x = mu_prime + (mu_prime.*mu_prime.*rand_y)./(2*lambda_prime)...
%             -(mu_prime./(2*lambda_prime)).*sqrt(4*lambda_prime.*mu_prime.*rand_y ...
%             + (mu_prime.*rand_y).*(mu_prime.*rand_y));
%         rand_z = rand(p-1,1);
%         
%         temp_logical = (rand_z <= (mu_prime)./(mu_prime + rand_x));
%         
%         %%% u_12 = rand_x.*(temp_logical) + (1-temp_logical).*(mu_prime.*mu_prime./rand_x);
%         %%% the above is buggy because there will be cases when rand_x is
%         %%% close to zero and temp_logical is 1. Then it results in 0 +
%         %%% 0*Inf which is NaN. Hence changing to for loop. 
%         
%         u_12 = rand_x;
%         u_12_else = (mu_prime.*mu_prime./rand_x);
%         
%         for temp_u_12_index = 1:length(u_12)
%             if temp_logical(temp_u_12_index,1) == 0
%                 u_12(temp_u_12_index,1) = u_12_else(temp_u_12_index,1);
%             end
%         end
%         
%         tau_12 = 1./u_12;
        
        %%% sampler for inverse-Gaussian from Generalized Inverse-Gaussian
        %%% This sampler is more efficient and stable than the sampler
        %%% commented above

        a_gig_tau = lambda_prime./(mu_prime.^2);
        b_gig_tau = lambda_prime;
        u_12 = zeros(p-1,1);
        for tau_idx = 1:p-1
            u_12(tau_idx,1) = gigrnd(-1/2,a_gig_tau(tau_idx,1), b_gig_tau,1);
        end
        tau_12 = 1./u_12;

        TAU(i,ind_noi) = tau_12;
        TAU(ind_noi,i) = tau_12;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iter > burnin
        omega_save(:,:,iter-burnin) = Omega;
        tau_save(:,:,iter-burnin) = TAU;
    end
    
end

end

