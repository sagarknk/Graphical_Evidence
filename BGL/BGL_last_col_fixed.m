function omega_save = ...
    BGL_last_col_fixed(S,n,burnin,nmc,lambda,fixed_last_col,...
    Matrix_2be_added_Gibbs, post_mean_omega, post_mean_tau_save)

%%% S: Sample covariance matrix
%%% n: sample size
%%% burnin: burn-in for MCMC
%%% nmc: number of samples to be saved after burn-in
%%% fixed_last_col: \omega_12^*

[p] = size(S,1);
omega_save = zeros(p,p,nmc);
tau_save = zeros(p,p,nmc);

%%% S_reduced is sample covariance matrix for first p-1 columns, as the
%%% lastone remains fixed

S_reduced = S(1:(p-1), 1:(p-1));
p_reduced = p-1;
Matrix_2be_added_Gibbs_reduced = Matrix_2be_added_Gibbs(1:p_reduced, 1:p_reduced);

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
TAU_reduced(1:p_reduced, 1:p_reduced) = ...
    post_mean_tau_save(1:p_reduced, 1:p_reduced);
omega_pp = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iter = 1:(burnin + nmc)

    %if(mod(iter,1000)==0)
    %fprintf('iter = %d \n',iter);
    %end

    %%% First we update omega_pp which is nothing but sampling \omega_22
    %%% with \omega_12^* held fixed

    inv_Omega_11 = inv(Omega_reduced);
    s_22 = S(p,p);

    %%%%% sample omega_22 %%%%%%%%%%%%%
    gamma_param = gamrnd(n/2 + 1,2/(s_22+lambda));
    omega_pp =  gamma_param + ...
        fixed_last_col'*inv_Omega_11*fixed_last_col ;

    if p_reduced~=1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Omega_reduced_tilda = Omega_reduced - ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');

        temp_matrix_2be_added = (1/omega_pp)*(fixed_last_col*fixed_last_col');

        %%%%% sample Omega_11_tilde %%%%%%%%%%%%%
        for i=1:p_reduced

            ind_noi = ind_noi_all(:,i);
            s_21_tilda = S_reduced(ind_noi,i); s_22_tilda = S_reduced(i,i);
            vec_2be_added_21 = Matrix_2be_added_Gibbs_reduced(ind_noi,i);

            temp_vec_2be_added_21 = temp_matrix_2be_added(ind_noi,i);
            tau_12  = TAU_reduced(ind_noi,i);
            gamma_param_tilda = gamrnd(n/2 + 1, 2/(s_22_tilda+lambda));

            tilda_W_11 = Omega_reduced_tilda(ind_noi, ind_noi);
            inv_Omega_11 = inv(tilda_W_11);

            inv_C = diag(1./tau_12) + (s_22_tilda + lambda)*inv_Omega_11;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            mu_i = -inv_C\(s_21_tilda + vec_2be_added_21./tau_12 ...
                + temp_vec_2be_added_21./tau_12);

            inv_C_chol = chol(inv_C);
            beta = mu_i+ inv_C_chol\randn(p_reduced-1,1);

            omega_12 = beta; omega_22 = gamma_param_tilda + beta'*inv_Omega_11*beta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Omega_reduced_tilda(i,ind_noi) = omega_12; Omega_reduced_tilda(ind_noi,i) = omega_12;
            Omega_reduced_tilda(i,i) = omega_22;

            %%% TAU
            %%%% This is inverse Gaussian with mu' = sqrt(lambda^2/omega_ij^2)
            %%%% and lambda' = lambda^2

            mu_prime = sqrt(lambda^2./...
                ((omega_12 + vec_2be_added_21 + temp_vec_2be_added_21).*...
                (omega_12 + vec_2be_added_21 + temp_vec_2be_added_21)));
            lambda_prime = lambda^2;

            %%% sampler for inverse-Gaussian from wiki
            %%% https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution

%             rand_nu = randn(p_reduced-1,1);
%             rand_y = rand_nu.*rand_nu;
%             rand_x = mu_prime + (mu_prime.*mu_prime.*rand_y)./(2*lambda_prime)...
%                 -(mu_prime./(2*lambda_prime)).*sqrt(4*lambda_prime.*mu_prime.*rand_y ...
%                 + (mu_prime.*rand_y).*(mu_prime.*rand_y));
%             rand_z = rand(p_reduced-1,1);
% 
%             temp_logical = (rand_z <= (mu_prime)./(mu_prime + rand_x));
% 
%             %u_12 = rand_x.*(temp_logical) + (1-temp_logical).*(mu_prime.*mu_prime./rand_x);
%             u_12 = rand_x;
%             u_12_else = (mu_prime.*mu_prime./rand_x);
% 
%             for temp_u_12_index = 1:length(u_12)
%                 if temp_logical(temp_u_12_index,1) == 0
%                     u_12(temp_u_12_index,1) = u_12_else(temp_u_12_index,1);
%                 end
%             end
%             tau_12 = 1./u_12;

            %%% sampler for inverse-Gaussian from Generalized Inverse-Gaussian
            %%% This sampler is more efficient and stable than the sampler
            %%% commented above

            a_gig_tau = lambda_prime./(mu_prime.^2);
            b_gig_tau = lambda_prime;
            u_12 = zeros(p_reduced-1,1);
            for tau_idx = 1:p_reduced-1
                u_12(tau_idx,1) = gigrnd(-1/2,a_gig_tau(tau_idx,1), b_gig_tau,1);
            end

            tau_12 = 1./u_12;

            TAU_reduced(i,ind_noi) = tau_12;
            TAU_reduced(ind_noi,i) = tau_12;
        end

        Omega_reduced = Omega_reduced_tilda + ...
            (1/omega_pp)*(fixed_last_col*fixed_last_col');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        s_22 = S_reduced(1,1);
        gamma_param =  gamrnd(n/2 + 1, 2/(s_22+lambda));
        Omega_reduced = gamma_param + fixed_last_col'*inv(omega_pp)*fixed_last_col;
    end



    if iter > burnin
        omega_save(1:p_reduced,1:p_reduced,iter-burnin) = Omega_reduced;
        omega_save(p, 1:p_reduced, iter-burnin) = fixed_last_col';
        omega_save(1:p_reduced, p , iter-burnin) = fixed_last_col;
        omega_save(p,p,iter-burnin) = omega_pp;

        tau_save(1:p_reduced,1:p_reduced,iter-burnin) = TAU_reduced;
    end

end

end

