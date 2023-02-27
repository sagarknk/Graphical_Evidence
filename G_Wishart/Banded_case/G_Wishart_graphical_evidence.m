%%%% Setting seed and dimensions of the precision matrix ('p' in the paper
%%%% which is 'q' here

rng(123456789)
q = 5;
n = 2*q;
delta = 2; %%% in paper delta is alpha 

%%%% Creating a banded G with \omega_{ij}~=0 if |i-j|<= banded_param

banded_param = 1;

G_mat_adj = eye(q);

for g_mat_row = 1:(q-banded_param)
    for g_mat_col = (g_mat_row + 1):(g_mat_row +banded_param)
        G_mat_adj(g_mat_row , g_mat_col) = 1;
        G_mat_adj(g_mat_col , g_mat_row) = 1;
    end
end

if banded_param >=2
    for g_mat_row = (q-banded_param+1):(q-1)
        for g_mat_col = (g_mat_row + 1):q
            G_mat_adj(g_mat_row , g_mat_col) = 1;
            G_mat_adj(g_mat_col , g_mat_row) = 1;
        end
    end
end


csvwrite(['./G_mat/G_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], G_mat_adj);

%%%%%%%%%%%%% Scale matrix V for G-Wishart

scale_matrix = (q)*eye(q);
%scale_matrix = (2*delta + max(temp_sum_G_adj))*eye(q);

csvwrite(['./Scale_matrix/Scale_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], scale_matrix);

%sum(eig(scale_matrix)>0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[True_Omega,~] = prior_sampling(q, 100, 100, G_mat_adj, scale_matrix, delta,0);
csvwrite(['./True_precision_mat/True_Omega_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], True_Omega);

xx_orig = mvnrnd(zeros(1,q),inv(True_Omega),n);
csvwrite(['./X_mat/X_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], xx_orig);

xx = xx_orig;
S = xx'*xx;
csvwrite(['./S_mat/S_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Computing the true marginal as G is banded and decomposition is
%%%% possible

num_cliques = q - banded_param;
num_seperators = num_cliques -1;

cliques_list = cell(1,num_cliques);
seperators_list = cell(1, num_seperators);

for cliques_index =1:num_cliques
    cliques_list{1,cliques_index} = cliques_index:(cliques_index+banded_param);
end

for seperators_index = 2:num_cliques
    seperators_list{1,seperators_index-1} = (seperators_index):(seperators_index+banded_param -1);
end

%%%%%%%%%%%%%%%%%%% Prior nomalizing constant %%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_modified_clique = (2*delta + banded_param +2)*0.5;
delta_modified_seperator = (2*delta + banded_param +1)*0.5;

log_num_prior_norm_const = 0;

for cliques_index = 1:num_cliques
    
    scale_mat_clique = scale_matrix(cliques_list{1,cliques_index}, cliques_list{1,cliques_index});
    
    log_num_prior_norm_const = log_num_prior_norm_const + delta_modified_clique*(banded_param+1)*log(2) + ...
        logmvgamma(delta_modified_clique, banded_param+1) - delta_modified_clique*log(det(scale_mat_clique));
    
end

log_denom_prior_norm_const = 0;

for seperators_index = 2:num_cliques
    
    scale_mat_seperator = scale_matrix(seperators_list{1,seperators_index-1}, seperators_list{1,seperators_index-1});
    
    log_denom_prior_norm_const = log_denom_prior_norm_const + delta_modified_seperator*(banded_param)*log(2) + ...
        logmvgamma(delta_modified_seperator, banded_param) - delta_modified_seperator*log(det(scale_mat_seperator));
    
end

log_prior_norm_constant = log_num_prior_norm_const - log_denom_prior_norm_const;

%%%%%%%%%%%%%%%%%%% Posterior nomalizing constant %%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_modified_clique = (2*delta + banded_param +2 + n)*0.5;
delta_modified_seperator = (2*delta + banded_param +1 + n)*0.5;

log_num_posterior_norm_const = 0;

for cliques_index = 1:num_cliques
    
    scale_mat_clique = S(cliques_list{1,cliques_index}, cliques_list{1,cliques_index}) + ...
        scale_matrix(cliques_list{1,cliques_index}, cliques_list{1,cliques_index});
    
    log_num_posterior_norm_const = log_num_posterior_norm_const + delta_modified_clique*(banded_param+1)*log(2) + ...
        logmvgamma(delta_modified_clique, banded_param+1) - delta_modified_clique*log(det(scale_mat_clique));
    
end

log_denom_posterior_norm_const = 0;

for seperators_index = 2:num_cliques
    
    scale_mat_seperator = S(seperators_list{1,seperators_index-1}, seperators_list{1,seperators_index-1})+...
        scale_matrix(seperators_list{1,seperators_index-1}, seperators_list{1,seperators_index-1});
    
    log_denom_posterior_norm_const = log_denom_posterior_norm_const + delta_modified_seperator*(banded_param)*log(2) + ...
        logmvgamma(delta_modified_seperator, banded_param) - delta_modified_seperator*log(det(scale_mat_seperator));
    
end

log_posterior_norm_constant = log_num_posterior_norm_const - log_denom_posterior_norm_const;

%%%%%%%%%%%%%%%%%%% True log marginal %%%%%%%%%%%%%%%%%%%%%%%%%%%

log_marginal_true = -(n*q*0.5)*log(2*pi)+log_posterior_norm_constant - log_prior_norm_constant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burnin = 2e3;
nmc = 1e4;

%%% total_num_rand_orders permutations of vector 1:q
total_num_rand_orders = 25;
Matrix_of_random_orders = zeros(total_num_rand_orders, q);
for num_rand_orders = 1:total_num_rand_orders
    Matrix_of_random_orders(num_rand_orders,:) = randperm(q);
end

csvwrite(['./Matrix_of_rand_order/Random_order_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'], Matrix_of_random_orders);

Matrix_of_log_ratio_likelihoods = zeros(total_num_rand_orders, q);
log_data_likelihood = zeros(total_num_rand_orders, q);
log_posterior_density = zeros(total_num_rand_orders, q);
log_prior_density = zeros(total_num_rand_orders, q);
log_prior_density_scalar_gamma = zeros(total_num_rand_orders, q);
log_prior_density_vector_normal = zeros(total_num_rand_orders, q);

Harmonic_mean_est_vec = zeros(1, total_num_rand_orders);
our_time = zeros(1, total_num_rand_orders);
HM_time = zeros(1, total_num_rand_orders);

parfor num_rand_orders = 1:total_num_rand_orders
   
    tStart = tic;
    random_order = Matrix_of_random_orders(num_rand_orders,:);
    
    fprintf("%dth random order is being computed \n", num_rand_orders);    
    log_ratio_of_liklelihoods = zeros(1,q);
    Harmonic_mean_est = zeros(1,1);
    
    Matrix_to_be_added = zeros(q,q);
    %%%% The above matrix keeps track of the linear shifts that we are
    %%%% supposed to do every step.
    
    G_mat_adj_this_order = G_mat_adj(random_order, random_order);
    scale_matrix_this_order = scale_matrix(random_order, random_order);
    
    %%%%%%%%%%%% Starting point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [start_point_first_gibbs, ~] = prior_sampling(q, 5, 5, G_mat_adj_this_order, scale_matrix_this_order, delta,0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for num_G_Wishart = 1:q
        
        if num_G_Wishart <= q-1
            %fprintf("%d th num_G_Wishart\n", num_G_Wishart);
            reduced_data_xx = xx(:,random_order(1,1:(q-num_G_Wishart+1)));
            
            [n,q_reduced] = size(reduced_data_xx);
            S = reduced_data_xx'*reduced_data_xx;
            
            Matrix_2be_added_Gibbs = ...
                Matrix_to_be_added(1:q_reduced, 1:q_reduced);
            
            G_mat_adj_this_order_reduced = G_mat_adj_this_order(1:q_reduced, 1:q_reduced);
            scale_matrix_this_order_reduced = scale_matrix_this_order(1:q_reduced, 1:q_reduced);
            
            if num_G_Wishart == 1

                tic;
                %%% Run the unrestricted sampler to get samples, which will
                %%% be used to approximate the Normal density in the
                %%% evaluation of the term IV_j

                omega_save = ...
                    G_wishart_Hao_wang(S,n,burnin,nmc,delta, ...
                    scale_matrix_this_order_reduced,  G_mat_adj_this_order_reduced, Matrix_2be_added_Gibbs,...
                    start_point_first_gibbs);

                post_mean_omega = mean(omega_save,3);
                fixed_last_col = post_mean_omega(1:(q_reduced-1),q_reduced);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% computing harmonic estimate when num_BGL = 1 or when the full
                %%%%% data matrix is considered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                temp_normal_log_likli_at_posterior_samples = zeros(1,nmc);
                for post_sample = 1:nmc
                    temp_normal_log_likli_at_posterior_samples(post_sample) = ...
                        sum(log((mvnpdf(reduced_data_xx,0, inv(omega_save(:,:,post_sample))))));
                end

                temp_numerical_1 = -temp_normal_log_likli_at_posterior_samples - ...
                    median(-temp_normal_log_likli_at_posterior_samples);

                mean_temp_numerical_1 = mean(exp(temp_numerical_1));

                Harmonic_mean_est(num_G_Wishart) = -median(-temp_normal_log_likli_at_posterior_samples) ...
                    -log(mean_temp_numerical_1);

                HM_time(1,num_rand_orders) = toc;
            else
                %%% Run the unrestricted sampler to get samples, which will
                %%% be used to approximate the Normal density in the
                %%% evaluation of the term IV_j

                omega_save = ...
                    G_wishart_Hao_wang(S,n,burnin,nmc,delta, ...
                    scale_matrix_this_order_reduced,  G_mat_adj_this_order_reduced, Matrix_2be_added_Gibbs,...
                    start_point_first_gibbs);

                post_mean_omega = mean(omega_save,3);
                fixed_last_col = post_mean_omega(1:(q_reduced-1),q_reduced);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            column_number = q - num_G_Wishart + 1;
            %%% Also to be noted that column_number = q_reduced
 
            %%%%%%%%%%% computing <4> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% Posterior likelihood at posterior mean %%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% ind_noi_all has indices of non-diagonal entries for each column
            %%% for first column it has numbers 2 to q-1. For second column it has
            %%% numbers 1, {3,...,q} and so on.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ind_noi_all = zeros(q_reduced-1,q_reduced);
            for i = 1:q_reduced
                if i==1
                    ind_noi = [2:q_reduced]';
                elseif i==q_reduced
                    ind_noi = [1:q_reduced-1]';
                else
                    ind_noi = [1:i-1,i+1:q_reduced]';
                end
                
                ind_noi_all(:,i) = ind_noi;
            end
            
            ind_noi = ind_noi_all(:,column_number);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            vec_log_normal_density = zeros(1,nmc);
            V_mat_22 = scale_matrix_this_order_reduced(column_number,column_number);
            V_mat_12 = scale_matrix_this_order_reduced(ind_noi,column_number);
            s_21 = S(ind_noi,column_number);
            s_22 = S(column_number,column_number);
            vec_2be_added_21 =  -1 * Matrix_2be_added_Gibbs(ind_noi,column_number);
            %%% Note the -1 above. This is done to make sure to respect the
            %%% crucial indicator function.
            %%%%%%%%%%%%%%%% finding which elements in beta to evaluate the %%%
            %%%%%%%%% likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            G_mat_current_col = G_mat_adj_this_order_reduced(ind_noi,column_number);
            logi_which_ones = (G_mat_current_col ==1);
            logi_which_zeros = (G_mat_current_col ==0);
            
            if sum(logi_which_ones)>=1
                V_mat_12_required = V_mat_12(logi_which_ones,1);
                s_21_required = s_21(logi_which_ones,1);
                
                if sum(logi_which_zeros)>=1
                    vec_2be_added_21_required = vec_2be_added_21(logi_which_zeros,1);
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for sample_index  = 1:nmc
                Omega_11 = omega_save(ind_noi,ind_noi, sample_index);
                inv_Omega_11 = inv(Omega_11);
                inv_C = (s_22 + V_mat_22)*inv_Omega_11;
                
                if sum(logi_which_ones)>=1
                    inv_C_required = inv_C(logi_which_ones, logi_which_ones);
                    
                    if sum(logi_which_zeros)>=1
                        
                        inv_C_not_required = inv_C(logi_which_zeros, logi_which_ones);
                        
                        vec_2be_added_21_required_modified = ...
                            vec_2be_added_21_required' * inv_C_not_required;
                        
                        vec_2be_added_21_required_modified = ...
                            vec_2be_added_21_required_modified' ;
                        
                        mu_i_reduced = -inv_C_required\...
                            (V_mat_12_required + s_21_required + vec_2be_added_21_required_modified);
                    else
                        mu_i_reduced = -inv_C_required\(V_mat_12_required + s_21_required);
                    end
                    
                    mean_vec = mu_i_reduced;
                    
                    vec_log_normal_density(1, sample_index)= ...
                        log(mvnpdf(post_mean_omega(column_number,logi_which_ones),...
                        mean_vec', inv(inv_C_required)));
                else
                    vec_log_normal_density(1, sample_index)= 0;
                end
            end
            
            %%% Run the restricted sampler to get samples, which will
            %%% be used to approximate the truncated Gamma density in the
            %%% evaluation of the term IV_j
            
            omega_save_2ndGibbs = ...
                 G_wishart_last_col_fixed(S,n,burnin, nmc,...
                delta,fixed_last_col,...
                scale_matrix_this_order_reduced,  G_mat_adj_this_order_reduced, Matrix_2be_added_Gibbs,...
                post_mean_omega);
            
            post_mean_omega_22_2ndGibbs = mean(omega_save_2ndGibbs(q_reduced, q_reduced, :));
            
            %%%%%%%%%%%%% Updating the start_point_first_gibbs %%%%%%%%%%%%
            start_point_first_gibbs = mean(omega_save_2ndGibbs,3);
            start_point_first_gibbs = start_point_first_gibbs(1:(q_reduced-1), 1:(q_reduced-1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            vec_log_gamma_density = ones(1,nmc);
            vec_log_gamma_density = -Inf.*vec_log_gamma_density;
            %%% We are starting with a vector of -Infinity because if the
            %%% indicator condition is not met, then the likelihood iz zero,
            %%% and the log-likelihood is -Infinity
            
            for sample_index = 1:nmc
                Omega_11 = omega_save_2ndGibbs(ind_noi,ind_noi, sample_index);
                inv_Omega_11 = inv(Omega_11);
                
                temp_gamma_val = post_mean_omega_22_2ndGibbs  - ...
                    fixed_last_col'*inv_Omega_11*fixed_last_col;
                
                if(temp_gamma_val > 0)
                    vec_log_gamma_density(1, sample_index) = ...
                        log(gampdf(temp_gamma_val,...
                        delta + n/2 + 1 , ...
                        2/(S(column_number,column_number)+...
                        scale_matrix_this_order_reduced(column_number, column_number))));
                else
                    % do nothing
                end
                
            end
            
            log_posterior_density(num_rand_orders, num_G_Wishart)= ...
                log(mean(exp(vec_log_normal_density)))+...
                log(mean(exp(vec_log_gamma_density)));
            
            %%%%%%%%%% Computing <1> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% Data likelihood (conditonal of multivatiate normal) %%%%%%%%%%%%
            %%%% computing mean, sd and then the data likelihood
            
            st_dev = sqrt(1/post_mean_omega_22_2ndGibbs);
            
            mean_vec = -1*reduced_data_xx(:,setdiff(1:q_reduced, column_number))*...
                post_mean_omega(column_number,setdiff(1:q_reduced,column_number))'...
                *st_dev*st_dev;
            
            log_data_likelihood(num_rand_orders, num_G_Wishart) = ...
                sum(log(normpdf(reduced_data_xx(:,column_number), mean_vec, st_dev)));
            
            %%%%%%%%%% Computing <3> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% Prior density at posterior mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            V_mat_22 = scale_matrix_this_order_reduced(column_number,column_number);
            V_mat_12 = scale_matrix_this_order_reduced(ind_noi,column_number);
            vec_2be_added_21 =  -1 * Matrix_2be_added_Gibbs(ind_noi,column_number);
            
            %%% Note the -1 above. This is done to make sure to respect the
            %%% crucial indicator function.
            %%%%%%%%%%%%%%%% finding which elements in beta to evaluate the %%%
            %%%%%%%%% likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            G_mat_current_col = G_mat_adj_this_order_reduced(ind_noi,column_number);
            logi_which_ones = (G_mat_current_col ==1);
            logi_which_zeros = (G_mat_current_col ==0);
            
            if sum(logi_which_ones)>=1
                V_mat_12_required = V_mat_12(logi_which_ones,1);
                
                if sum(logi_which_zeros)>=1
                    vec_2be_added_21_required = vec_2be_added_21(logi_which_zeros,1);
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            inv_C = (post_mean_omega_22_2ndGibbs^(-1))*scale_matrix_this_order_reduced(ind_noi,ind_noi);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if sum(logi_which_ones)>=1
                inv_C_required = inv_C(logi_which_ones, logi_which_ones);
                
                if sum(logi_which_zeros)>=1
                    
                    inv_C_not_required = inv_C(logi_which_zeros, logi_which_ones);
                    
                    vec_2be_added_21_required_modified = ...
                        vec_2be_added_21_required' * inv_C_not_required;
                    
                    vec_2be_added_21_required_modified = ...
                        vec_2be_added_21_required_modified' ;
                    
                    mu_i_reduced = -inv_C_required\...
                        (V_mat_12_required + vec_2be_added_21_required_modified);
                else
                    mu_i_reduced = -inv_C_required\(V_mat_12_required + s_21_required);
                end
                
                mean_vec = mu_i_reduced;
                
                log_prior_density_vector_normal(num_rand_orders, num_G_Wishart)= ...
                    log(mvnpdf(post_mean_omega(column_number,logi_which_ones),...
                    mean_vec', inv(inv_C_required)));
                    
            else
                log_prior_density_vector_normal(num_rand_orders, num_G_Wishart)=0;
            end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% Computing the GIG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (sum(logi_which_zeros)==0)
                log_prior_density_scalar_gamma(num_rand_orders, num_G_Wishart) = ...
                log(gampdf(post_mean_omega_22_2ndGibbs,...
                delta  + sum(logi_which_ones)/2 + 1 ,2/V_mat_22));
           
            else
                GIG_a = V_mat_22;
                
                diag_V_mat = diag(scale_matrix_this_order_reduced);
                diag_V_mat_required = diag_V_mat(logi_which_zeros,1);
                vec_2be_added_21 =  -1 * Matrix_2be_added_Gibbs(ind_noi,column_number);
                vec_2be_added_21_required = vec_2be_added_21(logi_which_zeros,1);
                
                GIG_b = sum(diag_V_mat_required.*vec_2be_added_21_required.* vec_2be_added_21_required);
                
                if GIG_b == 0
                    log_prior_density_scalar_gamma(num_rand_orders, num_G_Wishart) = ...
                        log(gampdf(post_mean_omega_22_2ndGibbs,...
                        delta  + sum(logi_which_ones)/2 + 1 ,2/V_mat_22));
                else
                    GIG_p =  delta  + sum(logi_which_ones)/2 + 1;
                    log_prior_density_scalar_gamma(num_rand_orders, num_G_Wishart) = ...
                        (GIG_p/2)*log(GIG_a/GIG_b) - log(2)...
                        -log(besselk(GIG_p, sqrt(GIG_a*GIG_b))) ...
                        +(GIG_p - 1)*log(post_mean_omega_22_2ndGibbs) ...
                        -0.5*(GIG_a*post_mean_omega_22_2ndGibbs + (GIG_b/post_mean_omega_22_2ndGibbs));
                    
                end
            end
   
            log_prior_density(num_rand_orders, num_G_Wishart) = ...
                log_prior_density_vector_normal(num_rand_orders, num_G_Wishart) + ...
                log_prior_density_scalar_gamma(num_rand_orders, num_G_Wishart);
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            log_ratio_of_liklelihoods(1,num_G_Wishart) = ...
                log_data_likelihood(num_rand_orders, num_G_Wishart) + ...
                log_prior_density(num_rand_orders, num_G_Wishart) ...
                - log_posterior_density(num_rand_orders, num_G_Wishart);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Matrix_to_be_added(1:(q-num_G_Wishart), 1:(q-num_G_Wishart)) = ...
                Matrix_to_be_added(1:(q-num_G_Wishart), 1:(q-num_G_Wishart)) + ...
                (1/post_mean_omega_22_2ndGibbs)*(fixed_last_col*fixed_last_col');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            Harmonic_mean_est_vec(num_rand_orders) = Harmonic_mean_est;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% To evaluate log f(y_1) or the marginal likelihood of the
            %%% remaining first column, it's like
            %%% imposing a univariate G-wishart prior on the precision of a
            %%% univariate normal random variable
            
            xx_reduced = xx(:,random_order(1,1));
            S = xx_reduced'*xx_reduced;
            q_reduced = 1;
            
            V_mat_22 = scale_matrix_this_order_reduced(q_reduced,q_reduced);
            LOG_marginal_first_col_only = -(n*q_reduced/2)*log(pi) + ...
                logmvgamma(delta+n/2+1, q_reduced) - ...
                (delta + n/2 + 1)*log(det(V_mat_22 + S)) - logmvgamma(delta +1,q_reduced) + ...
                (delta + 1)*log(V_mat_22);
            
            log_ratio_of_liklelihoods(1,q) = LOG_marginal_first_col_only;
            
            Matrix_of_log_ratio_likelihoods(num_rand_orders,:) = log_ratio_of_liklelihoods;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
    
    our_time(1, num_rand_orders) = toc(tStart);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% True log marginal %%%%%%%%%%%%%%%%%%%%%%%%%%%

log_marginal_true %#ok<NOPTS>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_sum = sum(Matrix_of_log_ratio_likelihoods,2);
temp_sum = temp_sum(~isinf(temp_sum));
[temp_sum, rows_removed] = rmoutliers(temp_sum);
length(temp_sum)

mean(temp_sum) %%% mean of non-outliers log-marginal estimates (of our proposed procedure)
std(temp_sum) %%% std of non-outliers log-marginal estimates (of our proposed procedure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean(Harmonic_mean_est_vec) %%% The mean of log marginal computed 
%%%% for 25 different permutations
std(Harmonic_mean_est_vec) %%% The std of log marginal computed 
%%%% for 25 different permutations

mean(our_time)
mean(HM_time)
