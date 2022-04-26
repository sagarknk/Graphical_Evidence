%%%% Setting seed and dimensions of the precision matrix ('p' in the paper
%%%% which is 'q' here

rng(123456789)
q = 5;
n = 2*q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% degree of freedom for wishart
alpha = q+ceil(0.3*q); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_matrix = eye(q);
for scale_col = 1:(q-1)
    scale_matrix(scale_col, scale_col+1) = 0.25;
    scale_matrix(scale_col+1, scale_col) = 0.25;
end

%sum(eig(scale_matrix)>0)
scale_matrix = (1/(alpha))*scale_matrix;
csvwrite(['./Scale_mat/Scale_mat_q_',num2str(q),'_n_',num2str(n),'_alpha_',num2str(alpha),'.csv'], scale_matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

True_Omega = wishrnd(scale_matrix, alpha);

%%%% Original data from a multivariate normal 
xx_orig = mvnrnd(zeros(1,q),inv(True_Omega),n);
%%%%% Scaled data
xx = xx_orig*sqrtm(scale_matrix);
S = xx'*xx;

csvwrite(['./X_mat/X_mat_q_',num2str(q),'_n_',num2str(n),'_alpha_',num2str(alpha),'.csv'], xx);
csvwrite(['./S_mat/S_mat_q_',num2str(q),'_n_',num2str(n),'_alpha_',num2str(alpha),'.csv'], S);

%%%% log of marginal density
%%%% Nothing but  \int \pi(\Y|\Omega)\pi(\Omega) d\Omega
%%%% with scaled data and Wishart(I, alpha) prior

log_marginal = -(n*q/2)*log(pi) + logmvgamma((alpha+n)/2, q) - ...
    ((alpha+n)/2)*log(det(eye(q) + S)) - logmvgamma(alpha/2,q);

%%%% orginal_marginal 
%%%% with unscaled data and Wishart(V, alpha) prior

log_orig_marginal =  -(n*q/2)*log(pi) + logmvgamma((alpha+n)/2, q) - ...
    ((alpha+n)/2)*log(det(inv(scale_matrix) + xx_orig'*xx_orig))...
    -(alpha/2)*log(det(scale_matrix))- logmvgamma(alpha/2,q);

%%% which is same as 

log_re_transformed_marginal = log_marginal+ (n/2)*log(det(scale_matrix));

%%% total_num_rand_orders permutations of vector 1:q
total_num_rand_orders = 25;
Matrix_of_random_orders = zeros(total_num_rand_orders, q);
for num_rand_orders = 1:total_num_rand_orders
    Matrix_of_random_orders(num_rand_orders,:) = randperm(q);
end

Matrix_of_log_ratio_likelihoods = zeros(total_num_rand_orders, q);
log_data_likelihood = zeros(total_num_rand_orders, q);
log_posterior_density = zeros(total_num_rand_orders, q);

log_prior_density = zeros(total_num_rand_orders, q);
log_prior_density_scalar_gamma = zeros(total_num_rand_orders, q);
log_prior_density_vector_normal = zeros(total_num_rand_orders, q);

Harmonic_mean_est_vec = zeros(1, total_num_rand_orders);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
burnin_MC = 1e3;
num_save_MC = 5e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for num_rand_orders = 1 :total_num_rand_orders
    
    random_order = Matrix_of_random_orders(num_rand_orders,:);
    
    fprintf("%dth random order is being computed \n", num_rand_orders);
    %random_order
    
    log_ratio_of_liklelihoods = zeros(1,q);
    %%%% temporary storage for I_{p-j+1} + III_{p-j+1} -IV_{p-j+1} for
    %%%% every j from 1 to p. 
    
    Harmonic_mean_est = zeros(1,1);
    %%% harmonic mean estimate for this given permutation. Initialised to a
    %%% zero vector of dimensions 1 x 1. It could have been a scalar zero.
    %%% But just to make a way around MATLAB's indexing rules and hence
    %%% make parfor work. 
    
    parfor num_Wishart = 1:(q-1)
        
        %fprintf("%d th num_Wishart\n", num_Wishart);
    
        reduced_data_xx = xx(:,random_order(1,1:(q-num_Wishart+1)));
        
        [~,q_reduced] = size(reduced_data_xx);
        S = reduced_data_xx'*reduced_data_xx;
        
        %%% Run the unrestricted sampler to get samples, which will
        %%% be used to approximate the Normal density in the
        %%% evaluation of the term IV
                
        omega_save = ...
            Wishart_Hao_wang(S,n,burnin_MC,num_save_MC,alpha + (1-num_Wishart));
        
        post_mean_omega = mean(omega_save,3);
        fixed_last_col = post_mean_omega(1:(q_reduced-1),q_reduced);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% computing harmonic estimate when num_BGL = 1 or when the full
        %%%%% data matrix is considered %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if num_Wishart == 1
            temp_normal_log_likli_at_posterior_samples = zeros(1,num_save_MC);
            for post_sample = 1:num_save_MC
                temp_normal_log_likli_at_posterior_samples(post_sample) = ...
                    sum(log((mvnpdf(reduced_data_xx,0, inv(omega_save(:,:,post_sample))))));
            end
            
            temp_numerical_1 = -temp_normal_log_likli_at_posterior_samples - ...
                median(-temp_normal_log_likli_at_posterior_samples);
            
            mean_temp_numerical_1 = mean(exp(temp_numerical_1));
            
            Harmonic_mean_est(num_Wishart) = -median(-temp_normal_log_likli_at_posterior_samples) ...
                -log(mean_temp_numerical_1);
  
        else
            % do nothing
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        column_number = q - num_Wishart + 1;
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
        
        vec_log_normal_density = zeros(1,num_save_MC);
        
        for sample_index  = 1:num_save_MC
            Omega_11 = omega_save(ind_noi,ind_noi, sample_index);
            inv_Omega_11 = inv(Omega_11);
            inv_C = (S(column_number, column_number)+1)*inv_Omega_11;
            s_12 = S(column_number, ind_noi);
            
            %%% mean vector of the normal (density of beta)
            
            mean_vec = -1*inv_C\s_12';
            
            vec_log_normal_density(1, sample_index)= ...
                log(mvnpdf(post_mean_omega(column_number,ind_noi),...
                mean_vec', inv(inv_C)));
        end
        
        %%% Run the restricted sampler to get samples, which will
        %%% be used to approximate the truncated Gamma density in the
        %%% evaluation of the term IV
        
        omega_save_2ndGibbs = ...
            Wishart_last_col_fixed(S,n,burnin_MC,num_save_MC,...
            alpha + (1-num_Wishart),fixed_last_col);
        
        post_mean_omega_22_2ndGibbs = mean(omega_save_2ndGibbs(q_reduced, q_reduced, :));
        
        vec_log_gamma_density = ones(1,num_save_MC);
        vec_log_gamma_density = -Inf.*vec_log_gamma_density;
        %%% We are starting with a vector of -Infinity because if the
        %%% indicator condition is not met, then the likelihood iz zero,
        %%% and the log-likelihood is -Infinity
        
        for sample_index = 1:num_save_MC
            Omega_11 = omega_save_2ndGibbs(ind_noi,ind_noi, sample_index);
            inv_Omega_11 = inv(Omega_11);
            
            temp_gamma_val = post_mean_omega_22_2ndGibbs  - ...
                fixed_last_col'*inv_Omega_11*fixed_last_col;
       
            if(temp_gamma_val > 0)
                vec_log_gamma_density(1, sample_index) = ...
                    log(gampdf(temp_gamma_val,...
                    (alpha + (1-num_Wishart) + n - q_reduced +1)/2 , ...
                    2/(S(column_number,column_number)+1)));
            else
                % do nothing
            end
            
        end
        
        log_posterior_density(num_rand_orders, num_Wishart)= ...
            log(mean(exp(vec_log_normal_density)))+...
            log(mean(exp(vec_log_gamma_density)));
        
        %%%%%%%%%% Computing <1> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Data likelihood (conditonal of multivatiate normal) %%%%%%%%%%%%
        %%%% computing mean, sd and then the data likelihood
        
        st_dev = sqrt(1/post_mean_omega_22_2ndGibbs);
        mean_vec = -1*reduced_data_xx(:,setdiff(1:q_reduced, column_number))*...
            post_mean_omega(column_number,setdiff(1:q_reduced,column_number))'...
            *st_dev*st_dev;
        
        log_data_likelihood(num_rand_orders, num_Wishart) = ...
            sum(log(normpdf(reduced_data_xx(:,column_number), mean_vec, st_dev)));
        
        
        %%%%%%%%%% Computing <3> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Prior density at posterior mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mean_vec = zeros(1, q_reduced - 1);
        covar_matrix = post_mean_omega_22_2ndGibbs*eye(q_reduced - 1);
        log_prior_density_vector_normal(num_rand_orders, num_Wishart)= ...
            log(mvnpdf(post_mean_omega(column_number,ind_noi),...
            mean_vec, covar_matrix));
        
        log_prior_density_scalar_gamma(num_rand_orders, num_Wishart) = ...
            log(gampdf(post_mean_omega_22_2ndGibbs,...
            (alpha + (1-num_Wishart))/2,2));        
        
        log_prior_density(num_rand_orders, num_Wishart) = ...
            log_prior_density_vector_normal(num_rand_orders, num_Wishart) + ...
            log_prior_density_scalar_gamma(num_rand_orders, num_Wishart);
        
        log_ratio_of_liklelihoods(1,num_Wishart) = ...
            log_data_likelihood(num_rand_orders, num_Wishart) + ...
            log_prior_density(num_rand_orders, num_Wishart) ...
            - log_posterior_density(num_rand_orders, num_Wishart);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    Harmonic_mean_est_vec(num_rand_orders) = Harmonic_mean_est;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% evaluate log f(y_1) or the last row in the telescoping sum. 
    
    xx_reduced = xx(:,random_order(1,1));
    S = xx_reduced'*xx_reduced;
    q_reduced = 1;
    omega_save_direct = zeros(q_reduced,q_reduced,num_save_MC);
    for i = 1:num_save_MC
        %omega_save_direct(:,:,i) = wishrnd(inv(eye(q_reduced) + S), alpha+(1-q)+n);
        %%% The above is same as below
        omega_save_direct(:,:,i) = gamrnd(0.5*(alpha+(1-q)+n),2/(1 + S));
    end
    
    post_mean_omega_direct = mean(omega_save_direct,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LOG_data_density = sum(log(normpdf(xx_reduced,0, inv(sqrt(post_mean_omega_direct)))));
    
    
    LOG_prior_density = (((alpha + (1-q))-q_reduced-1)/2)*log(det(post_mean_omega_direct)) ...
       -0.5*trace(post_mean_omega_direct)...
       -((alpha + (1-q))*q_reduced/2)*log(2)-logmvgamma((alpha + (1-q))/2,q_reduced);
    
    LOG_posterior_density = ((alpha + (1-q)+n-q_reduced-1)/2)*log(det(post_mean_omega_direct)) ...
       -0.5*trace((eye(q_reduced)+S)*post_mean_omega_direct)...
       +1*((alpha + (1-q)+n)/2)*log(det(eye(q_reduced)+S)) ...
       -((alpha + (1-q)+n)*q_reduced/2)*log(2)...
       -logmvgamma((alpha + (1-q)+n)/2,q_reduced);
    
    log_data_likelihood(num_rand_orders,q) = LOG_data_density;
    log_posterior_density(num_rand_orders, q) = LOG_posterior_density;
    
    log_prior_density(num_rand_orders, q) = LOG_prior_density;
    log_ratio_of_liklelihoods(1,q) = LOG_data_density + LOG_prior_density ...
        -LOG_posterior_density;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Matrix_of_log_ratio_likelihoods(num_rand_orders,:) = log_ratio_of_liklelihoods;
    
end

log_orig_marginal %#ok<NOPTS> prints the true log_marginal 

mean(sum(Matrix_of_log_ratio_likelihoods,2)) + (n/2)*log(det(scale_matrix)) %#ok<NOPTS>
%%% Above is the the mean of log marginal computed for 25 different permutations
std(sum(Matrix_of_log_ratio_likelihoods,2))
%%% Above is the the mean of log marginal computed for 25 different permutations

mean(Harmonic_mean_est_vec) + (n/2)*log(det(scale_matrix)) %#ok<NOPTS> %%% The mean of harmonic mean estimates 
%%%% for 25 different permutations
std(Harmonic_mean_est_vec) %%% The std of harmonic mean estimates 
%%%% for 25 different permutations