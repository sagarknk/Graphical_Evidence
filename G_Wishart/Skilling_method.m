%%% Implemented from the pseudocode (Section 6) in : 
%%% https://projecteuclid.org/journals/bayesian-analysis/volume-1/issue-4/Nested-sampling-for-general-Bayesian-computation/10.1214/06-BA127.full

rng(123456789) % set accordingly (if running in diff. workspace)
q = 5;
n = 2*q;
delta = 2;

file_name = ['./G_mat/G_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'.csv'];
G_mat_adj =  csvread(file_name);


file_name = ['./X_mat/X_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'.csv'];
xx =  csvread(file_name);


file_name = ['./Scale_matrix/Scale_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'.csv'];
scale_matrix = csvread(file_name);

log_marginal_skilling = zeros(1,25);
time_taken = zeros(1,25);

parfor skilling_iter = 1:25
    
    tic;
    N = 1e4;
    
    burnin = 2e3;
    nmc = 1e4;

    [~,prior_samples] = prior_sampling(q, burnin,nmc, G_mat_adj, ...
        scale_matrix, delta, skilling_iter);
    
    likelihood_vector = zeros(1,N);
    
    for sample_index = 1:N
        
        likelihood_vector(1,sample_index)= ...
            sum(log(mvnpdf(xx,0, inv(prior_samples(:,:,sample_index)))));
        
    end
    
    if sum(likelihood_vector == -Inf)>0
        likelihood_vector = likelihood_vector(likelihood_vector~= -Inf);
    end
    
    N = length(likelihood_vector);
    
    J = 1e4;
    Z = 0;
    X =zeros(1,J+1);
    W = zeros(1,J);
    X(1,1) = 1; %% same as X_0 = 1
    Rejection_vector = zeros(1,J);
    skilling_likelihood_vector = zeros(1,J);
    
    for i = 1:J
        
        [lowest_likelihood, rank_low_ll] = min(likelihood_vector);
        
        skilling_likelihood_vector(1,i) = lowest_likelihood;
        X(1,i+1) = exp(-i/N);
        W(1,i) = X(1,i) - X(1,i+1);
        
        %%% Z will be later incremented
        %%% need to replace the prior sample with index rank_low_ll
        
        new_greater_likelihood = lowest_likelihood - 1;
        
        temp_omega_sample = prior_samples(:,:,rank_low_ll);
        
        while new_greater_likelihood < lowest_likelihood
            
            [temp_omega_sample] = prior_sampling_for_Neal_and_Skilling(q, G_mat_adj, ...
        scale_matrix, delta, temp_omega_sample);
            
            
            new_greater_likelihood = sum(log(mvnpdf(xx,0, inv(temp_omega_sample))));
            
            Rejection_vector(1,i) = Rejection_vector(1,i)+1;
        end
        
        prior_samples(:,:,rank_low_ll) = temp_omega_sample;
        likelihood_vector(1,rank_low_ll) = new_greater_likelihood;
        
        %if mod(i,100) == 0
            %fprintf("%d iterations done\n", i);
        %end
        
    end
    
    Z = Z + sum(exp(skilling_likelihood_vector).*W);
    Z = Z + X(1,J+1)*mean(exp(likelihood_vector));
    
    log_marginal_skilling(1,skilling_iter) = log(Z);
    fprintf("%d skilling_iter is done\n", skilling_iter);

    time_taken(1, skilling_iter) = toc;


end

temp_skilling = rmoutliers(log_marginal_skilling);
mean(temp_skilling)
std(temp_skilling)

mean(time_taken)
