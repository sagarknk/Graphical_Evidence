%%% Implemented from the pseudocode (Section 6) in : 
%%% https://projecteuclid.org/journals/bayesian-analysis/volume-1/issue-4/Nested-sampling-for-general-Bayesian-computation/10.1214/06-BA127.full

%%% The n x p data matrix "xx" is the input to this procedure
rng(123456789)
q = 5;
n = 10;
lambda = 1;

%%% In the codes, we use the hirearchy \omega_{ij} ~ N(0,\tau_{ij}^2/\lambda^2)
%%% followed by the Inverse-gamma scale mixture of half-Cauchy densities
%%% (proposed by Makalic and Schmidt, 2016) for the half-Cauchy prior on \tau_{ij}.
%%% Hence the value of \lambda mentioned in the paper, is equivalent to 1/lambda in the code
%%% Example, if in Table 3, lambda = 6, then set lambda = 1/6 above

xx= csvread(['./X_mat/X_mat_q_',num2str(q),'_n_',num2str(n),'_lambda_',num2str(lambda),'.csv']);
log_marginal_skilling = zeros(1,25);
time_taken = zeros(1, 25);

parfor skilling_iter = 1:25
    tic;
    try
        N = 5e3;
        burnin = 1e3;
        nmc = 5e3;
        [~,prior_samples, tau_sq_samples] = prior_sampling(q, burnin,nmc, lambda, skilling_iter);
        
        likelihood_vector = zeros(1,N);
        
        for sample_index = 1:N
            
            likelihood_vector(1,sample_index)= ...
                sum(log(mvnpdf(xx,0, inv(prior_samples(:,:,sample_index)))));
            
        end
        
        if sum(likelihood_vector == -Inf)>0
            likelihood_vector = likelihood_vector(likelihood_vector~= -Inf);
        end
        
        N = length(likelihood_vector);
        
        J = 5e3;
        Z = 0;
        X =zeros(1,J+1);
        W = zeros(1,J);
        X(1,1) = 1; %% same as X_0 = 1
        Rejection_vector = zeros(1,J);
        skilling_likelihood_vector = zeros(1,J);
        
        for i = 1:J
            
            %fprintf("%d\n",i);
            [lowest_likelihood, rank_low_ll] = min(likelihood_vector);
            
            skilling_likelihood_vector(1,i) = lowest_likelihood;
            X(1,i+1) = exp(-i/N);
            W(1,i) = X(1,i) - X(1,i+1);
            
            %%% Z will be later incremented
            %%% need to replace the prior sample with index rank_low_ll
            
            new_greater_likelihood = lowest_likelihood - 1;
            
            temp_omega_sample = prior_samples(:,:,rank_low_ll);
            tau_sq_for_omega_sample = tau_sq_samples(:,:,rank_low_ll);
            
            while new_greater_likelihood < lowest_likelihood
                
                [temp_omega_sample,tau_sq_for_omega_sample] = prior_sampling_for_Neal_and_Skilling(q, lambda, ...
                    temp_omega_sample, tau_sq_for_omega_sample);
                
                
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
    catch
        fprintf("%d skilling_iter is skipped\n", skilling_iter);
    end
    
    time_taken(1, skilling_iter) = toc;
    
end

mean(log_marginal_skilling(log_marginal_skilling~=0))
std(log_marginal_skilling(log_marginal_skilling~=0))

mean(time_taken)
