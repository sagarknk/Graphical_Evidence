%%% Implemented from the description of Annealed Importance Samplingfrom:
%%% https://arxiv.org/pdf/1511.02543.pdf (Algorithm -1)

%%% The sample covariance matrix "S" is the input to this procedure 
rng(123456789)
q = 5;
n = 10;
lambda = 1;

S= csvread(['./S_mat/S_mat_q_',num2str(q),'_n_',num2str(n),'_lambda_',num2str(lambda),'.csv']);
log_marginal_AIS = zeros(1,25);
time_taken = zeros(1,25);

parfor neal_iter = 1:25
    
    tic;
    log_weight_step_1 = 0; %%% because Z_1 = 1
    delta = 0.01; %% step size for sequence of densities
    
    N = 5e3; %%% Total samples
    T = 1e2+1; %%%% Total intermediatary steps
    
    log_weights_samples = zeros(1,N);
    acceptances_num = zeros(1,N);
    
    %%%%%%%%%%%%%%%%% generating initial samples %%%%%%%%%%%%%%%%%%%%
    
    %%% prior_sampling samples from the BGL prior and sets the initial value of 
    %%% the precision matrix and the latent parameters \tau = (\tau_{ij})
    
    [init_Omega, ~, tau_samples] = prior_sampling(q,100,100,lambda,neal_iter);
    init_tau_sample = tau_samples(:,:,100);
    
    try
        for sample_index = 1:N
            
            [omega_sample,tau_sample] = prior_sampling_for_Neal_and_skilling(q,lambda,init_Omega,...
                init_tau_sample);
            
            log_weight_this_sample = log_weight_step_1;
            
            for density_seq = 2:T
                
                log_weight_this_sample = log_weight_this_sample + ...
                    (n*delta/2)*log(det(omega_sample)) - 0.5*trace(delta*S*omega_sample)...
                    -(n*q*delta/2)*log(2*pi);
                
                %%% new sample for the precision matrix
                %%% Proposal again from BGL
                
                [omega_proposed,tau_proposed] = prior_sampling_for_Neal_and_skilling(q,lambda,omega_sample,...
                    tau_sample);
                
                log_acceptance_prob = 0.5*(density_seq - 1)*delta*(n*log(det(omega_sample\omega_proposed)) -...
                    trace(S*(omega_proposed - omega_sample)));
                
                if log(rand(1)) < log_acceptance_prob
                    omega_sample = omega_proposed;
                    tau_sample = tau_proposed;
                    acceptances_num(1,sample_index) = acceptances_num(sample_index)+1;
                else
                    %%%% do nothing
                end
                
            end
            
            log_weights_samples(1,sample_index) = log_weight_this_sample;
            
            %if mod(sample_index,100) ==0
            %    fprintf("%d weights collected\n", sample_index);
            %end
            
        end
        log_marginal_AIS(1,neal_iter) = log(mean(exp(log_weights_samples)));
        fprintf("%d neal_iter is done\n", neal_iter);
    catch
        fprintf("Skipping %d sampler\n", neal_iter);
    end
    time_taken(1, neal_iter) = toc;
end

mean(log_marginal_AIS(log_marginal_AIS~=0))
std(log_marginal_AIS(log_marginal_AIS~=0))

mean(time_taken)
