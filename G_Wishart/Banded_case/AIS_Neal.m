%%% Implemented from the description of Annealed Importance Samplingfrom:
%%% https://arxiv.org/pdf/1511.02543.pdf (Algorithm -1)

rng(123456789) 
q = 5;
n = 2*q;
delta = 2;
banded_param =1;

file_name = ['./G_mat/G_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'];
G_mat_adj =  csvread(file_name);

file_name = ['./X_mat/X_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'];
xx =  csvread(file_name);

S = xx' * xx;

file_name = ['./Scale_matrix/Scale_mat_q_',num2str(q),'_n_',num2str(n),'_delta_',num2str(delta),'_banded_param_',num2str(banded_param),'.csv'];
scale_matrix = csvread(file_name);

log_marginal_AIS = zeros(1,25);
time_taken = zeros(1,25);

parfor neal_iter = 1:25
    
    tic;
    log_weight_step_1 = 0; %%% because Z_1 = 1
    delta = 0.01; %% step size for sequence of densities
    
    N= 1e4;

    T = 1e2+1; %%%% Total intermediatary steps
        
    log_weights_samples = zeros(1,N);
    acceptances_num = zeros(1,N);
    
    %%%%%%%%%%%%%%%%% generating initial samples %%%%%%%%%%%%%%%%%%%%
    
    [init_Omega, ~] = prior_sampling(q,100,100,G_mat_adj, scale_matrix,delta,neal_iter);
    
    try
        for sample_index = 1:N
            
            [omega_sample] = prior_sampling_for_Neal_and_Skilling(q,G_mat_adj, scale_matrix,delta,...
                init_Omega);
            
            log_weight_this_sample = log_weight_step_1;
            
            for density_seq = 2:T
                
                log_weight_this_sample = log_weight_this_sample + ...
                    (n*delta/2)*log(det(omega_sample)) - 0.5*trace(delta*S*omega_sample)...
                    -(n*q*delta/2)*log(2*pi);
                
                %%% new sample for omega
                %%% Proposal again from G-Wishart
                
                [omega_proposed] = prior_sampling_for_Neal_and_Skilling(q,G_mat_adj, scale_matrix,delta,...
                omega_sample);
                
                log_acceptance_prob = 0.5*(density_seq - 1)*delta*(n*log(det(omega_sample\omega_proposed)) -...
                    trace(S*(omega_proposed - omega_sample)));
                
                if log(rand(1)) < log_acceptance_prob
                    omega_sample = omega_proposed;
                    acceptances_num(1,sample_index) = acceptances_num(sample_index)+1;
                else
                    %%%% do nothing
                end
                
            end
            
            log_weights_samples(1,sample_index) = log_weight_this_sample;
            
            %if mod(sample_index,100) ==0
                %fprintf("%d weights collected\n", sample_index);
            %end
            
        end
        log_marginal_AIS(1,neal_iter) = log(mean(exp(log_weights_samples)));
        fprintf("%d neal_iter is done\n", neal_iter);
        fprintf("%f neal_iter val\n", log_marginal_AIS(1,neal_iter));
    catch
        fprintf("Skipping %d sampler\n", neal_iter);
    end
    
    time_taken(1, neal_iter) = toc;
end

temp_AIS_estimate = rmoutliers(log_marginal_AIS);
mean(temp_AIS_estimate)
std(temp_AIS_estimate)

mean(time_taken)
