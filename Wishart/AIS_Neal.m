%%% Implemented from the description of Annealed Importance Samplingfrom:
%%% https://arxiv.org/pdf/1511.02543.pdf (Algorithm -1)

%%% The sample covariance matrix "S" is the input to this procedure 
rng(123456789)
q = 5;
n = 10;
alpha = 7; 

S= csvread(['./S_mat/S_mat_q_',num2str(q),'_n_',num2str(n),'_alpha_',num2str(alpha),'.csv']);
scale_matrix = csvread(['./Scale_mat/Scale_mat_q_',num2str(q),'_n_',num2str(n),'_alpha_',num2str(alpha),'.csv']);

log_marginal_AIS = zeros(1,25);
time_taken = zeros(1,25);

for neal_iter = 1:25
    tic;
    log_weight_step_1 = 0; %%% because Z_1 = 1
    delta = 0.01; %% step size for sequence of densities
    
    N = 5e3; %%% Total samples
    T = 1e2+1; %%%% Totoal intermediatary steps
    
    seq_t = linspace(0,1,T); %%% not necessary since we have step size delta
    
    log_weights_samples = zeros(1,N);
    acceptances_num = zeros(1,N);
    
    for sample_index = 1:N
        
        omega_sample = wishrnd(eye(q),alpha);
        
        log_weight_this_sample = log_weight_step_1;
        
        for density_seq = 2:T
            
            log_weight_this_sample = log_weight_this_sample + ...
                (n*delta/2)*log(det(omega_sample)) - 0.5*trace(delta*S*omega_sample)...
                -(n*q*delta/2)*log(2*pi);
            
            %%% new sample for omega
            %%% Proposal again from wishart
            
            omega_proposed = wishrnd((1/alpha)*omega_sample,alpha);
            
            %%% proposing a new omega from Wishart whose mean is at the
            %%% existing sample

            log_acceptance_prob = 0.5*(alpha + n*(density_seq-1)*delta -q-1)*log(det(omega_sample\omega_proposed))...
                -0.5*trace(((density_seq-1)*delta*S + eye(q))*(omega_proposed - omega_sample))+...
                0.5*(2*alpha-q-1)*log(det(omega_proposed\omega_sample)) - ...
                0.5*alpha*(trace(omega_proposed\omega_sample)-trace(omega_sample\omega_proposed));
            
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
    time_taken(1, neal_iter) = toc;
end

mean(log_marginal_AIS) + (n/2)*log(det(scale_matrix)) %#ok<NOPTS>
std(log_marginal_AIS)

mean(time_taken)

