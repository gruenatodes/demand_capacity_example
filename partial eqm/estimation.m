%setup data
load_data

capacU = hpfilter(capacityutil, 10000);
capacU = capacityutil - capacU + mean(capacityutil);

wwK = hpfilter(log(workweekcapital), 10000);
wwK = log(workweekcapital) - wwK;

output = hpfilter(log(industrialprod), 10000);
output = log(industrialprod) - output;

data = cat(3, output, capacU, wwK);

% specify priors on parameters to estimate - the other parameters need to
% be set in the file 'set_params.m'
prior.rho_A = @(x) unifpdf(x,0,1);
prior.sigma_A = @(x) lognpdf(x, -2, 1);
prior.sigma_b = @(x) lognpdf(x, -2, 1);
prior.chi = @(x) lognpdf(x, 0, 1);
prior.sigma_output = @(x) lognpdf(x, -4, 1/10); % measurement error
prior.sigma_capacU = @(x) lognpdf(x, -1, sqrt(2)); % measurement error, in pct points
prior.sigma_wwK = @(x) lognpdf(x, -1/2, 1); % measurement error

prior.all = @(p) log(prior.rho_A(p.rho_A)) + ...
    log(prior.sigma_A(p.sigma_A)) + ...
    log(prior.sigma_b(p.sigma_b)) + ...
    log(prior.chi(p.chi)) + ...
    log(prior.sigma_output(p.sigma_output)) + ...
    log(prior.sigma_capacU(p.sigma_capacU)) + ...
    log(prior.sigma_wwK(p.sigma_wwK));

% set initial values
param.rho_A = 0.5;
param.sigma_A = 0.2;
param.sigma_b = 0.5;
param.chi = 1;
param.sigma_output = 0.02;
param.sigma_capacU = 10;
param.sigma_wwK = 0.1;

param_names = fieldnames(param)';
N_param = length(param_names);

% length of sample length and sd of proposal steps
N = 1000; % number of particles
M = 100000; % length of parameter sample
llhs = zeros(M,1);
parameters = zeros(M,N_param);
for ii = 1:N_param
    parameters(1,ii) = param.(param_names{ii});
end

propsd.rho_A = 0.07;
propsd.sigma_A = 0.01;
propsd.sigma_b = 0.1;
propsd.chi = 0.2;
propsd.sigma_output = 0.01;
propsd.sigma_capacU = 1.5;
propsd.sigma_wwK = 0.01;

prop_sd = ones(1, N_param);
for ii = 1:N_param
    prop_sd(ii) = propsd.(param_names{ii});
end

% initial evaluation
log_prior = prior.all(param);
llh = model_llh(param, data, N);
llhs(1) = log_prior + llh;

% sample:
rng(2)
proposal_chance = log(rand(M,1));
prop_step = randn(M, N_param);
for m = 2:M
    % proposal draw:
    prop_param_vec = parameters(m-1,:) + prop_step(m,:) .* prop_sd;
    for ii = 1:N_param
        prop_param_struct.(param_names{ii}) = prop_param_vec(ii);
    end
    
    % evaluate prior and model with proposal parameters:
    prop_prior = prior.all(prop_param_struct);
    if prop_prior > -Inf % theoretically admissible proposal
        prop_llh = model_llh(prop_param_struct, data, N);
        llhs(m) = prop_prior + prop_llh;
        if llhs(m) - llhs(m-1) > proposal_chance(m)
            accept = 1;
        else
            accept = 0;
        end
    else % reject proposal since disallowed by prior
        accept = 0;
    end
    
    % update parameters (or not)
    if accept
        parameters(m,:) = prop_param_vec;
    else
        parameters(m,:) = parameters(m-1,:);
        llhs(m) = llhs(m-1);
    end
    
    if mod(m, 100) == 0
        save current_params parameters llhs
        if mod(m, 1000) == 0
            display(['Sample ', num2str(m)])
        end
    end
    
end

