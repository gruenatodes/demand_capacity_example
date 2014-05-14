function [LLH] = model_llh(params, data, N)
[p, pp] = set_params(); % get default parameters
% override defaults with specified parameters - have to be named identical
% of course:
param_names = fieldnames(params);
for ii = 1:length(param_names)
    p.(param_names{ii}) = params.(param_names{ii});
end
p.mu_b = -(p.sigma_b ^ 2 / 2);
p.mu_A = -(p.sigma_A ^ 2 / 2); 
pp.sigexp = (p.sigma - 1) / p.sigma; % short for exponent in DS function

T = length(data);
z = 1;

agg.Y = 1;
agg.P = 1;
agg.w = 2/3 * agg.Y;
agg.R = 1 / p.beta;

[A, pi.A] = TAUCHEN(pp.NA, p.mu_A, p.rho_A, p.sigma_A, pp.tauchen_param_A);
demand_shocks = exp(A);
lr_probs = (ones(1, pp.NA) / pp.NA) * pi.A ^ 25;
% solve model for each of the aggregate states.
pol_output = zeros(1, pp.NA); pol_capacU = zeros(1, pp.NA); pol_wwK = zeros(1, pp.NA);
for iA = 1:pp.NA
    result = simulation(p, pp, agg, z, demand_shocks(iA));
    pol_output(iA) = log(result.output);
    pol_capacU(iA) = result.prod.avg_prod / result.max_ys * 100;
    pol_wwK(iA) = log(result.capital_util.avg);
end
% pol_output = pol_output - mean(pol_output);
pol_output = pol_output - lr_probs * pol_output';
% pol_wwK = pol_wwK - mean(pol_wwK);
pol_wwK = pol_wwK - lr_probs * pol_wwK';

cdfs = cumsum(lr_probs);
rng(1) % seed 0 already used in simulation
initial_draw = rand(1,N);
particles = zeros(T, N);
llhs = zeros(T,1);
particles(1, initial_draw <= cdfs(1)) = 1;
for ii = 2:length(lr_probs)
    particles(1, initial_draw > cdfs(ii-1) & initial_draw <= cdfs(ii)) = ii; 
end
% define distributions of measurement equation error
lh_output = @(t, states) normpdf( ... 
    data(t,:,1), pol_output(states), p.sigma_output);
lh_capacU = @(t, states) normpdf( ... 
    data(t,:,2), pol_capacU(states), p.sigma_capacU);
lh_wwK = @(t, states) normpdf( ... 
    data(t,:,3), pol_wwK(states), p.sigma_wwK);

% log-likelihood under initial distribution
llhs(1) = log( mean( exp( ...
        log( lh_output(1, particles(1,:)) ) + ...
        log( lh_capacU(1, particles(1,:))) + ...
        log( lh_wwK(1, particles(1,:))) ...
        ) ) );

% predict, filter, update particles and collect the log-likelihood
cdfs_lt = cumsum(pi.A,2);
cdfs_gt = [ zeros(pp.NA, 1), cdfs_lt(:,1:end-1) ];

for t = 2:T
    %%% Prediction:
    shocks = rand(1,N);
    for fromA = 1:pp.NA
        for toA = 1:pp.NA
            ind = ( (particles(t-1,:) == fromA) & ...
                (shocks < cdfs_lt(fromA, toA) & shocks > cdfs_gt(fromA, toA)) );
            particles(t, ind) = toA;
        end
    end
    
    %%% Filtering:
    llh = ( log( lh_output(t, particles(t,:)) ) + ...
        log( lh_capacU(t, particles(t,:))) + ...
        log( lh_wwK(t, particles(t,:))) );
    lh = exp(llh);
    if sum(lh) > 0 
        weights = exp( llh - log( sum(lh) ) );
    else 
        weights = zeros(1, N) + 1 / N; % avoid weights being undefined
    end
    % store the log(mean likelihood)
    llhs(t) = log(mean(lh));
    
    %%% Sampling:
    particles(t,:) = datasample(particles(t,:), N, 'Weights', weights);
    
end

LLH = sum(llhs);