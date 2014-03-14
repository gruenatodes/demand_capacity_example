function [ p, pp ] = set_params( )
%set_params Set economic and numerical parameters
% Output:
% p: economic parameters
% pp: program parameters (discretization, tolerances, etc)

% economic parameters
p.alpha = 1/3;
p.sigma = 4; % elasticity of substitution between interm goods
p.delta = 0.1;
p.sigma_b = 0.3; % sd of shocks to idiosyncratic demand
p.mu_b = -(p.sigma_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2
p.rho_b = 0.9; % AR-1 persistence of idiosyncratic demand

p.chi = 1.5; % parameter of capital utilization costs

% program parameters
pp.sigexp = (p.sigma - 1) / p.sigma; % short for exponent in DS function

pp.tauchen_param_b = 2; % how many sds to span the discretization over

pp.Nb = 11;


end

