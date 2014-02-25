function [ p, pp ] = set_params( )

% economic parameters
p.alpha = 1/3;
p.sigma = 4; % elasticity of substitution between interm goods
p.delta = 0.1;
p.sigma_b = 0.1; % sd of shocks to idiosyncratic demand
p.mu_b = -(p.sigma_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2
p.rho_b = 0.9; % AR-1 persistence of idiosyncratic demand

p.cbar = 0.75; % maximum ratio of labor to capital: l/k <= cbar

% program parameters
pp.sigexp = (p.sigma - 1) / p.sigma; % short for exponent in DS function

pp.tauchen_param_b = 2; % how many sds to span the discretization over

pp.Nk = 100;
pp.Npr = 100;
pp.Nb = 11;


