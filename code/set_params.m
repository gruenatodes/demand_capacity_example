function [ p, pp ] = set_params( )
%set_params Set economic and numerical parameters
% Output:
% p: economic parameters
% pp: program parameters (discretization, tolerances, etc)

%%% economic parameters
% technology:
p.alpha = 1/3; % production elasticity of capital
p.sigma = 4; % elasticity of substitution between interm goods
p.delta = 0.1; % depreciation rate
% household preferences:
p.beta = 0.96; % rate of time preference
p.phi = 1; % stst value of shock to rate of time preference
p.tau = 2; % EIS: u(c) = c^(1 - tau) / (1 - tau)
p.epsi = 1/3; % Inv. of Frisch elasticity of labor: v(n) = n^(1+epsi) / (1+epsi)
p.varphi = 1; % weight on disutility of labor in VNM-utility u(C) - varphi * v(n)
% evolution of demand:
p.sigma_b = 0.3; % sd of shocks to idiosyncratic demand
p.mu_b = -(p.sigma_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2
p.rho_b = 0.9; % AR-1 persistence of idiosyncratic demand

p.chi = .5; % parameter of capital utilization costs

% program parameters
pp.sigexp = (p.sigma - 1) / p.sigma; % short for exponent in DS function

pp.tauchen_param_b = 2; % how many sds to span the discretization over

pp.Nb = 11;


end

