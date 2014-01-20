p.alpha = 1/3;
p.sigma = 4; p.sigexp = (p.sigma - 1) / p.sigma;
p.delta = 0.1;
p.sig_b = 1; 
p.mu_b = -(p.sig_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2

p.w = 0.6;
p.R = 1.04;
p.P = 1;
p.Y = 1;

p.cbar = 1.5; % maximum ratio of labor to capital: l/k <= cbar

p.PY = (p.P ^ (p.sigma - 1) * p.Y); % short for agg demand factors

z = 1;

demand_shocks = exp([-0.05, 0, 0.05])

close all
results_stst = simulation(p,z,demand_shocks(2))

results_low = simulation(p,z,demand_shocks(1))

results_high = simulation(p,z,demand_shocks(3))

