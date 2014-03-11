p.alpha = 1/3;
p.sigma = 8; p.sigexp = (p.sigma - 1) / p.sigma;
p.delta = 0.1;
p.sig_b = 0.75; 
p.mu_b = -(p.sig_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2

p.w = 0.6;
p.R = 1.04;
p.P = 1;
p.Y = 1;

% p.cbar = 0.75; % maximum ratio of labor to capital: l/k <= cbar

p.PY = (p.P ^ (p.sigma - 1) * p.Y); % short for agg demand factors

z = 1;

shock_dev = 0.1;
demand_shocks = exp([-shock_dev, 0, shock_dev])

close all
r.stst = simulation(p,z,demand_shocks(2));
r.stst
title('agg steady state')

r.low = simulation(p,z,demand_shocks(1));
r.low
title('low shock')

r.high = simulation(p,z,demand_shocks(3));
r.high
title('high shock')

