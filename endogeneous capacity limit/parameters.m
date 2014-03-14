p.chi = 1.5;
p.alpha = 1/3;
p.delta = 0.1;

p.sigma = 6;

p.sigma_b = 0.3;
p.mu_b = -(p.sigma_b ^ 2 / 2); % E[b] = exp(mu + sig^2/2) = 1 -> mu = -sig^2/2

agg.Y = 0.95;
agg.P = 1;
agg.w = 0.6;
agg.R = 1.04;