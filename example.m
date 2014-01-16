p.alpha = 1/3;
p.sigma = 4;
p.delta = 0.1;
p.mu_b = -1/2;
p.sig_b = 1; % mu_b = -.5 and sig_b = 1 => E[b] = 1

p.w = 0.6;
p.R = 1.04;
p.P = 1;
p.Y = 1;

p.cbar = 1.5;

p.PY = (p.P ^ (p.sigma - 1) * p.Y);

z = 1;

[x_opt, max_f] = fminsearch(@(x) -exp_profits(z,x(1),x(2),p), [1 1]);
max_f = -max_f;
k_opt = x_opt(1);
pr_opt = x_opt(2);

kv = linspace(0,5,100);
plot(kv, arrayfun(@(x) exp_profits(z,x,pr_opt,p), kv))
figure
prv = linspace(0,2,100);
plot(prv, arrayfun(@(x) exp_profits(z,k_opt,x,p), prv))


