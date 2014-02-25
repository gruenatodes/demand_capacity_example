clear, clc

% load parameters
[p, pp] = set_params(); % p = economic parameters, pp = program parameters

% discretize state space
[b, pi.b] = TAUCHEN(pp.Nb, p.rho_b, p.sigma_b, pp.tauchen_param_b);
b = exp(b);
sig_b = zeros(1, pp.Nb) + p. sigma_b; % could adjust this to keep absolute variance constant
z = 1; % could think about TFP shock, not considered for now.

% since policy doesn't affect evolution of the state space, we have a
% sequence of 1-shot maximization problems, so can just analyze this static
% problem
agg.w = 0.6;
agg.R = 1.04;
agg.P = 1;
agg.Y = 1;
agg.PY = (agg.P ^ (p.sigma - 1) * agg.Y); % short for agg demand factors

k_opt = zeros(1, pp.Nb);
pr_opt = zeros(1, pp.Nb);
for ib = 1:pp.Nb
    [x_opt, max_profit] = fminsearch(@(x) ...
        -exp_profits(z, [b(ib), sig_b(ib)], x(1), x(2), agg, p), [1 1]);
    max_profit = -max_profit;
    k_opt(ib) = x_opt(1);
    pr_opt(ib) = x_opt(2);
end

