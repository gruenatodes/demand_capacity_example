function [ res ] = simulation( p, pp, agg, z, mu_shock )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x_opt, ~] = fminsearch(@(x) ...
    -expected_profits(x(1), x(2), z, p.mu_b, p, agg), [1 1]);
% max_f = -max_f;
k_opt = x_opt(1);
pr_opt = x_opt(2);

big_brackets = ( ((2 - 2 * p.alpha) / p.alpha) ^ (p.alpha / (2 - p.alpha)) + ...
    (p.alpha / (2 - 2 * p.alpha)) ^ ((2 - 2 * p.alpha) / (2 - p.alpha)) );
b_bar = pr_opt ^ ((2 - p.alpha + p.alpha * p.sigma) / p.alpha) * ...
              agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha ) * ...
              z ^ (2 / p.alpha) * ...
              k_opt / p.chi * ...
              1 / ( agg.P ^ (p.sigma - 1) * agg.Y ) * ...
              ( (2 - p.alpha) / 2) ^ ((2 - p.alpha) / p.alpha) * ...
              big_brackets ^ (- (2 - p.alpha) / p.alpha);
max_prod = pr_opt ^ ( (2 - p.alpha) / p.alpha) * ...
        z ^ (2 / p.alpha) * ...
        agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha) * ...
        k_opt / p.chi * ...
        big_brackets ^ ( -(2 - p.alpha) / p.alpha ) * ...
        ((2 - p.alpha) / 2) ^ ((2 - p.alpha) / p.alpha);

% kv = linspace(0,5,100); kv(1) = [];
% plot(kv, arrayfun(@(x) exp_profits(z,x,pr_opt,p), kv))
% figure
% prv = linspace(0,3,100); prv(1) = [];
% plot(prv, arrayfun(@(x) exp_profits(z,k_opt,x,p), prv))

% generate lognormally distributed random numbers
rng(0) % fix the seed of random number generator
mu_with_shock = log(mu_shock) + p.mu_b; % this shifts the mean of
% the demand shock around such that E[b] = mu_shock instead of the
% expected E[b] = 1
b_draws = lognrnd(mu_with_shock, p.sigma_b, [100000, 1]);

% realized profits, labor, and production:
[rd_profit, rd_labor, rd_prod, rd_capital_util] = ...
    realized_profit(b_draws, k_opt, pr_opt, z, p, agg);

% share of firms at capacity constraint
F_bbar = logncdf(b_bar, mu_with_shock, p.sigma_b);

%%% gather statistics
% production:
prod_field = get_stats(rd_prod, rd_prod, b_draws');
% labor:
% labor_field = get_stats(rd_labor, rd_prod, b_draws');
% profits:
% profit_field = get_stats(rd_profit, rd_prod, b_draws');
% capital utilization:
capital_util_field = get_stats(rd_capital_util, rd_prod, b_draws');

% aggregate production:
b_weights = b_draws ./ sum(b_draws);
agg_output = (sum( (rd_prod(:).^ pp.sigexp) .* b_weights )) ^ (1/pp.sigexp);

% store return values in results field
res.k_opt = k_opt;
res.pr_opt = pr_opt;
res.bbar = b_bar;
res.share_unconstrained = F_bbar;
res.output = agg_output;
res.max_ys = max_prod;
res.prod = prod_field;
% res.labor = labor_field;
% res.profits = profit_field;
res.capital_util = capital_util_field;

end

