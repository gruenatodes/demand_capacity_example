function [ res ] = simulation( p, z, mu_shock )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x_opt, max_f] = fminsearch(@(x) -exp_profits(z,p.mu_b,x(1),x(2),p), [1 1]);
max_f = -max_f;
k_opt = x_opt(1)
pr_opt = x_opt(2)

% kv = linspace(0,5,100); kv(1) = [];
% plot(kv, arrayfun(@(x) exp_profits(z,x,pr_opt,p), kv))
% figure
% prv = linspace(0,3,100); prv(1) = [];
% plot(prv, arrayfun(@(x) exp_profits(z,k_opt,x,p), prv))

% generate lognormally distributed random numbers
rng(0) % fix the seed of random number generator
mu_with_shock = log(mu_shock) - p.sig_b ^ 2 / 2; % this shifts the mean of
% the demand shock around such that E[b] = mu_shock instead of the
% expected E[b] = 1
b_draws = lognrnd(mu_with_shock, p.sig_b, [1000000, 1]);

% realized profits, labor, and production:
[rd_prof, rd_labor, rd_prod] = ...
    realized_profit(b_draws,z,k_opt,pr_opt,p);
l_k_ratio = rd_labor ./ k_opt;
% note that the l/k-ratio at the cap is just cbar, by definition.

% share of firms at capacity constraint
[~, Fbbar, bbar] = exp_profits(z,mu_with_shock,k_opt,pr_opt,p);
disp([num2str(round((1-Fbbar)*1000)/10), '% of firms capacity constrained'])

%%% gather statistics
% production:
prod_field = get_stats(rd_prod, rd_prod, b_draws);
[cap_profit, ~, cap_prod] = realized_profit(bbar,z,k_opt,pr_opt,p);
prod_field.cap = cap_prod;
% labor:
labor_field = get_stats(rd_labor, rd_prod, b_draws);
cbar = ( ( ( 1 - p.alpha ) * pr_opt * z ) / p.w ) .^ (1 / p.alpha);
labor_field.cap = cbar * k_opt;
% profits:
profit_field = get_stats(rd_prof, rd_prod, b_draws);
profit_field.cap = cap_profit;
% l/k ratio
l_k_field = get_stats(rd_prof, rd_prod, b_draws);
l_k_field.cap = cbar;
% Marginal costs under quasi-fixed capital (and constant z) are determined 
% only by labor input, which in turn, given prices, depends only on the 
% demand shock
MC_fun = @(l_k) 1/z * p.w / (1 - p.alpha) * l_k .^ (p.alpha);
MC = MC_fun(l_k_ratio);
MC_field = get_stats(MC, rd_prod, b_draws);
MC_field.cap = MC_fun(cbar);
% Average costs, including costs of capital.
% This measures the extent of misallocation that is due to quasi-fixed
% capital (-> suboptimal input factor ratio for a given level of
% production) and pre-set prices
AC_fun = @(l_k) 1/z * ( (p.R - 1 + p.delta) * l_k .^ (p.alpha - 1) + ...
            p.w * l_k .^ (p.alpha) );
AC = AC_fun(l_k_ratio);
AC_field = get_stats(AC, rd_prod, b_draws);
AC_field.cap = AC_fun(cbar);

l_k_vec = linspace(0,cbar, 101); l_k_vec(1:4) = [];
vec2 = linspace(cbar, 2*cbar, 201); vec2(1:2) = [];
figure, hold on
plot(l_k_vec, arrayfun(AC_fun, l_k_vec), vec2, arrayfun(AC_fun, vec2))
line([cbar cbar],[0.7 0.9])
dots = l_k_ratio(10000:10000:end); ACdots = AC(10000:10000:end);
scatter(dots(ACdots<2), ACdots(ACdots<2),3)

% aggregate production:
b_weights = b_draws ./ sum(b_draws);
agg_output = (sum( (rd_prod(:).^ p.sigexp) .* b_weights )) ^ (1/p.sigexp);

% store return values in results field
res.k_opt = k_opt;
res.pr_opt = pr_opt;
res.bbar = bbar;
res.share_unconstrained = Fbbar;
res.output = agg_output;
res.prod = prod_field;
res.labor = labor_field;
res.profits = profit_field;
res.l_k = l_k_field;
res.MC = MC_field;
res.AC = AC_field;


end

