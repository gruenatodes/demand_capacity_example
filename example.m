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
k_opt = x_opt(1)
pr_opt = x_opt(2)

kv = linspace(0,5,100); kv(1) = [];
plot(kv, arrayfun(@(x) exp_profits(z,x,pr_opt,p), kv))
figure
prv = linspace(0,2,100); prv(1) = [];
plot(prv, arrayfun(@(x) exp_profits(z,k_opt,x,p), prv))

% generate lognormally distributed random numbers
b_draws = lognrnd(p.mu_b, p.sig_b, [1000000, 1]);

% realized profits, labor, and production:
[rd_prof, rd_labor, rd_prod] = ...
    realized_profit(b_draws,z,k_opt,pr_opt,p);
disp(mean(rd_prof) - max_f)

%%% gather statistics
% share of firms at capacity constraint
[~, Fbbar, bbar] = exp_profits(z,k_opt,pr_opt,p);
disp([num2str(round((1-Fbbar)*1000)/10), '% of firms capacity constrained'])
% production at full capacity, capacity utilization (measured in production):
[~, ~, cap_prod] = realized_profit(bbar,z,k_opt,pr_opt,p);
mean_prod = mean(rd_prod);
wt_mean_prod = sum(rd_prod .* (rd_prod / sum(rd_prod)));
disp(['Maximum production is ', num2str(cap_prod)])
disp(['Average production by firm is ', ...
    num2str(mean_prod), ' (', num2str(round(1000*mean_prod/cap_prod)/10), '%)'])
disp(['Average production weighted by demand is ',...
    num2str(wt_mean_prod), ' (', num2str(round(1000*wt_mean_prod/cap_prod)/10), '%).']);
% First and second moments for marginal costs. Marginal costs under
% quasi-fixed capital are determined only by labor input, which in turn,
% given prices, depends only on the demand shock

% First and second moments for average costs (includes costs of capital).
% This measures the extent of misallocation that is due to quasi-fixed
% capital (-> suboptimal input factor ratio for a given level of
% production)







