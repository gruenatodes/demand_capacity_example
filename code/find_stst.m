clc
%%% get parameters
[p, pp] = set_params();
big_brackets = ( ((2 - 2 * p.alpha) / p.alpha) ^ (p.alpha / (2 - p.alpha)) + ...
                (p.alpha / (2 - 2 * p.alpha)) ^ ((2 - 2 * p.alpha) / (2 - p.alpha)) );
z = 1;
agg.P = 1;

%%% discretize state space for demand shock:
[b, pi.b] = TAUCHEN(pp.Nb, p.rho_b, p.sigma_b, pp.tauchen_param_b);
b = b';
lr_dist_b = (zeros(1, pp.Nb) + 1 / pp.Nb) * (pi.b ^ 50); % long-run 
% distribution of idiosyncratic demand shocks
b = log( exp(p.rho_b * b) / mean(exp(p.rho_b * b))) / p.rho_b; % normalize 
% such that exp(p.rho_b * b) = 1
% policy functions:
pol_k = zeros(1, pp.Nb);
pol_pr = zeros(1, pp.Nb);

%%% initial guesses for K and Y 
K_guess = 1;
Y_guess = 1;
% w_guess = 2/3 * Y_guess;

%%% loop
KY_tol = 1; KY_ct = 0;
while KY_tol > 0.0001 && KY_ct < 50
%%% set implied aggregates:
agg.R = 1 / p.beta; % interest rate
agg.K = K_guess;
agg.Y = Y_guess;
% agg.w = w_guess;
agg.Inv = p.delta * agg.K;
agg.C = agg.Y - agg.Inv;
% agg.N = (agg.w / (p.varphi * agg.C ^ p.tau)) ^ (1 / p.epsi); % labor supply
agg.w = ( p.varphi * ((1 - p.alpha) * agg.Y) ^ p.epsi * agg.C ^ p.tau ) ^ ...
                (1 / (1 + p.epsi) );
agg.N = (1 - p.alpha) * agg.Y / agg.w;
% from FOC for C and n: ( varphi * v'(n) ) / u'(C) = w
% using functional forms v(n) = n^(1+eps)/(1+eps) and 
% u(c) = c^(tau-1)/(tau-1) ---> varphi * n^eps * c^tau = w
% aggregate relationship: n = (1 - alpha) * I / w = (1 - alpha) * Y * P / (A * w)
% combine: w^(1+eps) = varphi * ((1 - alpha) * Y)^eps * c^tau
agg

%%% derive policy functions for each of the states
for ib = 1:pp.Nb
    exp_demand = p.rho_b * b(ib);
    [pol_opt, exp_profit] = ...
        fminsearch(@(x) -expected_profits(x(1), x(2), z, exp_demand, p, agg), [1 2]);
    pol_k(ib) = pol_opt(1);
    pol_pr(ib) = pol_opt(2);
    exp_profit = -exp_profit; 
end

%%% find aggregate production.
outer_integral = 0;
for ib = 1:pp.Nb
    exp_demand = p.rho_b * b(ib);
    % is y^s < yhat?
    if pol_pr(ib) * z < agg.w ^ (1 - p.alpha) * (2 * p.chi) ^ p.alpha / ...
            ( p.alpha ^ p.alpha * (1 - p.alpha) ^ (1 - p.alpha) )

        b_bar_1 = pol_pr(ib) ^ ((2 - p.alpha + p.alpha * p.sigma) / p.alpha) * ...
          agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha ) * ...
          z ^ (2 / p.alpha) * ...
          pol_k(ib) / p.chi * ...
          1 / ( agg.P ^ (p.sigma - 1) * agg.Y ) * ...
          ( (2 - p.alpha) / 2) ^ ((2 - p.alpha) / p.alpha) * ...
          big_brackets ^ (- (2 - p.alpha) / p.alpha);
        [F_bbar, ~, c_e] = cond_exp(b_bar_1, exp_demand, p.sigma_b);
        pe_bbar1 = c_e * F_bbar; % \int_0^bbar b_i df(b_i)

        b_bar_1_exp = b_bar_1 ^ (1 / p.sigma );
        mu_b_exp = 1 / p.sigma * exp_demand;
        sigma_b_exp = 1 / p.sigma * p.sigma_b;
        [F_bbar_exp, c_e, ~] = cond_exp(b_bar_1_exp, mu_b_exp, sigma_b_exp);
        pe_bbar1_exp = c_e * (1 - F_bbar_exp); % \int_bbar^\infty b_i^{1/\sigma} df(b_i)

        inner_integral = (agg.P ^ (p.sigma - 1) * agg.Y / ...
            pol_pr(ib) ^ p.sigma) ^ ((p.sigma - 1) / p.sigma) * pe_bbar1 + ...
            (pol_pr(ib) ^ ((2 - p.alpha) / p.alpha) * ...
             agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha ) * ...
             z ^ (2 / p.alpha) * pol_k(ib) / z * ...
             ((2 - p.alpha) / 2) ^ ((2 - p.alpha) / p.alpha) * ...
             big_brackets ^ ( -(2 - p.alpha) / p.alpha) ) ^ ...
             ( (p.sigma - 1) / p.sigma) * ...
            pe_bbar1_exp;
    else
        b_bar_2 = pol_pr(ib) ^ ((1 - p.alpha + p.alpha * p.sigma) / p.alpha) / ...
            (agg.P ^ (p.sigma - 1) * agg.Y) * ...
            ( (1 - p.alpha) / agg.w ) ^ ( (1 - p.alpha) / p.alpha) * ...
            z ^ (1 / p.alpha) * pol_k(ib);
        [F_bbar, ~, c_e] = cond_exp(b_bar_2, exp_demand, p.sigma_b);
        pe_bbar2 = c_e * F_bbar;

        bbar2_exp = b_bar_2 ^ (1 / p.sigma);
        mu_b_exp = 1 / p.sigma * exp_demand;
        sigma_b_exp = 1 / p.sigma * p.sigma_b;
        [F_bbar2_exp, c_e, ~] = cond_exp(bbar2_exp, mu_b_exp, sigma_b_exp);
        pe_bbar2_exp = c_e * (1 - F_bbar2_exp);

        inner_integral = (agg.P ^ (p.sigma - 1) * agg.Y / ...
            pol_pr(ib) ^ p.sigma) ^ ((p.sigma - 1) / p.sigma) * pe_bbar2 + ...
            ((pol_pr(ib) / agg.w * (1 - p.alpha)) ^ ((1 - p.alpha) / p.alpha) * ...
             z ^ (1 / p.alpha) * pol_k(ib)) ^ ((p.sigma - 1) / p.sigma) * ...
            pe_bbar2_exp;
    end
    outer_integral = outer_integral + lr_dist_b(ib) * inner_integral;
end
aggregate_prod = outer_integral ^ (p.sigma / (p.sigma - 1));
aggregate_capital = pol_k * lr_dist_b';
% aggregate_labordemand = (1 - p.alpha) * agg.Y / agg.w;

%%% compare production and capital demand with guesses
tol_Y = agg.Y - aggregate_prod;
tol_K = agg.K - aggregate_capital;
% tol_excess_N_demand = aggregate_labordemand - agg.N;
KY_tol = max(abs(tol_Y), abs(tol_K))
KY_ct = KY_ct + 1;

[aggregate_prod, aggregate_capital]

adj_speed_K = 1/3;
adj_speed_Y = 1/3;
% adj_speed_w = 1/10;
K_guess = (1 - adj_speed_K) * agg.K + adj_speed_K * aggregate_capital
Y_guess = (1 - adj_speed_Y) * agg.Y + adj_speed_Y * aggregate_prod
% w_guess = agg.w + adj_speed_w * tol_excess_N_demand

end


