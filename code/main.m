% set parameters
[p, pp] = set_params();
big_brackets = ( ((2 - 2 * p.alpha) / p.alpha) ^ (p.alpha / (2 - p.alpha)) + ...
                (p.alpha / (2 - 2 * p.alpha)) ^ ((2 - 2 * p.alpha) / (2 - p.alpha)) );
z = 1;

% discretize state space for demand shock:
[b, pi.b] = TAUCHEN(pp.Nb, p.rho_b, p.sigma_b, pp.tauchen_param_b);
b = b';
lr_dist_b = (zeros(1, pp.Nb) + 1 / pp.Nb) * (pi.b ^ 50); % long-run 
% distribution of idiosyncratic demand shocks

% initial guesses for aggregate variables
agg.Y = 0.5;
agg.w = 2/3 * agg.Y;
agg.R = 1.04;
agg.P = 1;

pol_k = zeros(1, pp.Nb);
pol_pr = zeros(1, pp.Nb);

tol = 1; count = 0;
while tol > 0.001 % && count < 150
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
    
    %%% loop control
    tol = abs(agg.Y - aggregate_prod);
    agg.Y = 2/3 * agg.Y + 1/3 * aggregate_prod;
    agg.w = 2/3 * agg.Y;
    count = count + 1;
end

agg.Y
aggregate_prod
[pol_k; pol_pr; exp(p.rho_b * b) ./ pol_k]

% put consistency test for aggregate price index



