function [exp_profits] = expected_profits(k, pr, z, e_b, p, agg)

if k <= 0 || pr <= 0
    exp_profits = - 1000000;
    return
end

full_capacity = 0;
if pr * z > agg.w ^ (1 - p.alpha) * (2 * p.chi) ^ p.alpha / ...
        ( p.alpha ^ p.alpha * (1 - p.alpha) ^ (1 - p.alpha) )
    full_capacity = 1;
end

big_brackets = ( ((2 - 2 * p.alpha) / p.alpha) ^ (p.alpha / (2 - p.alpha)) + ...
                (p.alpha / (2 - 2 * p.alpha)) ^ ((2 - 2 * p.alpha) / (2 - p.alpha)) );

if ~full_capacity
    b_bar_1 = pr ^ ((2 - p.alpha + p.alpha * p.sigma) / p.alpha) * ...
              agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha ) * ...
              z ^ (2 / p.alpha) * ...
              k / p.chi * ...
              1 / ( agg.P ^ (p.sigma - 1) * agg.Y ) * ...
              ( (2 - p.alpha) / 2) ^ ((2 - p.alpha) / p.alpha) * ...
              big_brackets ^ (- (2 - p.alpha) / p.alpha);
    [F_bbar, ~, c_e] = cond_exp(b_bar_1, e_b, p.sigma_b);
    pe_bbar1 = c_e * F_bbar; % \int_0^bbar b_i df(b_i)
    
    b_bar_1_exp = b_bar_1 ^ (2 / (2 - p.alpha) );
    mu_b_exp = (2 / (2 - p.alpha) * e_b);
    sigma_b_exp = (2 / (2 - p.alpha) * p.sigma_b);
    [F_bbar_exp, ~, c_e] = cond_exp(b_bar_1_exp, mu_b_exp, sigma_b_exp);
    pe_bbar1_exp = c_e * F_bbar_exp; % \int_0^bbar b_i^{1/(2-\alpha} df(b_i)
    
    exp_profits = (agg.P / pr) ^ (p.sigma - 1) * agg.Y * pe_bbar1 - ...
        agg.w ^ ( (2 - 2 * p.alpha) / (2 - p.alpha)) * ...
            ( agg.P ^ (p.sigma - 1) * agg.Y / (z * pr ^ p.sigma) ) ^ ...
                ( 2 / (2 - p.alpha) ) * ...
            (p.chi / k) ^ ( p.alpha / (2 - p.alpha) ) * ...
            big_brackets * pe_bbar1_exp + ...
        (1 - F_bbar) * (pr * z) ^ (2 / p.alpha) * ...
            agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha) * ...
            k / p.chi * ...
            big_brackets ^ ( -(2 - p.alpha) / p.alpha ) * ...
            ((2 - p.alpha) / 2) ^ (2 / p.alpha) * p.alpha / (2 - p.alpha);
else
    bhat = pr ^ p.sigma / (agg.P ^ (p.sigma - 1) * agg.Y) * z * k * ...
        (2 * p.chi * (1 - p.alpha) / (p.alpha * agg.w)) ^ (1 - p.alpha);
    b_bar_2 = pr ^ ((1 - p.alpha + p.alpha * p.sigma) / p.alpha) / ...
        (agg.P ^ (p.sigma - 1) * agg.Y) * ...
        ( (1 - p.alpha) / agg.w ) ^ ( (1 - p.alpha) / p.alpha) * ...
        z ^ (1 / p.alpha) * k;
    
    [F_bhat, ~, c_e] = cond_exp(bhat, e_b, p.sigma_b);
    pe_bhat = c_e * F_bhat;
    
    bhat_exp1 = bhat ^ (2 / (2 - p.alpha) );
    mu_b_exp = (2 / (2 - p.alpha)) * e_b;
    sigma_b_exp = (2 / (2 - p.alpha)) * p.sigma_b;
    [F_bhat_exp1, ~, c_e] = cond_exp(bhat_exp1, mu_b_exp, sigma_b_exp);
    pe_bhat_exp1 = c_e * F_bhat_exp1;
    
    [F_bbar2, ~, c_e] = cond_exp(b_bar_2, e_b, p.sigma_b);
    pe_bbar2 = c_e * F_bbar2;
    
    bbar2_exp = b_bar_2 ^ (1 / (1 - p.alpha) );
    mu_b_exp = (1 / (1 - p.alpha)) * e_b;
    sigma_b_exp = (1 / (1 - p.alpha)) * p.sigma_b;
    [F_bbar2_exp, ~, c_e] = cond_exp(bbar2_exp, mu_b_exp, sigma_b_exp);
    pe_bbar2_exp = c_e * F_bbar2_exp;
    
    bhat_exp2 = bhat ^ (1 / (1 - p.alpha) );
    mu_b_exp = (1 / (1 - p.alpha)) * e_b;
    sigma_b_exp = (1 / (1 - p.alpha)) * p.sigma_b;
    [F_bhat_exp2, ~, c_e] = cond_exp(bhat_exp2, mu_b_exp, sigma_b_exp);
    pe_bhat_exp2 = c_e * F_bhat_exp2;
    
    exp_profits = (agg.P / pr) ^ (p.sigma - 1) * agg.Y * pe_bhat - ...
        agg.w ^ ( (2 - 2 * p.alpha) / (2 - p.alpha) ) * ...
            ( agg.P ^ (p.sigma - 1) * agg.Y / (z * pr ^ p.sigma) ) ^ ...
                (2 / (2 - p.alpha)) * ...
            (p.chi / k) ^ (p.alpha / (2 - p.alpha)) * ...
            big_brackets * pe_bhat_exp1 + ...
        (agg.P / pr) ^ (p.sigma - 1) * agg.Y * (pe_bbar2 - pe_bhat) - ...
        ( agg.P ^ (p.sigma - 1) * agg.Y  / (z * pr ^ p.sigma) ) ^ (1 / (1 - p.alpha)) * ...
            agg.w / (k ^ (p.alpha / (1 - p.alpha))) * ...
            (pe_bbar2_exp - pe_bhat_exp2) + ...
        (1 - F_bbar2) * (pr*z) ^ (1 / p.alpha) * ...
            ( (1 - p.alpha) / agg.w ) ^ ((1 - p.alpha) / p.alpha) * ...
            k * p.alpha - ...
        (1 - F_bhat) * p.chi * k;
    
end

exp_profits = ((1 - p.delta) - agg.R) * k + exp_profits;

end
