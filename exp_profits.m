function [ profit, F_bbar, bbar ] = exp_profits( z,mu_b,k,pr,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if k <= 0 || pr <= 0
    profit = 0;
    return
end

bbar = z .* k .* p.cbar ^ (1 - p.alpha) .* pr .^ p.sigma / p.PY;

[F_bbar, ~, c_e] = cond_exp(bbar, mu_b, p.sig_b);
c_e = c_e * F_bbar;
% partial expectation E[x | x < xbar] * F(xbar)

[F_bb_exp, ~, c_e_exp] = cond_exp(bbar ^ (1/(1-p.alpha)), ...
    mu_b/(1-p.alpha), p.sig_b/(1-p.alpha));
c_e_exp = c_e_exp * F_bb_exp;
% partial expectation E[x^(1/(1-alpha)) | x^(1/(1-alpha)) < xbar] *
% F(xbar)

profit = ((1 - p.delta) - p.R) * k + ...
    p.PY / pr .^ (p.sigma - 1) * c_e - ...
    (p.PY / (pr .^ p.sigma .* z .* k .^ p.alpha)) .^ (1 / (1-p.alpha)) .* ...
        p.w .* c_e_exp + ...
    (1 - F_bbar) * (pr .* z .* k .* p.cbar ^ (1 - p.alpha) - p.w * p.cbar * k );

end

