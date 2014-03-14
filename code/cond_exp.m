function [F_xbar, E_x_gt_xbar, E_x_lt_xbar] = cond_exp(xbar, mu, sig)
% computes the conditional expectations E[x|x > xbar] and E[x|x < xbar]
% where x is distributed lognormally and mu and sig are mean and standard
% deviation of the associated normal distribution.
% xbar can be an array

if any(xbar<0)
    error('no negative xbar allowed');
end

xbar_0 = (xbar==0);
F_xbar(xbar_0) = 0;
E_x_gt_xbar(xbar_0) = exp(sig^2/2 + mu);
E_x_lt_xbar(xbar_0) = 0;

F_xbar(~xbar_0) = logncdf(xbar(~xbar_0), mu, sig);
E_x_gt_xbar(~xbar_0) = exp(sig^2/2 + mu) * ...
    ( 1 + erf( (sig^2 + mu - log(xbar(~xbar_0))) ./ (sqrt(2)*sig) ) ) ./ ...
    ( 2 * (1 - F_xbar(~xbar_0)) );
E_x_lt_xbar(~xbar_0) = exp(sig^2/2 + mu) * ...
    ( 1 - erf( (sig^2 + mu - log(xbar(~xbar_0))) ./ (sqrt(2)*sig) ) ) ./ ...
    ( 2 * F_xbar(~xbar_0) );

% E_x_lt_xbar(~xbar_0) = (exp(mu + sig^2/2) - (1-F_xbar(~xbar_0)).*E_x_gt_xbar(~xbar_0)) ./ F_xbar(~xbar_0);
