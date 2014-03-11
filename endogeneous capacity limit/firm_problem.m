function [util_rate, labor, cost, prod, profit] = firm_problem(k, pr, z, p, agg)
% firm input choice problem given price, capital as well as structural
% parameters. 
% INPUTS:
% k - capital stock
% pr - firm's goods price
% z - firm's Cobb-Douglas 'technology' level
% p - (structure of) structural parameters
% agg - structure of aggregate variables w, P, Y
% OUTPUTS:
% util_rate - utilized capital as a fraction of k
% labor - labor hours hired
% cost - total cost of production (wage*labor + cost of utilizing capital)
% prod - CD-function of labor and capital
% profit - price*prod - cost

% cutoff_k = y / z * ...
%     ( p.alpha / (1 - p.alpha) * agg.w / (2 * p.chi) ) ^ (1 - p.alpha); % cutoff for full capital utilization
cutoff_y_lt_yhat = p.chi ^ p.alpha * agg.w ^ (1 - p.alpha) * (2 / 2 - p.alpha) * ...
    ( ((2 - 2 * p.alpha) / p.alpha) ^ p.alpha + ...
      (p.alpha / (2 - 2 * p.alpha)) ^ (1 - p.alpha) ); % is a cutoff for pr*z
cutoff_y_gt_yhat = (2 * p.chi) ^ p.alpha / ...
    ( agg.w ^ (1 + p.alpha) * p.alpha ^ p.alpha * (1 - p.alpha) ^ (1 - p.alpha) );

if pr * z > cutoff_y_gt_yhat % y > yhat
    prod = ( pr / agg.w * (1 - p.alpha) ) ^ ( (1 - p.alpha) / p.alpha ) * ...
        z ^ (1 / p.alpha) * k;
    cost = cost_c(prod, z, k, p, agg);
    k_tilde = k;
    util_rate = k_tilde / k;
else
    if pr * z < cutoff_y_lt_yhat % y < yhat
        prod = pr ^ ( (2 - p.alpha) / p.alpha ) * ...
               agg.w ^ ( -(2 - 2 * p.alpha) / p.alpha ) * ...
               z ^ (2 / p.alpha) * ...
               k / p.chi * ...
               ( (2 - p.alpha) / 2 ) ^ ( (2 - p.alpha) / p.alpha ) * ...
               ( ( ( (2 - 2 * p.alpha) / p.alpha ) ^ ( p.alpha / (2 - p.alpha) ) + ...
                 ( p.alpha / (2 - 2 * p.alpha) ) ^ ( (2 - 2 * p.alpha) / (2 - p.alpha) ) ) ) ^ ...
                 ( -(2 - p.alpha) / p.alpha );
        cost = cost_u(prod, z, k, p, agg);
        k_tilde = ( (prod / z) * ...
            (k * agg.w / p.chi * p.alpha / (2 - 2 * p.alpha)) ^ (1 - p.alpha) ) ^ ...
                (1 / (2 - p.alpha));
        util_rate = k_tilde / k;
    else % y = yhat
        prod = z * k * ...
            ( 2 * p.chi / agg.w * (1 - p.alpha) / p.alpha ) ^ (1 - p.alpha); % yhat
        cost = cost_c(prod, z, k, p, agg); % = cost_u(prod) -- can assert this
        k_tilde = k;
        util_rate = k_tilde / k;
    end
end
labor = (prod / (z * k_tilde ^ p.alpha)) ^ (1 / (1 - p.alpha));
profit = pr * prod - cost;
end

