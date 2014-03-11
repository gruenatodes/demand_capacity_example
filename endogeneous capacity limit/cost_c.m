function [cost] = cost_c(y, z, k, p, agg)
% cost function at capacity, ie k_tilde = k
cost = agg.w * (y / (z * k ^ p.alpha)) ^ (1 / (1 - p.alpha)) + p.chi;
end