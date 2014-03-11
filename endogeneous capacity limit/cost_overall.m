function [cost] = cost_overall(y, z, k, p, agg)
cutoff_k = y / z * ...
    ( p.alpha / (1 - p.alpha) * agg.w / (2 * p.chi) ) ^ (1 - p.alpha);
if k > cutoff_k
    cost = cost_u(y, z, k, p, agg);
else
    cost = cost_c(y, z, k, p, agg);
end
end