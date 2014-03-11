function [ profit ] = profits_using_cost(y, z, k, pr, p, agg)

profit = pr * y - cost_overall(y,z,k,p,agg);

end

