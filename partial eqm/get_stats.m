function [ res_field ] = get_stats( values, prod_weight, demand_weight )
%get_stats returns weighted averages of means and dispersion measures

% normalize weights so they sum to 1:
prod_weight = prod_weight / sum(prod_weight);
demand_weight = demand_weight / sum(demand_weight);
flat_weight = ones(size(values)) / length(values);

%%% means:
% simple average:
res_field.avg = mean(values);
% geometric average:
try
    res_field.geomn = geomean(values);
catch err
    if strcmp(err.identifier, 'stats:geomean:BadData') 
        res_field.geomn = 'not set, data contained negative values';
    else
        rethrow(err);
    end
end
% mean weighted by production:
res_field.avg_prod = sum(values .* prod_weight);
% mean weighted by demand:
res_field.avg_demand = sum(values .* demand_weight);

%%% dispersion measures:
% simple coefficient of variation:
res_field.cv = coeff_var(values, flat_weight);
% cv weighted by demand:
res_field.cv_demand = coeff_var(values, demand_weight);

end

