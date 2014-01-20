function CV = coeff_var(values, weights)
% values: 
% weights: have to sum to 1; usually a distribution (pdf)
values = values(weights>0);
weights = weights(weights>0);
CV = sqrt(var(values, weights)) / (values' * weights);