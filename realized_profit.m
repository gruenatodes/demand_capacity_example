function [ profit, labor ] = realized_profit( b,z,k,pr,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Supports vector input for b, scalar for the other variables
if any(k <= 0) || any(pr <= 0)
    error('k and pr must be strictly positive')
end

bbar = z .* k .* p.cbar ^ (1 - p.alpha) .* pr .^ p.sigma / p.PY;
unconstr = b < bbar;

profit = zeros(size(b)) + pr * z * k * p.cbar ^ (1 - p.alpha) - ...
    p.w * k * p.cbar; % this is the constrained default (does not depend on b)
labor = ;

profit(unconstr) = b(unconstr) .* p.PY / pr .^ (p.sigma - 1) - ...
    p.w * (b(unconstr) * p.PY ./ ...
        (pr .^ p.sigma .* z .* k .^ p.alpha)) .^ (1 / (1-p.alpha));
labor(unconstr) = ;
    
% subtract capital costs:
profit = profit + ((1 - p.delta) - p.R) * k;

end

