Nk = 100;
Npr = 120;
k_plot = linspace(0.01, 1, Nk);
pr_plot = linspace(0.01, 4, Npr);

z_plot = 1;

E_profits_plot = zeros(Nk, Npr);

for ik = 1:Nk
    for ipr = 1:Npr
        E_profits_plot(ik, ipr) = ...
            expected_profits(k_plot(ik), pr_plot(ipr), z, 0, p, agg);
    end
end

surf(pr_plot, k_plot, E_profits_plot)