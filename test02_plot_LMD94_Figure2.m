addpath('./KPP');

clear;

tau0 = 0.1;
H = 10;
h = 5;
z = - linspace(0, 1, 1000) * h;
d = -z;
sig = d / h;

kpp = KPPConstants();

% This subplot should look like LMD94 Figure 2. scalar

suptitle("Reproduce LMD94 Figure 2 scalar case (dashed line)");

subplot(121);
title("G (left), w_s (right)");
plot(kpp.calG(sig), sig);
set(gca, 'Ydir', 'reverse');

subplot(122);
h_over_L_arr = [1 0.1 0 -1 -5];
hold on;
legend_str = {};
for i=1:length(h_over_L_arr)
    h_over_L = h_over_L_arr(i);
    L_star_target = h / h_over_L;
    [ u_star, ~ ] = calMOSTscales(kpp, tau0, 1);
    wb_sfc = u_star^3 / ( - kpp.kappa * L_star_target );
    [ w_s, ~, u_star, L_star ] = calw_s(kpp, z, h, tau0, wb_sfc);
    plot(w_s ./ (kpp.kappa * u_star), sig, '--');
    legend_str{end+1} = sprintf('h/L_* = %.1f', h_over_L);
    fprintf('h_over_L(%d) = %f\n', i, h / L_star);
end
xlim([0 4]);
set(gca, 'Ydir', 'reverse');
legend(gca, legend_str, 'Location', 'south');
hold off;
