clear;
close all;
tau0 = 0.0;
H = 10;
h = 5;
z = - linspace(0, 1, 1000) * h;
d = -z;
sig = d / h;

kpp = KPP();

% This subplot should look like LMD94 Figure 2. scalar

suptitle("Reproduce LMD94 Figure 2 scalar case (dashed line)");

ax = {};

for i = 1:3
    ax{i} = subplot(1, 3, i);
    hold(ax{i}, 'on');
end

title(ax{1}, "G");
title(ax{2}, "w_s");
title(ax{3}, "w_m");
plot(ax{1}, kpp.c.calG(sig), sig);
h_over_L_arr = [1 0.1 0 -1 -5];
legend_str = {};
for i=1:length(h_over_L_arr)
    h_over_L = h_over_L_arr(i);
    L_star_target = h / h_over_L;
    [ u_star, ~ ] = kpp.calMOSTscales(tau0, 1);
    wb_sfc = u_star^3 / ( - kpp.c.kappa * L_star_target );
    [ w_s, ~, u_star, L_star ] = kpp.calw_x(kpp.c.SCALAR, z, h, tau0, wb_sfc);
    [ w_x, ~, u_star, L_star ] = kpp.calw_x(kpp.c.MOMENTUM, z, h, tau0, wb_sfc);

    plot(ax{2}, w_s ./ (kpp.c.kappa * u_star), sig, '--');
    plot(ax{3}, w_x ./ (kpp.c.kappa * u_star), sig, '-');
    
    legend_str{end+1} = sprintf('h/L_* = %.1f', h_over_L);
    fprintf('h_over_L(%d) = %f\n', i, h / L_star);
end
xlim(ax{2}, [0 4]);    
xlim(ax{3}, [0 4]);    

legend(ax{3}, legend_str, 'Location', 'south');

for i = 1:length(ax)
    set(ax{i}, 'Ydir', 'reverse');
    hold(ax{i}, 'off');
end
