close all;

time_ticks = (t_sim(1)):1:(t_sim(end));
time_labels = {};
empty_labels = {};
for i=1:length(time_ticks)
    time_labels{i} = datestr(time_ticks(i), "mm/dd hh:MM");
    empty_labels{i} = '';
end

%%%%% plot forcings %%%%%%
figure;
figsize = [ 12 16 ];
set(gcf, 'Units', 'Inches', 'Position', [0, 0, figsize(1), figsize(2)], 'PaperUnits', 'Inches', 'PaperSize', [figsize(1), figsize(2)])

ax = [];
N = 6;
for i=1:N
    ax(i) = subplot(N, 1, i);
end

axes(ax(1));
hold(gca, 'on');
plot(gca, t_sim, rec.T_a, 'r-');
plot(gca, t_sim_mid, forcing.T(1, :), 'k-');
plot(gca, t_sim, rec.T(1, :), 'k--');
title(gca, 'SST: observation (solid black), simulation (dashed black); T_a (red)');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel(gca, '[ degC ]');
ylim(gca, [20 32]);
hold(gca, 'off');

axes(ax(2));
hold(gca, 'on');
plot(gca, t_sim_mid, forcing.S(1, :), 'k-');
plot(gca, t_sim, rec.S(1, :), 'k--');
title(gca, 'SSS: observation (solid), simulation (dashed)');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel(gca, '[ PSU ]');
ylim(gca, [33 34]);
hold(gca, 'off');

axes(ax(3));
hold(gca, 'on');
plot(gca, t_sim, - rec.I_0,  'r--');
title(gca, 'Left axis: I\_0 (dashed black). Right axis: Hf\_sen (black), Hf\_lat (blue), Hf\_lw (red).');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel(gca, '[ W / m^2 / s ]');
ylim(gca, [0 1200]);

yyaxis right;
plot(gca, t_sim, - rec.Hf_sen, 'k-');
plot(gca, t_sim, - rec.Hf_lat, 'b-');
plot(gca, t_sim, - rec.Hf_lw,  'r-');
ylabel(gca, '[ W / m^2 / s ]');

hold(gca, 'off');

axes(ax(4));
hold(gca, 'on');
plot(gca, t_sim, rec.precip * 3600 * 1000, 'b-');
title(gca, 'Left axis: precipitation (blue). Right axis: evaporation (red).');
set(gca, 'ycolor', 'blue');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel(gca, '[ mm / hr ]');
left_ylim = ylim(gca);

yyaxis right;
plot(t_sim, rec.evap * 3600 * 1000 * 100, 'r-');
ylabel('[ \times 10^{-2} mm / hr]');
ylim(left_ylim);
set(gca, 'ycolor', 'red');
hold(gca, 'off');

axes(ax(5));
hold(gca, 'on');
plot(gca, t_sim, sqrt(rec.U10.^2 + rec.V10.^2), 'k-');
title(gca, 'U10');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel(gca, '[ m / s ]');
hold(gca, 'off');

axes(ax(6));
hold(gca, 'on');
plot(gca, t_sim, rec.h, 'k-');
title(gca, 'Mixed-layer depth');
xticks(gca, time_ticks);
xticklabels(gca, empty_labels);
ylabel('[ m ]');
set(gca, 'YDir', 'reverse');
hold(gca, 'off');

axes(ax(end))
xticklabels(gca, time_labels);

saveas(gcf, 'forcing_and_diagnoses.png');




%%%%% Plot hovemoeller diagram %%%%%

clevel_T = linspace(25, 32, 15);
clevel_S = linspace(33, 34, 11);
figure;
suptitle("T and S");

subplot(2, 2, 1);
hold on;
[C, h] = contourf(t_sim, m.grid.z_T, rec.T, clevel_T);
set(h, 'LineColor', 'None');
plot(t_sim, - rec.h, 'w-');
colorbar();
caxis('manual');
caxis([clevel_T(1) clevel_T(end)]);
hold off;

subplot(2, 2, 2);
[C, h] = contourf(t_sim_mid, m.grid.z_T, forcing.T, clevel_T);
set(h, 'LineColor', 'None');
colorbar();
caxis('manual');
caxis([clevel_T(1) clevel_T(end)]);


subplot(2, 2, 3);
hold on;
[C, h] = contourf(t_sim, m.grid.z_T, rec.S, clevel_S);
set(h, 'LineColor', 'None');
plot(t_sim, - rec.h, 'w-');
colorbar();
caxis('manual');
caxis([clevel_S(1) clevel_S(end)]);
hold off;

subplot(2, 2, 4);
[C, h] = contourf(t_sim_mid, m.grid.z_T, forcing.S, clevel_S);
set(h, 'LineColor', 'None');
colorbar();
caxis('manual');
caxis([clevel_S(1) clevel_S(end)]);

% Compare SST
figure;
subplot(2, 1, 1);
hold on;
plot(t_sim, rec.T(1, :), 'r-');
plot(t_sim_mid, forcing.T(1, :), 'b--');
title('SST. Observation (blue dashed), simulation (red solid)');
hold off;

subplot(2, 1, 2);
hold on;
plot(t_sim, rec.S(1, :), 'r-');
plot(t_sim_mid, forcing.S(1, :), 'b--');
title('SSS. Observation (blue dashed), simulation (red solid)');
hold off;