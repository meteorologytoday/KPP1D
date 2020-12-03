clear;

% This simulation is meant to reproduce the DC case in
% Van Roekel et al (2018) where DC stands for "Diurnal Cycle"

Kv_iso = 1e-3;
Nz = 240;
H  = 60; % m
total_time = 86400 * 5;
dt = 20*60;
f = 1e-4;
Q_0 = 75; % sensible heat loss (pos = upward)
Q_sw_max = 235.62;

T_sfc = 20; dTdz = 0.01;
S_sfc = 35; dSdz = 0;
pause_time = 0.01;

total_steps = total_time / dt;

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

% Profile A
m.state.T = T_sfc + dTdz * m.grid.z_T;
m.state.S = S_sfc + dSdz * m.grid.z_T;
m.update_b();

% Solar radiation
m.state.Hf_sen = 75;
F_sw_func = @(t) - Q_sw_max * max(0, cos(2 * pi * (t / 86400 - 0.5)));

m.state.taux0 = 0.05 * 0;

figure;

K_s = zeros(m.grid.W_pts, total_steps+1);
wT = zeros(m.grid.W_pts, total_steps+1);
wS = zeros(m.grid.W_pts, total_steps+1);
T = zeros(m.grid.T_pts, total_steps+1);
S = zeros(m.grid.T_pts, total_steps+1);

t(1) = 0;
h_sim(1) = m.state.h;
F_sw_sim(1) = F_sw_func(t(1)); 
Hf(1) = m.state.Hf_sen + m.state.Hf_lat;
T(:, 1) = m.state.T;
S(:, 1) = m.state.S;


hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    
    
    m.state.I_0 = F_sw_sim(end);    
    diag_kpp = m.stepModel();
    
    plot(m.state.T, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;
    h_sim(end+1) = m.state.h;
    F_sw_sim(end+1) = F_sw_func(t(end));    
    Hf(end+1) = m.state.Hf_sen + m.state.Hf_lat;
    K_s(:, step+1) = diag_kpp.K_s_ML + diag_kpp.K_s_INT;  
    wT(:, step+1) = diag_kpp.loc_flux_T + diag_kpp.nloc_flux_T;
    wS(:, step+1) = diag_kpp.loc_flux_S + diag_kpp.nloc_flux_S;
    T(:, step+1) = m.state.T;
    S(:, step+1) = m.state.S;
    pause(pause_time);
end

title('T');
hold off;

t = t/86400;

figure;
%suptitle("Van Roekel et al. 2018 Figure 6(c)");

subplot(2,1,1);
hold on;
plot(t, -F_sw_sim, 'r-');
plot(t, Hf, 'b-');
hold off;
title('Shortwave radiation F_{sw} (red), surface cooling (blue)');
ylabel('[W / m^2]');

subplot(2,1,2);
plot(t, h_sim, 'k-');
title('Mixed-layer Depth');
ylim([0 50]);
ylabel('Depth [m]');
xlabel('Time [days]');
set(gca, 'Ydir', 'reverse');

hold off;


figure;

%titles = {'T', 'S', 'K_s', "w'T'", "w'S'"};
%ylabels = {'degC', 'PSU', 'm^2/s', 'm/s degC', 'm/s PSU'};
%vars = {T, S, K_s, wT, wS};
%grids = {'T', 'T', 'W', 'W', 'W'};

titles = {'T', 'K_s', "w'T'"};
ylabels = {'degC', 'm^2/s', 'm/s degC'};
vars = {T, K_s, wT};
grids = {'T', 'W', 'W'};


ax = {};
for i=1:length(titles)
    ax{i} = subplot(length(titles), 1, i);
    hold(ax{i}, 'on');

    if (grids{i} == 'W')
        grid = m.grid.z_W;
    elseif (grids{i} == 'T')
        grid = m.grid.z_T;
    end
    
    contourf(t, grid, vars{i}, 'LineColor','none');
    ylabel('z [m]');
    h = colorbar();
    ylabel(h, ylabels{i});
    plot(t, -h_sim, 'LineWidth', 1, 'Color', 'white', 'LineStyle', '--');
    title(ax{i}, titles(i));
    
    ylim([-40 0]);
    
    hold(ax{i}, 'off');
end

xlabel(ax{end}, 'Time [ day ]');

