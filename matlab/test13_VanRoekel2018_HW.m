clear;

% This simulation is meant to reproduce the HW case in
% Van Roekel et al (2018) where HW stands for "Heating and wind"

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 86400*3;
dt = 10*60;
f = 1e-4;
pause_time = 0.0;

total_steps = total_time / dt;

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

% Profile A
T_sfc = 20; dTdz = 0.01;
S_sfc = 35; dSdz = 0;
m.state.T = T_sfc + dTdz * m.grid.z_T;
m.state.S = S_sfc + dSdz * m.grid.z_T;
m.update_b();

% downward turbulent heat flux
m.state.Hf_sen = -75;

% wind stress
m.state.taux0 = 0.1;
m.state.tauy0 = 0;

figure;
varnames = {'T', 'u', 'v'};
n = length(varnames);
ax = {};
t(1) = 0;
h_sim(1) = m.state.h;
for i=1:n
    ax{i} = subplot(1,n,i);
end

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    m.stepModel(m.SURFFLUX_SIMPLE);
    
    plot(ax{1}, m.state.T, m.grid.z_T);
    plot(ax{2}, m.state.u, m.grid.z_T);
    plot(ax{3}, m.state.v, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;

    h_sim(end+1) = m.state.h;
    
    pause(pause_time);
end
for i=1:n
    title(ax{i}, varnames{i});
end

hold off;

figure;
hold on;
plot(t / 86400, h_sim, 'k-');
xlabel('Time [hr]');
ylabel('Depth [m]');
set(gca, 'Ydir', 'reverse');
title("Van Roekel et al. 2018 Figure 4(b)");
hold off;