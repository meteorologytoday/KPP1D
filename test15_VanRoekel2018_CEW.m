clear;

% This simulation is meant to reproduce the CEW case in
% Roekel et al (2018) where CEW stands for
% "Surface Cooling, Evaporation, and Wind stress"

Kv_iso = 1e-3;
Nz = 150;
H  = 150; % m
total_time = 86400 * 8;
dt = 10 * 60;
f = 1e-4;
Hf_sen = 75;
evap = 1.37 * (1e-3 / 86400.0); % mm / day
tau0 = 0.1;
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

m.state.Hf_sen = Hf_sen;
m.state.evap = evap;
m.state.taux0 = tau0;

t(1) = 0;
h(1) = m.state.h;

figure;
varnames = {'T', 'S', 'u', 'v'};
N = length(varnames);
ax = {};

for i=1:N
    ax{i} = subplot(1,N,i);
end

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    [ diag_kpp ] = m.stepModel();
    t(end+1) = t(end) + m.dt;
    h(end+1) = m.state.h;
    
    plot(ax{1}, m.state.T, m.grid.z_T);
    plot(ax{2}, m.state.S, m.grid.z_T);
    plot(ax{3}, m.state.u, m.grid.z_T);
    plot(ax{4}, m.state.v, m.grid.z_T);
    
    
    pause(pause_time);
end
for i=1:N
    title(ax{i}, varnames{i});
end

hold off;

figure;
plot( t / 86400, h, 'k-' );
title('Roekel et al (2018) Figure 4(a)');
xlabel('Time [days]');
ylabel('h [m]');
set(gca, 'Ydir', 'reverse');