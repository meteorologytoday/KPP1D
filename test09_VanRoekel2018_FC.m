clear;

% This simulation is meant to reproduce the FC case in
% Roekel et al (2018) where FC stands for "free convection"

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 86400 * 3;
dt = 30 * 60;
f = 1e-4;
wT_0 = 75 / 1024 / 3996;
pause_time = 0.0;

total_steps = total_time / dt;


m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

slope_T = 0.01;
slope_S = 0;
T_0 = 20;
S_0 = 35;

m.state.T = T_0 + slope_T * m.grid.z_T;
m.state.S = S_0 + slope_S * m.grid.z_T;
m.update_b();


m.state.taux0 = 1e-2;
m.state.tauy0 = 0;
m.state.wT_0 = wT_0;

t(1) = 0;
h(1) = m.state.h;

figure;
varnames = {'T', 'S', 'Ri', 'K_s \gamma_T', 'K_s'};
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
    plot(ax{3}, m.state.Ri, m.grid.z_T);
    plot(ax{4}, diag_kpp.nloc_flux_T, m.grid.z_W);
    plot(ax{5}, diag_kpp.K_s_ML + diag_kpp.K_s_INT, m.grid.z_W);
    

    ylim(ax{3}, [-H 0]);    
    
    pause(pause_time);
end
for i=1:N
    title(ax{i}, varnames{i});
end

hold off;

figure;
plot( t / 86400, h, 'k-' );
title('h vs time');
xlabel('Time [days]');
ylabel('h [m]');
set(gca, 'Ydir', 'reverse');