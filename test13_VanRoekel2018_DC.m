clear;

% This simulation is meant to reproduce the WNF case in
% Van Roekel et al (2018) where WNF stands for "Wind stress with no Coriolis"

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 86400;
dt = 10*60;
f = 1e-4;
Q_0 = 75;
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
F_sw = @(Q, t) - max(0, cos(t / 86400 - 0.5));

m.state.taux0 = 0.1;
m.state.tauy0 = 0;

figure;

t(1) = 0;
h_sim(1) = m.state.h;

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    m.stepModel();
    
    plot(m.state.T, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;
    h_sim(end+1) = m.state.h;
    
    pause(pause_time);
end

title('T');
hold off;


figure;
subplot(2,1,1);
plot(t / 86400, F_sw, 'k-');
plot(t / 3600, h_sim, 'k-');
legend(gca, {'h_{KP}', 'h_{sim}'}, 'Location', 'south');
xlabel('Time [hr]');
ylabel('Depth [m]');
set(gca, 'Ydir', 'reverse');
title("Van Roekel et al. 2018 Figure 5(a)");
hold off;