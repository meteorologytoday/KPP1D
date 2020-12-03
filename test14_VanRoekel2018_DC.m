clear;

% This simulation is meant to reproduce the WNF case in
% Van Roekel et al (2018) where WNF stands for "Wind stress with no Coriolis"

Kv_iso = 1e-3;
Nz = 150;
H  = 150; % m
total_time = 86400 * 8;
dt = 10*60;
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

figure;

t(1) = 0;
h_sim(1) = m.state.h;
F_sw_sim(1) = F_sw_func(t(1)); 

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    
    
    m.state.I_0 = F_sw_sim(end);    
    m.stepModel();
    
    plot(m.state.T, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;
    h_sim(end+1) = m.state.h;
    F_sw_sim(end+1) = F_sw_func(t(end));    
    
    pause(pause_time);
end

title('T');
hold off;


figure;
suptitle("Van Roekel et al. 2018 Figure 6(c)");

subplot(2,1,1);
plot(t / 86400, F_sw_sim, 'k-');
title('F_{sw}');
ylabel('Depth [m]');

subplot(2,1,2);
plot(t / 86400, h_sim, 'k-');
title('h');
ylim([0 50]);
ylabel('Depth [m]');
xlabel('Time [days]');
set(gca, 'Ydir', 'reverse');

hold off;