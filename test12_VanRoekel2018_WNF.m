clear;

% This simulation is meant to reproduce the WNF case in
% Van Roekel et al (2018) where WNF stands for "Wind stress with no Coriolis"

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 86400;
dt = 10*60;
f = 0;
pause_time = 0.01;

total_steps = total_time / dt;

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

% Profile D
dTdz = 0.05;
m.state.T = 20 + dTdz * m.grid.z_T;
m.state.S(:) = 35;
m.update_b();

N = sqrt(m.grid.sop.W_ddz_T * m.state.b);
N = mean( N(2:end-1) );


m.state.taux0 = 0.1;
m.state.tauy0 = 0;

figure;
varnames = {'u', 'v'};
n = length(varnames);
ax = {};
t(1) = 0;
u_0_sim(1) = m.state.u(1);
v_0_sim(1) = m.state.v(1);
h_sim(1) = m.state.h;
for i=1:n
    ax{i} = subplot(1,n,i);
end

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    m.stepModel();
    
    plot(ax{1}, m.state.u, m.grid.z_T);
    plot(ax{2}, m.state.v, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;
    u_0_sim(end+1) = m.state.u(1);
    v_0_sim(end+1) = m.state.v(1);
    h_sim(end+1) = m.state.h;
    
    pause(pause_time);
end
for i=1:n
    title(ax{i}, varnames{i});
end

hold off;

% Van Roekel et al (2018) equation (27)
[ u_star, ~ ] = m.kpp.calMOSTscales(m.state.tau0, 0);
h_KP = 1.05 * u_star * sqrt(t / N);

figure;
hold on;
plot(t / 3600, h_KP, 'r-');
plot(t / 3600, h_sim, 'k-');
legend(gca, {'h_{KP}', 'h_{sim}'}, 'Location', 'south');
xlabel('Time [hr]');
ylabel('Depth [m]');
set(gca, 'Ydir', 'reverse');
title("Van Roekel et al. 2018 Figure 5(a)");
hold off;