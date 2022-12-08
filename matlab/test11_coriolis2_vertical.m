clear;

% This simulation is meant to reproduce the FC case in
% Roekel et al (2018) where FC stands for "free convection"

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 3600 * 6;
dt = 10*60;
f = 1e-4;
wT_0 = 75 / 1024 / 3996;
pause_time = 0.1;

total_steps = total_time / dt;


m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();


shear_bnd = [-5 -15];
x = (m.grid.z_T - shear_bnd(1)) / (shear_bnd(2) - shear_bnd(1));
U = 1 - 3 * x.^2 + 2 * x.^3;
U(m.grid.z_T > shear_bnd(1)) = 1;
U(m.grid.z_T <= shear_bnd(2)) = 0;

angle = 30 * pi / 180;
m.state.u = U * cos(angle);
m.state.v = U * sin(angle);

figure;
varnames = {'u', 'v'};
N = length(varnames);
ax = {};
t(1) = 0;
u_0_sim(1) = m.state.u(1);
v_0_sim(1) = m.state.v(1);
for i=1:N
    ax{i} = subplot(1,N,i);
end

hold on;

for step = 1:total_steps
    fprintf('Step %d\n', step);
    m.stepModel_momentum();
    
    plot(ax{1}, m.state.u, m.grid.z_T);
    plot(ax{2}, m.state.v, m.grid.z_T);
    
    t(end+1) = t(end) + m.dt;
    u_0_sim(end+1) = m.state.u(1);
    v_0_sim(end+1) = m.state.v(1);
    
    pause(pause_time);
end
for i=1:N
    title(ax{i}, varnames{i});
end

hold off;

figure;
hold on;
plot(t / 3600, u_0_sim, 'k-');
plot(t / 3600, v_0_sim, 'k--');
legend(gca, {'u', 'v'}, 'Location', 'south');
xlabel('Time [hr]');
ylabel('Velocity [m/s]');
hold off;