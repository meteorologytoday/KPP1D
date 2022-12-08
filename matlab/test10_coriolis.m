clear;

% This simulation only test the "calRK4Coriolis" method in Model.m

Kv_iso = 1e-3;
Nz = 50;
H  = 50; % m
total_time = 86400;
dt = 30 * 60;
f = 1e-4;
pause_time = 0.0;

total_steps = total_time / dt;

m = Model(H, Nz, Kv_iso, dt, f);
m.showModelInfo();

u = 1;
v = 0;

t(1) = 0;
u_sim(1) = u;
v_sim(1) = v;

for step = 1:total_steps
    
    fprintf('Step %d\n', step);
    [ u, v ] = m.calRK4Coriolis(m.f, m.dt, u, v);
    t(end+1) = t(end) + m.dt;
    u_sim(end+1) = u;
    v_sim(end+1) = v;
   
end

% analytical solution
u_true =   u_sim(1) * cos( m.f * t ) + v_sim(1) * sin( m.f * t );
v_true = - u_sim(1) * sin( m.f * t ) + v_sim(1) * cos( m.f * t );

hold on;
plot( t / 86400, u_sim, 'r-' );
plot( t / 86400, v_sim, 'r--');
plot( t / 86400, u_true, 'k-' );
plot( t / 86400, v_true, 'k--');
title('u, v vs time');
legend(gca, {'u_{sim}', 'v_{sim}', 'u_{true}', 'v_{true}'}, 'Location', 'south');
xlabel('Time [days]');
ylabel('velocity [m / s]');