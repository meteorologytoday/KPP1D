clear;

varnames = { "T", "b", "N"};

Kv_iso = 1e-3;
Nz = 100;
H  = 25; % m
I0 = - 1000;
%I0 = 1000;
total_time = 86400;
dt = 10*60;
total_steps = total_time / dt;


m = Model(H, Nz, Kv_iso, dt);
m.showModelInfo();

slope_T = 5 / H;
slope_S = 1 / H;
T_sfc = 30;
S_sfc = 35;

m.state.T = T_sfc + slope_T * m.grid.z_T;
m.state.S = S_sfc + slope_S * m.grid.z_T;
m.state.I0 = I0;
m.update_b();

[ F, Q ] = m.rad.calRadiation(I0);

ax{1} = subplot(1,2,1);
plot(F, m.grid.z_W);
title("Irradiance Flux");

ax{2} = subplot(1,2,2);
plot(Q, m.grid.z_T);
title("Irradiance Flux convergence");

linkaxes([ax{1} ax{2}], 'y');


%%%% simulation %%%%

figure;

ax1 = subplot(1, 2, 1);
ax2 = subplot(1, 2, 2);

title(ax1, "T");
title(ax2, "b");
ylabel(ax1, 'z (m)');

t(1) = 0;
int_T(1) = sum( m.grid.dz_T .* m.state.T);
int_b(1) = sum( m.grid.dz_T .* m.state.b);

for step = 1:total_steps
        
    if mod(step, 1) == 0
        
        hold(ax1, 'on');
        plot(ax1, m.state.T, m.grid.z_T);
        hold(ax1, 'off');
        
        hold(ax2, 'on');
        plot(ax2, m.state.b, m.grid.z_T);
        hold(ax2, 'off');
        
        pause(.1);
        

    end
    
    %m.stepModel_radiation();
    m.stepModel();
    
    t(end+1) = t(end) + dt;
    int_T(end+1) = sum( m.grid.dz_T .* m.state.T);
    int_b(end+1) = sum( m.grid.dz_T .* m.state.b);
end


figure;
subplot(2,1,1);
plot(t, int_T);
title("Integrated T");

subplot(2,1,2);
plot(t, int_b);
title("Integrated b");
xlabel('time (s)');


%%%% cal total change of T %%%%
fprintf("Expected input energy = (-1) * total time * I0 = %f (J)\n", - total_time * m.state.I0);
fprintf("Simulated change = %f (J) \n", ( int_T(end) - int_T(1) ) * m.c.rho_sw * m.c.cp_sw);




