clear;

omega = 2 * pi / (86400 * (360 / 365) );
lat = 20.7024;
f = 2 * omega * sin( pi * lat / 180);
beg_sim_date = datenum(2019, 5, 28, 12,0,0);
end_sim_date = datenum(2019, 5, 29, 12,0,0);

dt = 60 * 10; % 2 min

momentum_nudging_timescale = 30*60;

%%%% Initialization %%%%
disp( datestr(beg_sim_date, 'YYYY-mm-dd HH:MM:SS') );
disp( datestr(end_sim_date, 'YYYY-mm-dd HH:MM:SS') );

total_time = (end_sim_date - beg_sim_date) * 86400;
total_steps = total_time / dt;
t_sim = ( linspace(0, total_time, total_steps+1) ) / 86400 + beg_sim_date;
t_sim_mid = ( t_sim(2:end) + t_sim(1:end-1) ) / 2;

z_W = - linspace(0, 15, 101);
m = Model(z_W, dt, f);
m.momentum_nudging_timescale = momentum_nudging_timescale;
m.showModelInfo();

%%%%% forcing data %%%%%
sf = SurfaceFlux();
forcing.t      = t_sim_mid;
forcing.z_T    = m.grid.z_T;
forcing.solar  = zeros(length(t_sim_mid), 1);
forcing.wind_u = 0 * forcing.solar + 1e-5;
forcing.wind_v = 0 * forcing.solar;
forcing.T_a    = 0 * forcing.solar;
forcing.q_a    = 0 * forcing.solar;
forcing.precip = 0 * forcing.solar;
forcing.Ux     = zeros(m.grid.Nz, length(t_sim_mid));
forcing.Uy     = zeros(m.grid.Nz, length(t_sim_mid));
forcing.T      = zeros(m.grid.Nz, length(t_sim_mid));
forcing.S      = zeros(m.grid.Nz, length(t_sim_mid));

Q_sw_max = 340*4; %235.62;

% Assign forcing here
for i=1:length(t_sim_mid)
    t_now = t_sim_mid(i);
    forcing.solar(i) = - Q_sw_max * max(0, cos(2 * pi * (t_now - 0.5)));
    forcing.T_a(i) = 10;
    forcing.q_a(i) = 0 / 1000;
    forcing.wind_u(i) = 1.0;
    
    slope_T = 0.01;
    slope_S = 0;
    T_0 = 20;
    S_0 = 35;
    forcing.T(:, i) = T_0 + slope_T * m.grid.z_T;
    forcing.S(:, i) = S_0 + slope_S * m.grid.z_T;

end

%%%%% assign initial condition %%%%%
m.state.I_0 = forcing.solar(1);
m.state.U10 = forcing.wind_u(1);
m.state.V10 = forcing.wind_v(1);
m.state.T_a = forcing.T_a(1);
m.state.q_a = forcing.q_a(1);
m.state.precip = forcing.precip(1);
m.state.u_nudging = forcing.Ux(:, 1);
m.state.v_nudging = forcing.Uy(:, 1);

m.state.T(:) = forcing.T(:, 1);
m.state.S(:) = forcing.S(:, 1);
m.state.u(:) = forcing.Ux(:, 1);
m.state.v(:) = forcing.Uy(:, 1);
m.state.u(isnan(m.state.u)) = 0;
m.state.v(isnan(m.state.v)) = 0;
m.update_b();

figure;

time_ticks = (t_sim(1)):1:(t_sim(end));
time_labels = {};
empty_labels = {};
for i=1:length(time_ticks)
    time_labels{i} = datestr(time_ticks(i), "mm/dd hh:MM");
    empty_labels{i} = '';
end

met_varnames = {'solar', 'wind_u', 'wind_v', 'T_a', 'q_a', 'precip'};
for i=1:length(met_varnames)
    varname = met_varnames{i};
    ax(i) = subplot(length(met_varnames), 1, i);
    plot( ax(i), forcing.t, forcing.(varname) );
    
    title(ax(i), varname);
    xticks(ax(i), time_ticks);
    xticklabels(ax(i), empty_labels);
end
xticklabels(ax(end), time_labels);

figure;
ax = [];
for i=1:5
    ax(i) = subplot(1, 5, i);
end

plot(ax(1), m.state.T, m.grid.z_T );
plot(ax(2), m.state.S, m.grid.z_T );
plot(ax(3), m.state.b, m.grid.z_T );
plot(ax(4), m.state.u, m.grid.z_T );
plot(ax(5), m.state.v, m.grid.z_T );
    title(ax(1), 'T');
    title(ax(2), 'S');
    title(ax(3), 'b');
    title(ax(4), 'u');
    title(ax(5), 'v');
%xlim(ax(1), [27, 30]);
%xlim(ax(2), [33.5, 34]);
xlim(ax(4), [-0.5, 0.5]);
xlim(ax(5), [-0.5, 0.5]);
pause();

%%%%% Ready to run model %%%%%
rec = Recorder(m.grid.Nz, total_steps+1);
rec.record(1, m);

for step=1:total_steps
    
    % assign forcing
    % U10, V10, T_a, q_a, T_o, I_0
    m.state.I_0 = forcing.solar(step);
    m.state.U10 = forcing.wind_u(step);
    m.state.V10 = forcing.wind_v(step);
    m.state.T_a = forcing.T_a(step);
    m.state.q_a = forcing.q_a(step);
    m.state.precip = forcing.precip(step);
    m.state.u_nudging(:) = nan;% = forcing.Ux(:, step);
    m.state.v_nudging(:) = nan;% = forcing.Uy(:, step);
    
    m.stepModel(m.SURFFLUX_REALISTIC);
    rec.record(step+1, m);
    
    plot(ax(1), m.state.T, m.grid.z_T );
    plot(ax(2), m.state.S, m.grid.z_T );
    plot(ax(3), m.state.b, m.grid.z_T );
    plot(ax(4), m.state.u, m.grid.z_T );
    plot(ax(5), m.state.v, m.grid.z_T );
    
    %xlim(ax(1), [27, 35]);
    %xlim(ax(2), [33.5, 34]);
    %xlim(ax(4), [-0.5, 0.5]);
    %xlim(ax(5), [-0.5, 0.5]);
    
    title(ax(1), 'T');
    title(ax(2), 'S');
    title(ax(3), 'b');
    title(ax(4), 'u');
    title(ax(5), 'v');
    
    suptitle(sprintf('%s (h=%.1f, I=%.1f, prec=%.1f, evap=%.1f, Hf_sen=%.1f, Hf_lat=%.1f, Hf_lw=%.1f)', datestr(t_sim(step), "mm/dd hh:MM"), m.state.h, m.state.I_0, m.state.precip*1000*3600, m.state.evap*1000*3600, m.state.Hf_sen, m.state.Hf_lat, m.state.Hf_lw));
    
    %pause(1);
end

%save('result.mat', '-struct', 'rec');
figure;
ax = [];
N = 5;
for i=1:N
    ax(i) = subplot(N, 1, i);
end

hold(ax(1), 'on');
plot(ax(1), t_sim, rec.T_a, 'k-');
plot(ax(1), t_sim, rec.T(1, :), 'r-');
title(ax(1), 'T_a (black), SST (red)');
hold(ax(1), 'off');


hold(ax(2), 'on');
plot(ax(2), t_sim, rec.Hf_sen, 'k-');
plot(ax(2), t_sim, rec.Hf_lat, 'b-');
plot(ax(2), t_sim, rec.Hf_lw,  'r-');
plot(ax(2), t_sim, rec.I_0,  'r--');
title(ax(2), 'Hf\_sen (black), Hf\_lat (blue), Hf\_lw (red), I\_0 (red-dashed)');
hold(ax(2), 'off');

hold(ax(3), 'on');
plot(ax(3), t_sim, rec.h, 'k-');
title(ax(3), 'Mixed-layer depth');
hold(ax(3), 'off');

hold(ax(4), 'on');
plot(ax(4), t_sim, rec.precip, 'b-');
plot(ax(4), t_sim, rec.evap, 'r-');
title(ax(4), 'Precipitation (blue) and evaporation (red)');
hold(ax(4), 'off');

figure;
hold on;
contourf(t_sim, m.grid.z_T, rec.T);
plot(t_sim, - rec.h, 'w-');
hold off;