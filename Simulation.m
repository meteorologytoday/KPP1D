clear;

omega = 2 * pi / (86400 * (360 / 365) );
lat = 20.7024;
f = 2 * omega * sin( pi * lat / 180);
beg_sim_date = datenum(2019, 5,28,12,0,0);
end_sim_date = datenum(2019, 6, 1,12,0,0);
%end_sim_date = datenum(2019,6, 5, 0,0,0);
datafile = 'data/processed-WWE_binned_0p05m.mat';
dt = 120; % 2 min
momentum_nudging_timescale = 30*60;

%%%% Initialization %%%%
data = load(datafile);

disp( datestr(data.t(1), 'YYYY-mm-dd HH:MM:SS') );
disp( datestr(data.t(end), 'YYYY-mm-dd HH:MM:SS') );

if (data.t(1) > beg_sim_date || data.t(end) < end_sim_date)
    error('Simulation time out of possible range.');
end

total_time = (end_sim_date - beg_sim_date) * 86400;
total_steps = total_time / dt;
t_sim = ( linspace(0, total_time, total_steps+1) ) / 86400 + beg_sim_date;
t_sim_mid = ( t_sim(2:end) + t_sim(1:end-1) ) / 2;

z_T = data.z;
z_W = zeros(length(z_T)+1, 1);
z_W(2:end-1) = (z_T(1:end-1) + z_T(2:end) ) / 2;
z_W(1) = 0;
z_W(end) = z_W(end-1) - ( z_W(end-2) - z_W(end-1) );

m = Model(z_W, dt, f);
m.momentum_nudging_timescale = momentum_nudging_timescale;
m.showModelInfo();

%%%%% interpolate meteorological data %%%%%
sf = SurfaceFlux();

% Unit conversion
data.solar  = - data.solar * 1000; % kW / m^2 => W / m^2 and positive upward
data.precip = data.precip * 0.001 / 3600; % mm / hr => m / s
data.RH     = data.RH / 100; %  percentage => none

% Calculate air humidity
data.q_a = sf.calSaturatedSpecificHumidity(data.T_a) .* data.RH;

met_varnames = {'solar', 'wind_u', 'wind_v', 'T_a', 'RH', 'q_a', 'precip'};

forcing.t = t_sim_mid;
forcing.z_T = m.grid.z_T;
for i = 1:length(met_varnames)
    varname = met_varnames{i};
    tmp = griddedInterpolant(data.t, data.(varname));
    forcing.(varname) = tmp(t_sim_mid);
end

%%%%% interpolate ocean profile %%%%%
ocn_varnames = { 'Ux', 'Uy', 'T', 'S' };
[ tt, zz ] = meshgrid(data.t, data.z);
[ tt_interp, zz_interp ] = meshgrid(t_sim(1:end-1), data.z);

for i = 1:length(ocn_varnames)
    varname = ocn_varnames{i};

    fprintf('Interpolating ocean variable: %s\n', varname);
    
    % structured interpolation
    d = data.(varname);
    d(d==-999) = nan;
    forcing.(varname) = interp2(tt, zz, d, tt_interp, zz_interp);
end


%%%%% assign initial condition %%%%%
m.state.I_0 = forcing.solar(1);
m.state.U10 = forcing.wind_u(1);
m.state.V10 = forcing.wind_v(1);
m.state.T_a = forcing.T_a(1);
m.state.q_a = forcing.q_a(1);
m.state.u_nudging = forcing.Ux(:, 1);
m.state.v_nudging = forcing.Uy(:, 1);

m.state.T(:) = forcing.T(:, 1);
m.state.S(:) = forcing.S(:, 1);
m.state.u(:) = forcing.Ux(:, 1);
m.state.v(:) = forcing.Uy(:, 1);
m.state.u(isnan(m.state.u)) = 0;
m.state.v(isnan(m.state.v)) = 0;

figure;

time_ticks = (t_sim(1)):1:(t_sim(end));
time_labels = {};
empty_labels = {};
for i=1:length(time_ticks)
    time_labels{i} = datestr(time_ticks(i), "mm/dd hh:MM");
    empty_labels{i} = '';
end


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
for i=1:4
    ax(i) = subplot(1, 4, i);
end

plot(ax(1), m.state.T, m.grid.z_T );
plot(ax(2), m.state.S, m.grid.z_T );
plot(ax(3), m.state.u, m.grid.z_T );
plot(ax(4), m.state.v, m.grid.z_T );
xlim(ax(1), [27, 30]);
xlim(ax(2), [33.5, 34]);
xlim(ax(3), [-0.5, 0.5]);
xlim(ax(4), [-0.5, 0.5]);
pause();

%%%%% Ready to run model %%%%%
rec = Recorder(m.grid.Nz, total_steps+1);
rec.record(1, m);

for step=1:total_steps
    
    % assign forcing
    % U10, V10, T_a, q_a, T_o, I_0
    m.state.I_0 = -2000;%forcing.solar(step);
    m.state.U10 = 1; %forcing.wind_u(step);
    m.state.V10 = 0; %forcing.wind_v(step);
    m.state.T_a = 30; %forcing.T_a(step);
    m.state.q_a = 0.02; % forcing.q_a(step);
    m.state.u_nudging(:) = nan;% = forcing.Ux(:, step);
    m.state.v_nudging(:) = nan;% = forcing.Uy(:, step);
    
    m.stepModel(m.SURFFLUX_SIMPLE);
    rec.record(step+1, m);
    
    plot(ax(1), m.state.T, m.grid.z_T );
    plot(ax(2), m.state.S, m.grid.z_T );
    plot(ax(3), m.state.u, m.grid.z_T );
    plot(ax(4), m.state.v, m.grid.z_T );
    
    xlim(ax(1), [27, 30]);
    xlim(ax(2), [33.5, 34]);
    xlim(ax(3), [-0.5, 0.5]);
    xlim(ax(4), [-0.5, 0.5]);
    
    pause(1);
end

save('result.mat', '-struct', 'rec');

