clear;

varnames = { "T", "S", "b", "Ri", "b_{sl} - b", "u_{sl} - u", "V_t^2", "N"};
wb_sfc = 0.001;
Kv_iso = 1e-3;
dt = 0.1;
total_steps = 10 / dt;
Nz = 200;
H  = 25; % m

m = Model(H, Nz, Kv_iso, dt);
m.showModelInfo();

m.state.wb_sfc = wb_sfc;

slope_T = 5 / H;
slope_S = 1 / H;
T_sfc = 30;
S_sfc = 35;

m.state.T = T_sfc + slope_T * m.grid.z_T;
m.state.S = S_sfc + slope_S * m.grid.z_T;

dT = 2; dT_width = 5; dT_cent = -H/2;
dS = 0.5; dS_width = 3; dS_cent = -H/4;

gau = @(mag, cent, width) mag * exp(-( (m.grid.z_T - cent) / width ).^2 );

%m.state.T = m.state.T + gau(dT, dT_cent, dT_width);
%m.state.S = m.state.S + gau(dS, dS_cent, dS_width);
m.update_b();

figure;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [10 20]);

ax = {};

n=length(varnames);
for i=1:n
    ax{i} = subplot(1, n, i);
end


for i=1:n
    title(ax{i}, varnames{i});
end

ylabel(ax{1}, 'z (m)');

t(1) = 0;
int_S(1) = sum( m.grid.dz_T .* m.state.S);
int_T(1) = sum( m.grid.dz_T .* m.state.T);
int_b(1) = sum( m.grid.dz_T .* m.state.b);

for step = 1:total_steps

    
    t(end+1) = t(end) + dt;
    int_S(end+1) = sum( m.grid.dz_T .* m.state.S);
    int_T(end+1) = sum( m.grid.dz_T .* m.state.T);
    int_b(end+1) = sum( m.grid.dz_T .* m.state.b);
        
    if mod(step, 1) == 0
        [ h, Ri, db, du_sqr, Vt_sqr ] = m.update_ML();
        
        for i=1:n
            hold(ax{i}, 'on');
        end
        
        plot(ax{1}, m.state.T, m.grid.z_T);
        plot(ax{2}, m.state.S, m.grid.z_T);
        plot(ax{3}, m.state.b, m.grid.z_T);
        plot(ax{4}, Ri, m.grid.z_T);
        plot(ax{5}, db, m.grid.z_T);
        plot(ax{6}, du_sqr, m.grid.z_T);
        plot(ax{7}, Vt_sqr, m.grid.z_T);
        plot(ax{8}, m.grid.sop.T_interp_W * m.grid.sop.W_ddz_T * m.state.b, m.grid.z_T);
        
        for i=1:n
            hold(ax{i}, 'off');
        end
        
        
        fprintf('Mixed-layer k : %d\n', m.state.h_k)
        fprintf('Mixed-layer thickness: %f\n', m.state.h)
        
        return;
        pause(1);
        

    end
    
    m.stepModel();
    
    
end


figure;
subplot(3,1,1);
plot(t, int_T);
title("Integrated T");

subplot(3,1,2);
plot(t, int_S);
title("Integrated S");

subplot(3,1,3);
plot(t, int_b);
title("Integrated b");
xlabel('time (s)');
