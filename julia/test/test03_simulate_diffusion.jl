include(joinpath(@__DIR__, "..", "src", "KPP.jl"))
using .KPP
using Formatting

println(" *** This is a simulation to test convective adjustment. *** ")

H  = 25.0   # m
Nz = 1000
f = 1e-4
Kv_iso = 1e-3
dt = 0.1
total_time = 20.0
total_steps = total_time / dt

z_W = - collect(Float64, range(0, 1, length=1000)) * H


ev = KPP.Env(KPP.Grid(z_W = z_W))
m = KPP.Model(ev)

slope_T = 2.0 / H
slope_S = -1.0 / H
T_sfc = 30
S_sfc = 35
dT = 2.0; dT_width = 5.0; dT_cent = -H/2;
dS = 0.5; dS_width = 3.0; dS_cent = -H/4;


gaussian(mag, cent, width) = mag .* exp.(-( (m.ev.gd.z_T .- cent) ./ width ).^2 )

@. m.st.T = T_sfc + slope_T * m.ev.gd.z_T
@. m.st.S = S_sfc + slope_S * m.ev.gd.z_T

m.st.T .+= gaussian(dT, dT_cent, dT_width)
m.st.S .+= gaussian(dS, dS_cent, dS_width)

KPP.updateBuoyancy!(m)

#=
println("Loading PyPlot... ")
using PyPlot
plt=PyPlot
println("Done.")

figure;

ax1 = subplot(1, 3, 1);
ax2 = subplot(1, 3, 2);
ax3 = subplot(1, 3, 3);

title(ax1, "T");
title(ax2, "S");
title(ax3, "b");
ylabel(ax1, 'z (m)');

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
        
        hold(ax1, 'on');
        plot(ax1, m.state.T, m.grid.z_T);
        hold(ax1, 'off');
        
        hold(ax2, 'on');
        plot(ax2, m.state.S, m.grid.z_T);
        hold(ax2, 'off');
        
        hold(ax3, 'on');
        plot(ax3, m.state.b, m.grid.z_T);
        hold(ax3, 'off');
        
        pause(.01);
        

    end
    
    m.stepModel(m.SURFFLUX_SIMPLE);
    
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
=#
