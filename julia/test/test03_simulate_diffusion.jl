include(joinpath(@__DIR__, "..", "src", "KPP.jl"))
using .KPP
using Formatting

println(" *** This is a simulation to test convective adjustment. *** ")

H  = 100.0   # m
Nz = 1000
f = 1e-4
Kv_iso = 1e-3
dt = 0.1    # second
total_time = 3600.0
total_steps = Int64(total_time / dt)

plot_interval_physical_time = 10 * 60.0
plot_interval_steps = floor(Int64, plot_interval_physical_time / dt)


z_W = - collect(Float64, range(0, 1, length=1000)) * H


ev = KPP.Env(
    gd = KPP.Grid(z_W = z_W),
    f  = f,
)
m = KPP.Model(ev)

m.fo.I_0 = -200.0

m.fo.surf_flux_calculation_type = :NOCALCULATION
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


# =====

t    = zeros(Float64, total_steps + 1)
∫Sdz = zeros(Float64, total_steps + 1)
∫Tdz = zeros(Float64, total_steps + 1)
∫bdz = zeros(Float64, total_steps + 1)

# =====

println("Loading PyPlot... ")
using PyPlot
plt=PyPlot
println("Done.")


fig, ax = plt.subplots(1, 3, sharey=true)

ax[1].set_title("T")
ax[2].set_title("s")
ax[3].set_title("b")

ax[1].set_ylabel("z [m]")

t[1] = 0.0

gd = m.ev.gd
∫Sdz[1] = sum( gd.Δz_T .* m.st.S)
∫Tdz[1] = sum( gd.Δz_T .* m.st.T)
∫bdz[1] = sum( gd.Δz_T .* m.st.b)

plt.ion()
plt.show()

function plot(step :: Int64)

    ax[1].plot(m.st.T, gd.z_T)
    ax[2].plot(m.st.S, gd.z_T)
    ax[3].plot(m.st.b, gd.z_T)

    fig.suptitle(format("t = {:.1f} s", t[step]))
        
    plt.draw()
end


plot(1)

for step = 1:total_steps

    println(format("t = {:.2f} s", t[step]))
    
   
    KPP.stepModel!(m, dt)

    t[step+1] = t[step] + dt
    ∫Sdz[step+1] = sum( gd.Δz_T .* m.st.S)
    ∫Tdz[step+1] = sum( gd.Δz_T .* m.st.T)
    ∫bdz[step+1] = sum( gd.Δz_T .* m.st.b)

    if mod(step, plot_interval_steps) == 0
        
        plot(step+1)
        sleep(0.1)

    end
 
end


#

fig, ax = plt.subplots(3, 1, sharex=true)

ax[1].plot(t, ∫Tdz)
ax[1].set_title("\$ \\int T \\, \\mathrm{d}z \$")

ax[2].plot(t, ∫Sdz)
ax[2].set_title("\$ \\int S \\, \\mathrm{d}z \$")

ax[3].plot(t, ∫bdz)
ax[3].set_title("\$ \\int b \\, \\mathrm{d}z \$")



ax[3].set_xlabel("time [s]")

plt.show(block=true)
