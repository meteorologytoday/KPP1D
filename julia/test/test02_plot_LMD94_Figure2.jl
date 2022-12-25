include(joinpath(@__DIR__, "..", "src", "KPP.jl"))
using .KPP
using Formatting

println(" *** Reproduce LMD94 Figure 2. Scalar = dashed. Momentum = solid. *** ")

H = 10
z_W = - collect(Float64, range(0, 1, length=1000)) * H

τ0 = 0.0
h = 5.0
d = - z_W

σ = d / h

println("Loading PyPlot... ")
using PyPlot
plt=PyPlot
println("Done.")


fig, ax = plt.subplots(1, 3, sharey=true)

fig.suptitle("Reproduce LMD94 Figure 2. Scalar = dashed. Momentum = solid.")

ax[1].set_title("\$G \\left(\\sigma\\right)\$")
ax[2].set_title("\$w_s / \\left(\\kappa u_*\\right)\$")
ax[3].set_title("\$w_m / \\left(\\kappa u_*\\right)\$")

ax[1].plot(KPP.calG(σ), σ, "k-")
h_over_L_arr = [1, 0.1, 0.0, -1.0, -5.0]
legend_str = []

for (i, h_over_L) in enumerate(h_over_L_arr)
    
    L_star_target = h / h_over_L

    u_star, _ = KPP.calMOSTscales(τ0, 1.0)
    wb_sfc = u_star^3 / ( - KPP.κ * L_star_target )
    w_s, _, u_star, L_star = KPP.calw_x(:SCALAR,   z_W, h, τ0, wb_sfc)
    w_x, _, u_star, L_star = KPP.calw_x(:MOMENTUM, z_W, h, τ0, wb_sfc)

    ax[2].plot( w_s ./ (KPP.κ * u_star), σ, "k--")
    ax[3].plot( w_x ./ (KPP.κ * u_star), σ, "k-")
    
    push!(legend_str, format("\$ h / L_* = {:.1f}\$", h_over_L))
    println(format("h_over_L[{:d}] = {:f}", i, h / L_star))
end

ax[1].set_xlim([0, 0.2])
ax[2].set_xlim([0, 4])    
ax[3].set_xlim([0, 4])    

ax[3].legend()
#ax{3}, legend_str, 'Location', 'south');

ax[1].set_ylim([1, 0])
#ax[1].invert_yaxis()


plt.show()
