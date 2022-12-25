include(joinpath(@__DIR__, "..", "src", "KPP.jl"))
using .KPP
using Formatting

println(" *** This subplot should look like LMD94 Figure B1. scalar *** ")

H = 10
z_W = - collect(Float64, range(0, 1, length=1000)) * H

h = 5.0
d = - z_W
τ0 = 0.1

u_star, L_star = KPP.calMOSTscales(τ0, 0.1 * 0.0001)

println("u_star = ", u_star)
println("L_star = ", L_star)

ϕ_s_conv = KPP.calϕ_s(d,   L_star)
ϕ_s_stab = KPP.calϕ_s(d, - L_star)

ϕ_m_conv = KPP.calϕ_m(d,   L_star)
ϕ_m_stab = KPP.calϕ_m(d, - L_star)


println("Loading PyPlot... ")
using PyPlot
plt=PyPlot
println("Done.")


fig, ax = plt.subplots(1, 1)

# A scale to magnified positive x-axis 
# so that the figures are comparable
scale = 10.0

ax.set_title(format("\$u_* = {:.2f}\$, \$L_* = {:.2f}\$", u_star, L_star))



ax.plot(  d / L_star, ϕ_s_conv, "k--")
ax.plot(- d / L_star * scale, ϕ_s_stab, "k--")

ax.plot(  d / L_star, ϕ_m_conv, "k-")
ax.plot(- d / L_star * scale, ϕ_m_stab, "k-")

ax.text(-1,  1, "Unstable", fontsize=15, va="center", ha="center")
ax.text( 1, 1, "Stable", fontsize=15, va="center", ha="center")

ax.text( -1.5, 0.50, "\$\\phi_m\$", fontsize=15, va="center", ha="center")
ax.text( -1.0, 0.10, "\$\\phi_s\$", fontsize=15, va="center", ha="center")
ax.text(  1.0, 1.25, "\$\\phi_m = \\phi_s\$", fontsize=15, va="center", ha="center")

ax.axvline(x=0, color="r", linestyle="dashed")

ax.set_xlim([-2, 2])
ax.set_ylim([0, 2])

ax.set_xticks([-2, -1, 0, 1, 2])
ax.set_xticklabels(["-2", "-1", "0", ".1", ".2"])

plt.show()
