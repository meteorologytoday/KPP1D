include(joinpath(@__DIR__, "..", "src", "KPP.jl"))


using .KPP


gd = KPP.loadGrid("zgrid.nc")
ev = KPP.Env(gd)
model = KPP.Model(
    ev,
)





