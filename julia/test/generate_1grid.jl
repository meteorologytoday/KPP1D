using NCDatasets
using ArgParse


function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table! s begin

        "--cent-lat"
            help = "an option with an argument"
            arg_type = Float64
            default = 0.0

        "--dlat"
            help = "The size of box in lat"
            arg_type = Float64
            default = 1.0

        "--nlon"
            help = "Number of boxes in longitude. Minimum = 2"
            arg_type = Int64
            default = 2

        "--output"
            help = "Name of output file"
            arg_type = String
            default = "output.nc"
    end

    args = parse_args(s)

    println("Received Parameters:")
    for (arg, val) in args
        println("  $arg  =>  $val")
    end

    return args

end

args = parse_commandline()

Nx = args["nlon"]
Ny = 1
Nv = 4

xc = zeros(Float64, Nx, Ny)
yc = zeros(Float64, Nx, Ny)

xv = zeros(Float64, Nv, Nx, Ny)
yv = zeros(Float64, Nv, Nx, Ny)

area = zeros(Float64, Nx, Ny)
frac = zeros(Float64, Nx, Ny)
mask = zeros(Int64, Nx, Ny)


dlon = 360.0 / Nx
dlat = args["dlat"]
cent_lat = args["cent-lat"]

for i = 1:Nx
    
    xc[i, 1] = dlon * (i - 0.5)
    yc[i, 1] = cent_lat
 
    xv[1, i, 1] = xc[i, 1] - dlon / 2.0
    xv[2, i, 1] = xc[i, 1] + dlon / 2.0
    xv[3, i, 1] = xv[2, i, 1]
    xv[4, i, 1] = xv[1, i, 1]

    yv[1, i, 1] = yc[i, 1] - dlat / 2.0
    yv[2, i, 1] = yv[1, i, 1]
    yv[3, i, 1] = yc[i, 1] + dlat / 2.0
    yv[4, i, 1] = yv[3, i, 1]

    mask[i, 1] = 1

    area[i, 1] = 2π * ( sin(deg2rad(yv[3, i, 1])) - sin(deg2rad(yv[2, i, 1])) )
    frac[i, 1] = 1.0

end

println("Output file: ", args["output"])
Dataset(args["output"], "c") do ds

    defDim(ds, "ni", Nx)
    defDim(ds, "nj", Ny)
    defDim(ds, "nv",  4)

    for (varname, vardata, vardim, attrib) in [
        ("xc",  xc, ("ni", "nj",), Dict()),
        ("yc",  yc, ("ni", "nj",), Dict()),
        ("xv",  xv, ("nv", "ni", "nj",), Dict()),
        ("yv",  yv, ("nv", "ni", "nj",), Dict()),
        ("mask", mask, ("ni", "nj",), Dict()),
        ("area", area, ("ni", "nj",), Dict()),
        ("frac", frac, ("ni", "nj",), Dict()),
    ]

        println("Doing var: ", varname)

        var = defVar(ds, varname, eltype(vardata), vardim)

        if eltype(vardata) <: Union{Float64, Float32}
            var.attrib["_FillValue"] = 1e20
        end

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end
