using NCDatasets
using ArgParse


function parse_commandline()

    s = ArgParseSettings()

    @add_arg_table! s begin

        "--z_W"
            help = "z_W values"
            arg_type = Float64
            nargs = '+'
            default = [0.0, -10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0, -90.0, -100.0]

        "--output"
            help = "Name of output file"
            arg_type = String
            default = "zgrid.nc"
    end

    args = parse_args(s)

    println("Received Parameters:")
    for (arg, val) in args
        println("  $arg  =>  $val")
    end

    return args

end

args = parse_commandline()

z_W = args["z_W"]
Nz = length(z_W) - 1

println("Output file: ", args["output"])
Dataset(args["output"], "c") do ds

    defDim(ds, "Nzp1", Nz+1)

    for (varname, vardata, vardim, attrib) in [
        ("z_W",  z_W, ("Nzp1", ), Dict()),
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
