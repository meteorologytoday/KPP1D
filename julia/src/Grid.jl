using NCDatasets
using LinearAlgebra: ⋅, normalize, normalize!, norm

mutable struct Grid

    Nz    :: Integer

    z_T :: AbstractArray{Float64, 1}
    z_W :: AbstractArray{Float64, 1}

    Δz_T  :: AbstractArray{Float64, 1}
    Δz_W  :: AbstractArray{Float64, 1}
  
    function Grid(;
        z_W    :: AbstractArray{Float64, 1}, 
    )

        if ! all(isfinite.(z_W))
            throw(ErrorException("At least one element of z_W is not finite."))
        end



        if z_W[1] != 0
            throw(ErrorException("First element of z_W must be 0"))
        end


        Nz = length(z_W) - 1

        z_T = ( z_W[1:end-1] + z_W[2:end] ) / 2.0

        Δz_T = z_W[1:end-1] - z_W[2:end]

        if any(Δz_T .< 0)
            throw(ErrorException("z_W should be monotically decreasing"))
        end



        Δz_W = zeros(Float64, Nz+1)
        Δz_W[2:Nz] = (Δz_T[1:end-1] + Δz_T[2:end]) / 2.0
        Δz_W[1] = Δz_W[2]
        Δz_W[end] = Δz_W[end-1]

        return new(
            Nz,
            z_T,
            z_W,
            Δz_T,
            Δz_W,
        )
 
    end
end



function loadGrid(
    filename :: String,
    z_W_varname :: String = "z_W"
)

    ds = Dataset(filename, "r")
    z_W = replace(ds[z_W_varname][:], missing=>NaN)
    close(ds)

    gd = Grid(
        z_W = z_W,
    )

    return gd
end

