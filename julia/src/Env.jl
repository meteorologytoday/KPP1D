mutable struct Env

    gd :: Grid
    f  :: Float64
    mask_T :: AbstractArray{Float64, 1}

    function Env(;
        gd :: Grid,
        f  :: Float64, # Coriolis parameter
        mask_T :: Union{Nothing, AbstractArray} = nothing,
    )

        if mask_T == nothing
            mask_T = ones(Float64, length(gd.z_T))
        end

        return new(
            gd,
            f,
            mask_T,
        )
        
    end
end
