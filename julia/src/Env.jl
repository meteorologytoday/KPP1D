mutable struct Env

    gd :: Grid

    function Env(
        gd :: Grid,
        f  :: Float64, # Coriolis parameter
    )
        
        return new(
            gd,
        )
        
    end
end
