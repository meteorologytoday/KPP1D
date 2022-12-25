mutable struct State

    T  :: AbstractArray{Float64, 1}
    S  :: AbstractArray{Float64, 1}
    b  :: AbstractArray{Float64, 1}
    u  :: AbstractArray{Float64, 1}
    v  :: AbstractArray{Float64, 1}
    Ri :: AbstractArray{Float64, 1}
    h   :: Float64
    h_k :: Integer

    function State(
        ev :: Env,
    )
        Nz = ev.gd.Nz
 
        T   = zeros( Float64, Nz ) 
        S   = zeros( Float64, Nz ) 
        b   = zeros( Float64, Nz ) 
        Ri  = zeros( Float64, Nz ) 
        u   = zeros( Float64, Nz ) 
        v   = zeros( Float64, Nz ) 
    

        h_k = 1
        h   = ev.gd.z_W[h_k + 1]

        return new(
            T,
            S,
            b,
            u,
            v,
            Ri,
            h,
            h_k,
        )    
    end
end



