Base.@kwdef mutable struct Forcing

    B_f :: Float64   = 0.0

    wT_0 :: Float64  = 0.0
    wS_0 :: Float64  = 0.0
    wb_0 :: Float64  = 0.0
    
    wu_0 :: Float64  = 0.0
    wv_0 :: Float64  = 0.0
    
    U10  :: Float64  = 0.0 # surface 10m wind speed on x direction
    V10  :: Float64  = 0.0 # surface 10m wind speed on y direction
    
    I_0  :: Float64     = 0.0 # solar radiation,     positive = upward
    Hf_sen  :: Float64  = 0.0 # sensible heat flux,  positive = upward
    Hf_lat  :: Float64  = 0.0 # latent heat flux,    positive = upward
    Hf_lw   :: Float64  = 0.0 # longwave radiation heat flux. This is the radiation from the blackbody radiation rather than solar irradiance I_0

    precip :: Float64 = 0.0 # m / s
    evap   :: Float64 = 0.0 # m / s
    τ0     :: Float64 = 0.0 # Pa
    τx0    :: Float64 = 0.0 # Pa
    τy0    :: Float64 = 0.0 # Pa
    albedo :: Float64 = 0.06 #
    T_a    :: Float64 = 0.0
    q_a    :: Float64 = 0.0


    surf_flux_calculation_type :: Symbol = :REALISTIC

end
