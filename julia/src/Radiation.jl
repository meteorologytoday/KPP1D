# As suggested by LMD94, we use the same coefficient
# given by Paulson & Simpson (1977)

# The caculation of Raidation follows the description of
# Appendix A of Large et al. (1994)

mutable struct Radiation

    N     :: Integer                    # bands of radiation
    r     :: AbstractArray{Float64, 1}  # fraction of radiation in total irradiance
    μ_inv :: AbstractArray{Float64, 1}  # penetration depth

    coe_flux_W           :: AbstractArray{Float64, 2}
    coe_fluxconv_T       :: AbstractArray{Float64, 2}
    coe_total_flux_W     :: AbstractArray{Float64, 1}
    coe_total_fluxconv_T :: AbstractArray{Float64, 1}
    coe_turbulent_flux_T :: AbstractArray{Float64, 1}

    function Radiation(;
        amo :: AdvancedMatrixOperators,
    )
        r     = [0.58, 0.42]  # This array should sum to 1.0
        μ_inv = [0.35, 23.0]
        
        #r     = [0.00, 1.00] # for testing
        #μ_inv = [0.35, 10.0] # for testing

        N     = length(r)
        
        gd = amo.gd
        coe_flux_W     = zeros(Float64, amo.bmo.W_pts, N)
        coe_fluxconv_T = zeros(Float64, amo.bmo.T_pts, N)
        
        for n = 1:N
            coe_flux_W[:, n]     = r[n] .* exp.( gd.z_W ./ μ_inv[n] )
            coe_flux_W[end, n]   = 0.0 # no flux penetrates into the ground
            coe_fluxconv_T[:, n] = - amo.T_DIVz_W * coe_flux_W[:, n]
        end

        # Total flux across different type of radiation (i.e. longwave and shortwave)
        coe_total_flux_W     = sum(coe_flux_W,     dims=2)[:, 1]
        coe_total_fluxconv_T = sum(coe_fluxconv_T, dims=2)[:, 1]

        # Equation (A4) of Large et al. 1994
        # Basically this is a pre-calculated fraction of radiation
        # that is used as the surface non-local forcing within
        # the ocean boundary layer h
        coe_turbulent_flux_T = coe_total_flux_W[1] .- coe_total_flux_W[2:end]
      
        return new(
            N,
            r,
            μ_inv,
            coe_flux_W,
            coe_fluxconv_T,
            coe_total_flux_W,
            coe_total_fluxconv_T,
            coe_turbulent_flux_T,
        ) 
    end
end
 
function calRadiation(
    rad :: Radiation,
    I0  :: Float64,
)

    F = I0 * rad.coe_total_flux_W
    Q = I0 * rad.coe_total_fluxconv_T
    
    return F, Q

end
