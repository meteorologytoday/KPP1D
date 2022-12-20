function calMOSTscale(
    τ0 :: Float64,
    B_f :: Float64,
)

    if τ0 == 0
        τ0 = 1e-10
    elseif τ0 < 0
        throw(ErrorException("Error: τ0 must be non-negative."))
    end
    
    u_star = ( abs(τ0) / ρ0 )^0.5
    
    if B_f == 0 # assume very weak wb_0 > 0
        L_star = - Inf;
    else
        L_star = u_star^3 / ( - κ * B_f )
    end

    return u_star, L_star

end

function calUnresolvedShear(
    b    :: AbstractArray{Float64},
    B_f  :: Float64,
    bmo  :: BasicMatrixOperator,
    gd   :: Grid,
)

    N_osbl = bmo.W_ddz_T * b
    N_osbl[N_osbl .< 0] = 0

    @. N_osbl = sqrt(N_osbl)

    d = - gd.z_W

    if B_f > 0
        w_star = (d * B_f).^(1/3)
    else
        w_star = d * 0
    end

    Vtr_sqr = β_T * C_v / (Ri_c * kappa^(2/3) * (c_s * ϵ)^(1/6)) * (d .* N_osbl .* w_star)

    Vtr_sqr[Vtr_sqr .< min_Vtr_sqr] = min_Vtr_sqr

    return Vtr_sqr

end

# LMD1994 Appendex B for scalars
function calϕ_s(
    d      :: Union{AbstractArray{Float64}, Float64},
    L_star :: Float64,
)

    ζ = d ./ L_star

    if L_star >= 0
        ϕ_s = 1 + 5 * ζ
    else
        # c_s = 98.96 in KPPconstants
        ϕ_s = (1 .- 16 * ζ).^(-1/2) .* (ζ .> -1.0) +
            (-28.86 .- c_s * ζ).^(-1/3) .* (ζ .<= -1.0)
    end

    return ϕ_s
end

# LMD1994 Appendex B for momentum
function calϕ_m(
    d      :: Union{AbstractArray{Float64}, Float64},
    L_star :: Float64,
)    
    ζ = d ./ L_star;
    if (L_star >= 0)
        ϕ_m = 1 + 5 * ζ
    else
        ϕ_m = (1 .- 16*ζ).^(-1/4) .* (ζ .> -0.2) + 
            (1.26 .- 8.38 * ζ).^(-1/3) .* (ζ .<= -0.2)
    end
end

#   function [ w_x, sig, u_star, L_star ] = calw_x(kpp, x, z, h, tau0, B_f)
function calw_x(
    x :: Symbol,
    z :: AbstractArray{Float64},
    h :: Float64,
    τ0 :: Float64,
    B_f :: Float64,
) 

    u_star, L_star = calMOSTscales(τ0, B_f)
    d = - z

    σ = d ./ h
           
    if x == :SCALAR
        ϕ_x = calϕ_s(d, L_star)
        w_x = κ * u_star ./ ϕ_x;
        w_x_bnd = κ * u_star ./ calϕ_s(h * ϵ, L_star) 
    elseif x == :MOMENTUM
        ϕ_x = calϕ_m(d, L_star)
        w_x = κ * u_star ./ ϕ_x
        w_x_bnd = κ * u_star ./ calϕ_m(h * ϵ, L_star)
    else
        throw(ErrorException("[calw_x] Unknown input $(x)"))
    end
    
    
    if L_star < 0 # In the case of convective condition, w_s is topped at σ = ϵ
        for i=1:length(σ)
            if σ[i] >= ϵ
                w_x[i:end] .= w_x_bnd; 
                break
            end
        end
    end

    return w_x, σ, u_star, L_star
end

function calBulkRichardsonNumber(
    B_f :: Float64,
    b   :: AbstractArray{Float64},
    u   :: Float64,
    v   :: Float64,
    bmo :: BasicMatrixOperator,
    gd  :: Grid,
)
    db = b[1] .- b
    du_sqr = (u[1] .- u).^2 + (v[1] .- v).^2
    Vt_sqr = bmo.T_interp_W * calUnresolvedShear(b, B_f, bmo, gd)

    Ri = (- gd.z_T) .* db ./ ( du_sqr + Vt_sqr )

    return Ri, db, du_sqr, Vt_sqr
end

# LMD94 equation (27)
# Input on T-grid, output on W-grid
function calGradientRichardsonNumber(
    b   :: AbstractArray{Float64},
    u   :: AbstractArray{Float64},
    v   :: AbstractArray{Float64},
    bmo :: BasicMatrixOperator,
    gd  :: Grid,
)

    N_sqr = bmo.W_ddz_T * b
    dudz = bmo.W_ddz_T * u
    dvdz = bmo.W_ddz_T * v
    gradU_sqr = dudz.^2 .+ dvdz.^2

    Ri_g = bmo.W_imask_W * ( N_sqr ./ gradU_sqr )

    Ri_g[gradU_sqr == 0] = 1e20 # or inifinity

    return Ri_g
end

# k is the index of deepest layer of mixed-layer 
# h = grid.h_W[k+1];
function calMixedLayerDepth(
    Ri :: AbstractArray{Float64},
    u_star :: Float64,
    L_star :: Float64,
    f      :: Float64,
    bmo :: BasicMatrixOperator,
    gd  :: Grid,
)
    
    # find Richardson number first exceeds Ri_c
    k = -1
    for i = 1:grid.Nz
        if Ri[i] > Ri_c
            
            if i == 1
                k = 1  # mixed-layer thickness cannot be zero
            else
                k = i - 1;
            end
            break
        end
    end
    
    # Whole except of the last layer becomes mix layer
    if k == -1  # cannot find the bottom of mix layer
        k = grid.Nz - 1
    end
    
    h = grid.d_W[k+1]
    
    
    # When stable forcing (B_f < 0 or L_star > 0) h cannot exceeds 
    # either L_star or Ekman depth as described in  LMD94 equation (24)
    if L_star > 0

        if f == 0.0
            h_E = Inf
        else
            h_E = 0.7 * u_star / abs(f)
        end
        
        h_max = min(h_E, L_star)

        if h > h_max
            for i=2:grid.Nz  # start from 2 because h is at least one level
                if grid.d_W[i+1] > h_max
                    k = i-1
                    h = grid.d_W[k+1]
                    break
                end
            end
        end
    end
    
    #fprintf('Mixed-layer depth: %f  , h_E =  %f \n', h, h_E);

    return h, k
end

# Calculate scalar shape function in LMD equation (28)
function shapeInterior(
    x :: AbstractArray{Float64},
)
    
    S = (1 .- (x ./ 0.7).^2).^3
    S[x < 0]    .= 1.0
    S[x >= 0.7] .= 0.0

    return S    
end
    
# Calculate the shear related interior diffusivity given by
# LMD equation (28). Ri_g is given in equation (27)
function calInteriorK_sh(
    b :: AbstractArray{Float64},
    u :: AbstractArray{Float64},
    v :: AbstractArray{Float64},
    bmo  :: BasicMatrixOperator,
    grid :: Grid,
)
    Ri_g = calGradientRichardsonNumber(b, u, v, bmo, grid)
    K_sh = 50e-4 * shapeInterior(Ri_g)
    return K_sh
end

function d0(v)
    return spdiagm(0 => view(v, :))
end

function calK_x(
    x :: Symbol,
    h_k :: Integer,
    τ0  :: Float64,
    B_f :: Float64,
    b :: AbstractArray{Float64},
    u :: AbstractArray{Float64},
    v :: AbstractArray{Float64},
    bmo  :: BasicMatrixOperator,
    grid :: Grid,
)
    
    h = grid.d_W[h_k + 1]
    w_x, σ, _, _ = calw_x(x, grid.z_W, h, τ0, B_f)
    
    G = calG(σ)
    ML_mask_arr = convert(Array{Float64}, σ <= 1)
    W_ML_mask_W = d0( ML_mask_arr )         # Mixed-layer  (σ <= 1)
    W_INT_mask_W = d0( 1.0 - ML_mask_arr )  # Interior     (σ >  1)
    
    K_x_ML  = W_ML_mask_W  * (h * G .* w_x) 
    K_x_INT = W_INT_mask_W * calInteriorK_sh(b, u, v, bmo, grid)

    return K_x_ML, K_x_INT
end

# LMD94 equation (19) and (20) for scalar
# Also the same as RAD18 equation (19)
function calNonLocalFlux_s(
    h_k  :: Integer,
    B_f  :: Float64,
    ws_0 :: Float64,
    bmo  :: BasicMatrixOperator,
    grid :: Grid,
)
    
    if B_f > 0 # unstable case, nonlocal flux is nonzero
        h    = grid.d_W[h_k+1]
        sig  = grid.d_W / h
        flux = C_s * calG.(sig) * ws_0
        flux[sig .> 1] .= 0.0
    else
        flux = zeros(Float64, grid.W_pts, 1)
    end
    
    return flux
end

# LMD94 equation (19) and (20) for momentum
# Also the same as RAD18 equation (19)
function calNonLocalFlux_m(
    h_k  :: Integer,
    tau0 :: Float64,
    bmo :: BasicMatrixOperator,
    grid :: Grid,
)
    flux = zeros(Float64, length(grid.W_pts), 1)

    return flux
end

    


