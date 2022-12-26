mutable struct AdvancedMatrixOperators

    bmo       :: BasicMatrixOperators
    gd        :: Grid

    T_DIVz_W    :: AbstractArray{Float64, 2}
    
    W_∂z_T      :: AbstractArray{Float64, 2}
    T_∂z_W      :: AbstractArray{Float64, 2}
    
    W_interp_T :: AbstractArray{Float64, 2}  # interpolation of T grid onto W grid
    T_interp_W :: AbstractArray{Float64, 2}  # interpolation of W grid onto T grid

    T_mask_T       :: AbstractArray{Float64, 2}
    W_mask_W       :: AbstractArray{Float64, 2}
    #W_imask_W       :: AbstractArray{Float64, 2}  # internal mask

    T_Δz_T :: AbstractArray{Float64, 2}
    W_Δz_W :: AbstractArray{Float64, 2}
 
    T_invΔz_T :: AbstractArray{Float64, 2}
    W_invΔz_W :: AbstractArray{Float64, 2}
 
    
    function AdvancedMatrixOperators(;
        gd             :: Grid,
        mask_T         :: AbstractArray{Float64, 1},
        bmo :: Union{Nothing, BasicMatrixOperators} = nothing,
    )

        # define a converter to make 2D variable repeat in z direction for Nz times
        cvt_diagm = (x,) -> spdiagm( 0 => view(x, :) )

        Nz = gd.Nz
        if bmo == nothing 
            println("Construct BMO")
            @time bmo = BasicMatrixOperators(Nz=Nz)
        end
 
        if length(mask_T) != bmo.T_pts
            throw(ErrorException("Length of topo.mask_T does not conform"))
        end

        onW_if_unblocked_up_onT    = bmo.W_DN_T * mask_T
        onW_if_unblocked_dn_onT    = bmo.W_UP_T * mask_T

        W_mask = onW_if_unblocked_up_onT .* onW_if_unblocked_dn_onT

        T_mask_T = spdiagm(0 => mask_T)
        W_mask_W = spdiagm(0 => W_mask)

        # ===== [ BEGIN ] =====

        T_Δz_T = (gd.Δz_T |>  cvt_diagm)
        W_Δz_W = (gd.Δz_W |>  cvt_diagm)
 
        T_invΔz_T = (gd.Δz_T.^(-1) |>  cvt_diagm)
        W_invΔz_W = (gd.Δz_W.^(-1) |>  cvt_diagm)

        # ===== [ END ] =====
        
        #println("Making derivatives") 

        # ===== [ BEG making matrix ] =====
        # MAGIC!!

        T_DIVz_W = T_mask_T * T_invΔz_T * (bmo.T_DN_W - bmo.T_UP_W )    ; dropzeros!(T_DIVz_W);
        W_∂z_T  = W_mask_W * W_invΔz_W * (bmo.W_DN_T - bmo.W_UP_T)      ; dropzeros!(W_∂z_T);
        T_∂z_W  = T_mask_T * T_invΔz_T * ( bmo.T_DN_W - bmo.T_UP_W )    ; dropzeros!(T_∂z_W);

        function selfDivision(m, ones_vec)
            local wgts = m * ones_vec
            m_t = transpose(m) |> sparse
            for (i, wgt) in enumerate(wgts)
                if wgt != 0
                    _beg = m_t.colptr[i]
                    _end = m_t.colptr[i+1]-1
                    m_t.nzval[_beg:_end] ./= wgt
                end
            end
          
            return dropzeros(transpose(m_t) |> sparse)
        end

        #println("Making interpolations part 1") 
        ones_T  = ones(Float64, bmo.T_pts)
        W_interp_T = (bmo.W_DN_T + bmo.W_UP_T) * T_mask_T
        W_interp_T = selfDivision(W_interp_T, ones_T)

        ones_W  = ones(Float64, bmo.W_pts)
        T_interp_W = (bmo.T_DN_W + bmo.T_UP_W) * W_mask_W
        T_interp_W = selfDivision(T_interp_W, ones_W)


        return new(

            bmo,
            gd,

            T_DIVz_W,
            
            W_∂z_T,
            T_∂z_W,
            
            W_interp_T,
            T_interp_W,

            T_mask_T,
            W_mask_W,

            T_Δz_T,
            W_Δz_W,
         
            T_invΔz_T,
            W_invΔz_W,
            
        )
    end
end
