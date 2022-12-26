using SparseArrays

@inline function speye(dtype, n)
    return spdiagm(0=>ones(dtype, n))
end

# Assuming x-direction is periodic
struct BasicMatrixOperators

    T_dim
    W_dim

    T_pts
    W_pts

    # Nomenclature:
    #
    # [new-grid][direction][old-grid]
    #
    # U_W_T : sending variable westward from T grid to U grid

    T_I_T
    W_I_W

    T_UP_T
    T_DN_T
 
    T_UP_W
    T_DN_W

    W_UP_T
    W_DN_T

    function BasicMatrixOperators(;
        Nz             :: Int64,
    )

        T_dim = (Nz, )
        W_dim = (Nz+1, )

        T_pts = reduce(*, T_dim)
        W_pts = reduce(*, W_dim)

        T_I_T = speye(Float64, T_pts)
        W_I_W = speye(Float64, W_pts)

        T_I_T_expand = vcat(T_I_T, zeros(Float64, 1, T_pts))
        W_I_W_expand = vcat(W_I_W, zeros(Float64, 1, W_pts))


        num_T = zeros(Int64, T_dim...)
        num_W = zeros(Int64, W_dim...)

        num_T[:] = 1:length(num_T)
        num_W[:] = 1:length(num_W)
        
        T = num_T * 0
        W = num_W * 0

        #smb = SparseMatrixBuilder(Nx*(Ny+1)*(Nz+1)*4)
        function build!(id_mtx, idx; wipe=:none)
           #println("Build!")
            local result
            rows = size(id_mtx)[1]
            if wipe == :t
                idx[1] = rows
            elseif wipe == :b
                idx[end] = rows
            elseif wipe != :none
                throw(ErrorException("Wrong keyword"))
            end
           
            # using transpose speeds up by 100 times 
            tp = transpose(id_mtx) |> sparse
            result = transpose(tp[:, view(idx, :)]) |> sparse
            #result = id_mtx[view(idx, :), :]
            #dropzeros!(result)

            idx .= 0 # clean so that debug is easier when some girds are not assigned
            return result
        end

        # upward, downward passing mtx
        T[1:Nz-1] = view(num_T, 2:Nz,  );    T_UP_T = build!(T_I_T_expand, T; wipe=:b)
        T[2:Nz  ] = view(num_T, 1:Nz-1,);    T_DN_T = build!(T_I_T_expand, T; wipe=:t)

        T[:     ] = view(num_W, 2:Nz+1,);    T_UP_W = build!(W_I_W_expand, T)
        T[:     ] = view(num_W, 1:Nz  ,);    T_DN_W = build!(W_I_W_expand, T)

        # inverse directions
        W_DN_T = T_UP_W' |> sparse
        W_UP_T = T_DN_W' |> sparse

        return new(
            T_dim, W_dim,
            T_pts, W_pts,
            T_I_T, W_I_W,

            T_UP_T, T_DN_T,
            T_UP_W, T_DN_W,
            W_UP_T, W_DN_T,
        )

    end
end
