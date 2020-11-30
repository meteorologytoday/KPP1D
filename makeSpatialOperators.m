function [ sop ] = makeSpatialOperators(grid)

    d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));

    W_ddz_T = grid.W_imask_W * d0(1./grid.dz_W) * (grid.W_DN_T - grid.W_UP_T);
    T_ddz_W = d0(1./grid.dz_T) * (grid.T_DN_W - grid.T_UP_W);
    
    T_ddz2_T = T_ddz_W * W_ddz_T;
    W_interp_T = grid.W_imask_W * (grid.W_DN_T + grid.W_UP_T) / 2.0;
    T_interp_W = (grid.T_DN_W + grid.T_UP_W) / 2.0;
    
    
    sop.W_ddz_T = W_ddz_T;
    sop.T_ddz_W = T_ddz_W;
    sop.T_ddz2_T = T_ddz2_T;
    sop.W_interp_T = W_interp_T;
    sop.T_interp_W = T_interp_W;
end
