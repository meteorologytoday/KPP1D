function [ op ] = mkOp_diffz(grid, sop, b, Kv_iso, Kv_cva)
    
    d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));
    
    Kv_W = zeros(grid.W_pts, 1) + Kv_iso;

    db = (grid.W_DN_T - grid.W_UP_T) * b;
    Kv_W( db < 0 ) = Kv_cva;
    Kv_W( db >= 0 ) = Kv_iso;
    
    op = sop.T_ddz_W * grid.W_imask_W * d0(Kv_W) * sop.W_ddz_T;
    %op = Kv_iso * sop.T_ddz_W * grid.W_mask_W * sop.W_ddz_T;
    
end