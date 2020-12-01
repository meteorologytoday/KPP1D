function [ grid ] = makeGrid(H, Nz)

    d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));

    T_pts = Nz;
    W_pts = Nz+1;

    dz = H / Nz;
    z_W = 0:-dz:-H;
    z_W = z_W(:); % make z_W into column form
    h_W = -z_W;
    
    if (any(h_W < 0))
        error('[Grid] h_W < 0 for at least one grid')
    end
    
    z_T = (z_W(1:end-1) + z_W(2:end)) / 2;
    dz_W = z_W .* 0 + dz;
    dz_T = z_T .* 0 + dz;
    
    % grid operators
    mask_T = ones(T_pts, 1);
    mask_W = ones(W_pts, 1);
    
    T_mask_T = d0(mask_T);
    W_mask_W = d0(mask_W);
    
    imask_W = mask_W * 1;  imask_W(1) = 0; imask_W(end) = 0;
    W_imask_W = d0(imask_W);
    
    T_num = 1:T_pts;
    W_num = 1:W_pts;
    
    T_I_T = speye(T_pts);
    W_I_W = speye(W_pts);
    
    % send upward
    T_UP_W = W_I_W([2:W_pts], :);
    W_UP_T = T_I_T([1:T_pts], :);   W_UP_T(end+1, :) = 0;
    
    % send downward
    W_DN_T = T_UP_W';
    T_DN_W = W_UP_T';
    
    T_UP_T = T_UP_W * W_UP_T;
    W_UP_W = W_UP_T * T_UP_W;
    
    T_DN_T = T_DN_W * W_DN_T;
    W_DN_W = W_DN_T * T_DN_W;
    
    grid.Nz    = Nz;
    grid.T_pts = T_pts;
    grid.W_pts = W_pts;
    
    grid.dz = dz;
    grid.z_W = z_W;
    grid.z_T = z_T;
    grid.dz_W = dz_W;
    grid.dz_T = dz_T;
    
    grid.h_W = h_W;
    
    grid.T_I_T = T_I_T;
    grid.W_I_W = W_I_W;
    grid.T_mask_T = T_mask_T;
    grid.W_mask_W = W_mask_W;
    grid.W_imask_W = W_imask_W;
    
    grid.T_UP_W = T_UP_W;
    grid.W_UP_T = W_UP_T;
    grid.W_DN_T = W_DN_T;
    grid.T_DN_W = T_DN_W;
    
    grid.T_UP_T = T_UP_T;
    grid.W_UP_W = W_UP_W;
    grid.T_DN_T = T_DN_T;
    grid.W_DN_W = W_DN_W;
    
    grid.sop = makeSpatialOperators(grid);
end
