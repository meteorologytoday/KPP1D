function [ BRi ] = calBulkRichardsonNumber(kpp, grid, sop, zt, b_T, u_T, v_T)

    
    db = b(1) - b;
    du = sqrt((u(1) - u).^2 + (v(1) - v).^2);
    Vtr_sqr = calUnresolvedShear(kpp, grid, sop, b, wb_sfc);
    
    U_T = (u_T.^2 + v_T.^2).^0.5;
    
    dUdz_W = sop.W_ddz_T * U_T(:);
end





