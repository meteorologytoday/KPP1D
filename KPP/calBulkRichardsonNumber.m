

function [ BRi ] = calBulkRichardsonNumber(sop, zt, b_T, u_T, v_T)
    U_T = (u_T.^2 + v_T.^2).^0.5;
    dUdz_W = sop.W_ddz_T * U_T(:);
end





