function [ Ri, db, du_sqr, Vt_sqr ] = calBulkRichardsonNumber(kpp, grid, sop, wb_sfc, b, u, v)

    db = b(1) - b;
    du_sqr = (u(1) - u).^2 + (v(1) - v).^2;
    Vt_sqr = sop.T_interp_W * calUnresolvedShear(kpp, grid, sop, b, wb_sfc);
    
    %Ri = (- grid.z_T) .* db ./ ( du_sqr + Vt_sqr );
    Ri = (- grid.z_T) .* db./ ( du_sqr + Vt_sqr );

end





