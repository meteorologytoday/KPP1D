function [ Vtr_sqr ] = calUnresolvedShear(kpp, grid, sop, b, wb_sfc)
    
    N_osbl = sop.W_ddz_T * b;
    d = - grid.zw;
    w_star = (- d * abs(wb_sfc)).^(1/3) * sign(wb_sfc) ;

    Vtr_sqr = (-kpp.beta_T) * kpp.C_v / (kpp.Ri_c * kpp.kappa^(2/3) * (kpp.c_s * kpp.eps)^(1/6)) * (d .* N_osbl .* w_star);
        
    Vtr_sqr(Vtr_sqr < kpp.min_Vtr_sqr) = kpp.min_Vtr_sqr;
    
end
