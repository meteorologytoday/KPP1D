function [ Ks ] = calKPP(kpp, z, h, tau0, wb_sfc)
    [ w_s, sig, ~, ~ ]  = calw_s(kpp, z, h, tau0, wb_sfc);
    G = kpp.calG(sig);
    Ks = h * ( w_s .* G );
end



