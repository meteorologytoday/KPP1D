function [ w_s, sig, u_star, L_star ] = calw_s(kpp, z, h, tau0, wb_sfc)
    [ u_star, L_star ] = calMOSTscales(kpp, tau0, wb_sfc);
    d = -z;
    sig = d ./ h;
    
    phi_s = calPhi_s(d, L_star);
    w_s = kpp.kappa * u_star ./ phi_s;
    
    if (L_star < 0) % convective case, w_s topped at sig = eps
        for i=1:length(sig)
            if (sig(i) >= kpp.eps)
                phi_s = calPhi_s(h * kpp.eps, L_star);
                w_s(i:end) = kpp.kappa * u_star / phi_s; 
                break
            end
        end
    end
end