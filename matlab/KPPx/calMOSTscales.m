function [ u_star, L_star ] = calMOSTscales(kpp, tau0, wb_sfc)
    u_star = ( abs(tau0) / kpp.rho0 )^0.5;
    L_star = u_star^3 / ( - kpp.kappa * wb_sfc );
end