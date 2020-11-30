function [ o ] = KPPConstants()

    o.kappa = 0.4; % von Karman constant
    o.eps = 0.1;  % epsilon in LMD1994
    o.c_s = 1;
    o.calG = @(x) x .* ( 1 - x ).^2;
    o.rho0 = 1026;
    o.beta_T = -0.2;
    o.c_s = 98.96;

end