function [ o ] = KPPConstants()

    o.Ri_c  = 0.3; % Critical Richardson number
    o.kappa = 0.4; % von Karman constant
    o.eps = 0.1;  % epsilon in LMD1994
    o.c_s = 98.96;
    o.calG = @(x) x .* ( 1 - x ).^2;
    o.rho0 = 1026;
    o.beta_T = -0.2;
    o.C_v = 1.0;
    o.min_Vtr_sqr = 1e-10;

end