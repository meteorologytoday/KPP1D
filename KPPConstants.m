function [ o ] = KPPConstants()

    c = Const();
    o.Ri_c  = 0.3; % Critical Richardson number
    o.kappa = 0.4; % von Karman constant
    o.eps = 0.1;  % epsilon in LMD1994
    o.c_s = 98.96;  % LMD94 equation (B1) and (B2)
    o.c_m = 8.38;  % LMD94 equation (B1) and (B2)
    

    o.calG = @(x) x .* ( 1 - x ).^2;
    o.rho0 = c.rho_sw;  % seawater density
    o.beta_T = -0.2;
    o.C_v = 1.6;
    
    % Mininum unresolved shear.
    % However, although it is seen in the Cvmix project source code but
    % I cannot find the actual value anywhere in either LMD94 or Van Roekel et al. (2018)
    o.min_Vtr_sqr = 1e-4; 

    o.C_star = 10; % LMD94 equation (20)
    o.C_s = o.C_star * o.kappa * (o.c_s * o.kappa * o.eps)^(1/3); % LMD94 equation (20)

    o.SCALAR = 1;
    o.MOMENTUM = 2;
end