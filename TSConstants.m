% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.
function [ c ] = TSConstants()
    c.g      = 9.80616;    % m / s^2      copied from models/csm_share/shr/shr_const_mod.F90

    c.T_ref =   13.5;
    c.S_ref =   32.6;
    c.rho_ref = 1024.458;

    c.rho1 = -.20134    / c.rho_ref;
    c.rho2 =  .77096    / c.rho_ref;
    c.rho3 = -.49261e-2 / c.rho_ref;
    c.rho4 =  .46092e-3 / c.rho_ref;
    c.rho5 = -.20105e-2 / c.rho_ref;
    c.rho6 =  .36597e-4 / c.rho_ref;
    c.rho7 =  .47371e-5 / c.rho_ref;
    c.rho8 =  .37735e-4 / c.rho_ref;
    c.rho9 =  .65493e-5 / c.rho_ref;
end