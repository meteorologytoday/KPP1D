% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.
function [ b ] = TS2b(T, S)

    c = constants();

    dT = T - c.T_ref;
    dS = S - c.S_ref;
    b = - c.g * ( c.rho1*dT + c.rho2*dS +  c.rho3*(dT.^2) + c.rho4*(dS.^2) ...
        + c.rho5* (dT.*dS) + c.rho6*(dT.^3) + c.rho7*(dS.^2).*dT + c.rho8*(dT.^2).*dS + c.rho9*(dS.^3) );
    
end

% Thermal expansion coefficient £\ := (?rho/?T)_S
function [ alpha ] = TS2alpha(T, S)
    dT = T - T_ref;
    dS = S - S_ref;
    alpha = - ( rho1 + 2*rho3*dT + rho5*dS + 3*rho6*(dT^2) + rho7*(dS^2) + 2*rho8*dS*dT );
end

% Salinity expansion coefficient £] := - (?rho/?S)_T  (notice the negative sign)
function [ beta ] = TS2beta(T, S)
    dT = T - T_ref;
    dS = S - S_ref;
    beta = rho2 + 2*rho4*dS + rho5*dT + 2*rho7*dS*dT + rho8*(dT^2) + 3*rho9*(dS^2);
end

function [ c ] = constants()
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