% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.

% Salinity expansion coefficient £] := - (drho/dS)_T  (notice the negative sign)
function [ beta ] = TS2beta(T, S)
    c = TSConstants();
    dT = T - c.T_ref;
    dS = S - c.S_ref;
    beta = c.rho2 + 2*c.rho4*dS + c.rho5*dT + 2*c.rho7.*dS.*dT + c.rho8*(dT.^2) + 3*c.rho9*(dS.^2);
end