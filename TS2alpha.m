% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.

% Thermal expansion coefficient £\ := (drho/dT)_S
function [ alpha ] = TS2alpha(T, S)
    c = TSConstants();
    dT = T - c.T_ref;
    dS = S - c.S_ref;
    alpha = - ( c.rho1 + 2*c.rho3*dT + c.rho5*dS + 3*c.rho6*(dT.^2) + c.rho7*(dS.^2) + 2*c.rho8.*dS.*dT );
end
