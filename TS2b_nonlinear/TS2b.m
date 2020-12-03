% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.
function [ b ] = TS2b(T, S)

    c = TSConstants();

    dT = T - c.T_ref;
    dS = S - c.S_ref;
    b = - c.g * ( c.rho1*dT + c.rho2*dS +  c.rho3*(dT.^2) + c.rho4*(dS.^2) ...
        + c.rho5* (dT.*dS) + c.rho6*(dT.^3) + c.rho7*(dS.^2).*dT + c.rho8*(dT.^2).*dS + c.rho9*(dS.^3) );
    
    
end