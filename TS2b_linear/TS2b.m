% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.
function [ b ] = TS2b(T, S)

    c = TSConstants();

    dT = T - c.T_ref;
    dS = S - c.S_ref;
    b = - c.g * (1 - c.alpha * dT + c.beta * dS);
    
end