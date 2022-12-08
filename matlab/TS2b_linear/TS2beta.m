% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.

% Salinity expansion coefficient £] := - (drho/dS)_T  (notice the negative sign)
function [ beta ] = TS2beta(T, S)
    c = TSConstants();
    beta = c.beta;
end