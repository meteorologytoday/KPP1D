% Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
% Here we use Table 3, and arbitrary pick Z=0m as the formula.

% Thermal expansion coefficient £\ := (drho/dT)_S
function [ alpha ] = TS2alpha(T, S)
    c = TSConstants();
    alpha = c.alpha;
end
