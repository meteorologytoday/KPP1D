% Using Van Roekel et al. (2018) equation (22)~(23)
function [ c ] = TSConstants()
    c.g      = 9.80616;    % m / s^2      copied from models/csm_share/shr/shr_const_mod.F90

    c.T_ref =  298.15 - 273.15;
    c.S_ref =  35.0;
    c.alpha = 2e-4;
    c.beta  = 8e-4;
end