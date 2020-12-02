classdef Const < handle 
    properties
        rho_sw = 1026
        cp_sw  = 3996
        g      = 9.80616;    % m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
    end
    methods

        function o = Const()

        end
       
    end
end