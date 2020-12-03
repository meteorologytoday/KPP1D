classdef Const < handle 
    properties
        rho_sw = 1026
        cp_sw  = 3996
        g      = 9.80616;    % m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
        cprho_sw
    end
    methods

        function o = Const()
            o.cprho_sw = o.rho_sw * o.cp_sw;
        end
       
    end
end