classdef Const < handle 
    properties
        g        = 9.80616;    % m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
        Lv_w     = 2.257e6;    % Heat of evaporation of water in J/kg
        rho_fw   = 1000;
        rho_sw   = 1026
        cp_sw    = 3996
        rho_a    = 1.22
        cp_a     = 1004
        sigma_sb = 5.670374419e-8 % Stefan-Boltzmann constant
        Cel_Kel_offset = 273.15
        cprho_sw
        cprho_a
        
    end
    methods

        function o = Const()
            o.cprho_sw = o.rho_sw * o.cp_sw;
            o.cprho_a = o.rho_a * o.cp_a;
        end
       
    end
end