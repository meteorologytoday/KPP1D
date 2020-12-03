% As suggested by LMD94, we use the same coefficient
% given by Paulson & Simpson (1977)
classdef Radiation < handle 
    properties
        N     % bands of radiation
        r      % fraction of radiation in total irradiance
        mu_inv % penetration depth
        coe_flux_W
        coe_fluxconv_T
        coe_total_flux_W
        coe_total_fluxconv_T
        coe_turbulent_flux_T
    end
    methods

        function o = Radiation(grid)
            r = [0.58 ; 0.42];
            mu_inv = [0.35 ; 23];
            N = length(r);
            
            coe_flux_W    = zeros(grid.W_pts, N);
            coe_fluxconv_T = zeros(grid.T_pts, N);
            
            for n=1:N
                coe_flux_W(:, n) = r(n) * exp( grid.z_W / mu_inv(n) );
                coe_flux_W(end, n) = 0; % no flux penetrates into the ground
                coe_fluxconv_T(:, n) = - grid.sop.T_ddz_W * coe_flux_W(:, n);
            end
            
            o.r = r;
            o.N = N;
            o.mu_inv = mu_inv;
            o.coe_flux_W    = coe_flux_W;
            o.coe_fluxconv_T = coe_fluxconv_T;
            o.coe_total_flux_W     = sum(coe_flux_W, 2);
            o.coe_total_fluxconv_T = sum(coe_fluxconv_T, 2);
            o.coe_turbulent_flux_T = o.coe_total_flux_W(1) - o.coe_total_flux_W(2:end);
        end
        
        function [ F, Q ] = calRadiation(rad, I0)
            F = I0 * rad.coe_total_flux_W;
            Q = I0 * rad.coe_total_fluxconv_T;
        end
    end
end