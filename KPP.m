classdef KPP < handle 
    properties
        c
    end
    methods

        function kpp = KPP()
            kpp.c = KPPConstants();
        end
        
        function [ u_star, L_star ] = calMOSTscales(kpp, tau0, wb_sfc)
            u_star = ( abs(tau0) / kpp.c.rho0 )^0.5;
            L_star = u_star^3 / ( - kpp.c.kappa * wb_sfc );
        end
        
        function [ Vtr_sqr ] = calUnresolvedShear(kpp, grid, sop, b, wb_sfc)

            N_osbl = sop.W_ddz_T * b;
            d = - grid.z_W;

            if (wb_sfc > 0)
                w_star = (d * wb_sfc).^(1/3);
            else
                w_star = d * 0;
            end

            Vtr_sqr = (-kpp.c.beta_T) * kpp.c.C_v / (kpp.c.Ri_c * kpp.c.kappa^(2/3) * (kpp.c.c_s * kpp.c.eps)^(1/6)) * (d .* N_osbl .* w_star);
            Vtr_sqr(Vtr_sqr < kpp.c.min_Vtr_sqr) = kpp.c.min_Vtr_sqr;

        end
        
        % LMD1994 Appendex B for scalars
        function [ phi_s ] = calPhi_s(kpp, d, L_star)
            zeta = d / L_star;
            if (L_star >= 0)
                phi_s = 1 + 5 * zeta;
            else
                % kpp.c.c_s = 98.96 in KPPconstants
                phi_s = (1 - 16*zeta).^(-1/2) .* (zeta > -1.0) + ...
                    (-28.86 - kpp.c.c_s * zeta).^(-1/3) .* (zeta <= -1.0);
            end
        end
        
        function [ w_s, sig, u_star, L_star ] = calw_s(kpp, z, h, tau0, wb_sfc)
            
            [ u_star, L_star ] = kpp.calMOSTscales(tau0, wb_sfc);
            d = -z;
            sig = d ./ h;

            phi_s = kpp.calPhi_s(d, L_star);
            w_s = kpp.c.kappa * u_star ./ phi_s;

            if (L_star < 0) % convective case, w_s topped at sig = eps
                for i=1:length(sig)
                    if (sig(i) >= kpp.c.eps)
                        phi_s = kpp.calPhi_s(h * kpp.c.eps, L_star);
                        w_s(i:end) = kpp.c.kappa * u_star / phi_s; 
                        break
                    end
                end
            end
        end
        
        
        function [ Ri, db, du_sqr, Vt_sqr ] = calBulkRichardsonNumber(kpp, grid, sop, wb_sfc, b, u, v)
            db = b(1) - b;
            du_sqr = (u(1) - u).^2 + (v(1) - v).^2;
            Vt_sqr = sop.T_interp_W * calUnresolvedShear(kpp, grid, sop, b, wb_sfc);

            Ri = (- grid.z_T) .* db ./ ( du_sqr + Vt_sqr );
        end
        
        % LMD94 equation (27)
        % Input on T-grid, output on W-grid
        function [ Ri_g ] = calGradientRichardsonNumber(kpp, grid, b, u, v)
            N = grid.sop.W_ddz_T * b;
            dudz = grid.sop.W_ddz_T * u;
            dvdz = grid.sop.W_ddz_T * v;
            
            Ri_g = N.^2 ./ (dudz.^2 + dvdz.^2);
        end
        
        % k is the index of deepest layer of mixed-layer 
        % h = grid.h_W(k+1);
        function [ h, k ] = calMixedLayerDepth(kpp, grid, Ri)
            
            % find Richardson number first exceeds kpp.c.Ri_c
            k = 0;
            for i=1:grid.Nz
                if (Ri(i) > kpp.c.Ri_c)
                    if (i == 1)
                        k = 1;  % min MLT needs at least one;
                    else
                        k = i - 1;
                    end
                    
                    break
                end
            end
            
            % Whole except of the last layer becomes ML
            if (k == 0)
                k = grid.Nz - 1;
            end
            
            h = grid.h_W(k+1);
        end
        
        % LMD94 equation (19) and (20)
        % Also the same as RAD18 equation (19)
        function flux = calNonlocalFlux_scalar(kpp, grid, h, tau0, wb_sfc, ws_sfc)
            flux = kpp.c.C_s * kpp.c.G(sig) * ws_sfc;
        end
        
        % Calculate scalar shape function in LMD equation (28)
        function S = shapeInterior(kpp, x)
            if x < 0
                S = 1;
            elseif x < 0.7
                S = (1 - (x / 0.7)).^3;
            else
                S = 0;
            end
        end
            
        % Calculate the interior diffusivity given by
        % LMD equation (28). Ri_g is given in equation (27)
        function K_s = calInteriorK_s(kpp, grid, b, u, v)
            Ri_g = calGradientRichardsonNumber(kpp, grid, b, u, v);
            K_s = 50e-4 * kpp.shapeInterior(Ri_g);
        end
        
        function K_s = calK_scalar(kpp, grid, h_k, tau0, wb_sfc)
            d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));
            h = grid.h_W(h_k + 1);
            
            
            
            [ w_s, sigma, ~, ~ ] = calw_s(kpp, grid.z_W, h, tau0, wb_sfc);
            
            G = kpp.c.calG(sigma);
            W_MLmask_W = d0( sigma <= 1 );
            
            K_s = W_MLmask_W * (h * G .* w_s) ;
        end
        
        function flux = calKPPFlux(kpp)
            flux = kpp.calNonlocalFlux() + kpp.calLocalFlux() + kpp.calInteriorFlux();
        end
        
    end
    
    
end