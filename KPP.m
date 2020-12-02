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
        
        function [ Vtr_sqr ] = calUnresolvedShear(kpp, grid, b, wb_sfc)

            N_osbl = grid.sop.W_ddz_T * b;
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
            
            % surface phi_s(1) = 0
            w_s(1) = 1;
            
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
        
        function [ Ri, db, du_sqr, Vt_sqr ] = calBulkRichardsonNumber(kpp, grid, wb_sfc, b, u, v)
            db = b(1) - b;
            du_sqr = (u(1) - u).^2 + (v(1) - v).^2;
            Vt_sqr = grid.sop.T_interp_W * calUnresolvedShear(kpp, grid, b, wb_sfc);

            Ri = (- grid.z_T) .* db ./ ( du_sqr + Vt_sqr );
        end
        
        % LMD94 equation (27)
        % Input on T-grid, output on W-grid
        function [ Ri_g ] = calGradientRichardsonNumber(kpp, grid, b, u, v)
            N = grid.sop.W_ddz_T * b;
            dudz = grid.sop.W_ddz_T * u;
            dvdz = grid.sop.W_ddz_T * v;
            gradU_sqr = dudz.^2 + dvdz.^2;
            Ri_g = grid.W_imask_W * ( N.^2 ./ gradU_sqr );
            Ri_g(gradU_sqr == 0) = 1e20; % or inifinity
        end
        
        % k is the index of deepest layer of mixed-layer 
        % h = grid.h_W(k+1);
        function [ h, k ] = calMixedLayerDepth(kpp, grid, Ri, u_star, f)
            
            % find Richardson number first exceeds kpp.c.Ri_c
            k = 0;
            for i=1:grid.Nz
                if (Ri(i) > kpp.c.Ri_c)
                    %fprintf('%d => %f\n', i, Ri(i));
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
            
            h = grid.d_W(k+1);
            
            % h cannot exceeds Ekman depth as in
            % LMD94 equation (24)
            h_E = (0.7 * u_star / f);
            if (h > h_E && h < 0)
                for i=2:grid.Nz  % start from 2 because h is at least one level
                    if grid.d_W(i+1) > h_E
                        k = i-1;
                        h = grid.d_W(k+1);
                        break;
                    end
                end
            end
            
            %fprintf('Mixed-layer depth: %f  , h_E =  %f \n', h, h_E);
        end
        
        % Calculate scalar shape function in LMD equation (28)
        function S = shapeInterior(kpp, x)
            
            S = (1 - (x / 0.7).^2).^3;
            S(x < 0) = 1;
            S(x >= 0.7) = 0;
            
        end
            
        % Calculate the interior diffusivity given by
        % LMD equation (28). Ri_g is given in equation (27)
        function K_s = calInteriorK_s(kpp, grid, b, u, v)
            Ri_g = calGradientRichardsonNumber(kpp, grid, b, u, v);
            K_s = 50e-4 * kpp.shapeInterior(Ri_g);
        end
        
        function [ K_s_ML, K_s_INT ]  = calK_s(kpp, grid, h_k, tau0, wb_sfc, b, u, v)
            d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));
            h = grid.d_W(h_k + 1);
           
            [ w_s, sigma, ~, ~ ] = calw_s(kpp, grid.z_W, h, tau0, wb_sfc);
            
            G = kpp.c.calG(sigma);
            W_ML_mask_W = d0( sigma <= 1 ); % Mixed-layer
            W_INT_mask_W = d0( sigma > 1 ); % Interior
            
            K_s_ML = W_ML_mask_W * (h * G .* w_s);     
            K_s_INT = W_INT_mask_W * calInteriorK_s(kpp, grid, b, u, v);
            
        end
        
        % LMD94 equation (19) and (20) for scalar
        % Also the same as RAD18 equation (19)
        function flux = calNonLocalFlux_s(kpp, grid, h_k, wb_sfc, ws_sfc)
            
            if (wb_sfc > 0) % unstable case, nonlocal flux is nonzero
                h = grid.d_W(h_k+1);
                sig = grid.d_W / h;
                flux = kpp.c.C_s * kpp.c.calG(sig) * ws_sfc;
                flux(sig > 1) = 0;
            else
                flux = zeros(grid.W_pts, 1);
            end
            
        end
        
        % LMD94 equation (19) and (20) for momentum
        % Also the same as RAD18 equation (19)
        function flux = calNonLocalFlux_m(kpp, grid, h_k, tau0)
            flux = zeros(length(grid.W_pts), 1);
        end
        
        
    end
    
    
end