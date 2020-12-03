classdef KPP < handle 
    properties
        c
    end
    methods

        function kpp = KPP()
            kpp.c = KPPConstants();
        end
        
        function [ u_star, L_star ] = calMOSTscales(kpp, tau0, B_f)
            
            if (tau0 == 0)
                tau0 = 1e-10;
            elseif (tau0 < 0)
                error('tau0 must be equal or greater than zero');
            end
            
            u_star = ( abs(tau0) / kpp.c.rho0 )^0.5;
            L_star = u_star^3 / ( - kpp.c.kappa * B_f );
            
            if (B_f == 0)  % assume very weak wb_b0 > 0
                L_star = - inf;
            end
            
        end
        
        function [ Vtr_sqr ] = calUnresolvedShear(kpp, grid, b, B_f)

            N_osbl = (grid.sop.W_ddz_T * b).^0.5;
            d = - grid.z_W;

            if (B_f > 0)
                w_star = (d * B_f).^(1/3);
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
        
        % LMD1994 Appendex B for momentum
        function [ phi_m ] = calPhi_m(kpp, d, L_star)
            zeta = d / L_star;
            if (L_star >= 0)
                phi_m = 1 + 5 * zeta;
            else
                phi_m = (1 - 16*zeta).^(-1/4) .* (zeta > -0.2) + ...
                    (1.26 - 8.38 * zeta).^(-1/3) .* (zeta <= -0.2);
            end
        end
        
        function [ w_x, sig, u_star, L_star ] = calw_x(kpp, x, z, h, tau0, B_f)
            
            [ u_star, L_star ] = kpp.calMOSTscales(tau0, B_f);
            d = -z;
            sig = d ./ h;
                   
            if (x == kpp.c.SCALAR)
                fprintf('L_star = %f\n', L_star);
                phi_x = kpp.calPhi_s(d, L_star);
                w_x = kpp.c.kappa * u_star ./ phi_x;
                w_x_bnd = kpp.c.kappa * u_star / kpp.calPhi_s(h * kpp.c.eps, L_star); 
            elseif (x == kpp.c.MOMENTUM)
                phi_x = kpp.calPhi_m(d, L_star);
                w_x = kpp.c.kappa * u_star ./ phi_x;
                w_x_bnd = kpp.c.kappa * u_star / kpp.calPhi_m(h * kpp.c.eps, L_star); 
            else
                disp(x);
                error('[calw_x] Unknown input');
            end
            
            
            if (L_star < 0) % convective case, w_s topped at sig = eps
                for i=1:length(sig)
                    if (sig(i) >= kpp.c.eps)
                        
                        w_x(i:end) = w_x_bnd; 
                        break
                    end
                end
            end
        end
        
        function [ Ri, db, du_sqr, Vt_sqr ] = calBulkRichardsonNumber(kpp, grid, B_f, b, u, v)
            db = b(1) - b;
            du_sqr = (u(1) - u).^2 + (v(1) - v).^2;
            Vt_sqr = grid.sop.T_interp_W * calUnresolvedShear(kpp, grid, b, B_f);

            Ri = (- grid.z_T) .* db ./ ( du_sqr + Vt_sqr );
        end
        
        % LMD94 equation (27)
        % Input on T-grid, output on W-grid
        function [ Ri_g ] = calGradientRichardsonNumber(kpp, grid, b, u, v)
            N_sqr = grid.sop.W_ddz_T * b;
            dudz = grid.sop.W_ddz_T * u;
            dvdz = grid.sop.W_ddz_T * v;
            gradU_sqr = dudz.^2 + dvdz.^2;
            Ri_g = grid.W_imask_W * ( N_sqr ./ gradU_sqr );
            Ri_g(gradU_sqr == 0) = 1e20; % or inifinity
        end
        
        % k is the index of deepest layer of mixed-layer 
        % h = grid.h_W(k+1);
        function [ h, k ] = calMixedLayerDepth(kpp, grid, Ri, u_star, L_star, f)
            
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
            
            
            % When stable forcing (B_f < 0 or L_star > 0)
            % h cannot exceeds L_star or Ekman depth as in
            % LMD94 equation (24)
            if (L_star > 0)
                if ( f == 0 )
                    h_E = Inf;
                else
                    h_E = (0.7 * u_star / abs(f));
                end
                
                h_max = min(h_E, L_star);
                
            
            
                if (h > h_max)
                    for i=2:grid.Nz  % start from 2 because h is at least one level
                        if grid.d_W(i+1) > h_max
                            k = i-1;
                            h = grid.d_W(k+1);
                            break;
                        end
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
            
        % Calculate the shear related interior diffusivity given by
        % LMD equation (28). Ri_g is given in equation (27)
        function K_sh = calInteriorK_sh(kpp, grid, b, u, v)
            Ri_g = calGradientRichardsonNumber(kpp, grid, b, u, v);
            K_sh = 50e-4 * kpp.shapeInterior(Ri_g);
        end
        
        function [ K_x_ML, K_x_INT ]  = calK_x(kpp, x, grid, h_k, tau0, B_f, b, u, v)
            d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));
            h = grid.d_W(h_k + 1);
           
            [ w_x, sigma, ~, ~ ] = kpp.calw_x(x, grid.z_W, h, tau0, B_f);
            
            G = kpp.c.calG(sigma);
            W_ML_mask_W = d0( sigma <= 1 ); % Mixed-layer
            W_INT_mask_W = d0( sigma > 1 ); % Interior
            
            K_x_ML = W_ML_mask_W * (h * G .* w_x);     
            K_x_INT = W_INT_mask_W * kpp.calInteriorK_sh(grid, b, u, v);
            
        end
        
        % LMD94 equation (19) and (20) for scalar
        % Also the same as RAD18 equation (19)
        function flux = calNonLocalFlux_s(kpp, grid, h_k, B_f, ws_0)
            
            if (B_f > 0) % unstable case, nonlocal flux is nonzero
                h = grid.d_W(h_k+1);
                sig = grid.d_W / h;
                flux = kpp.c.C_s * kpp.c.calG(sig) * ws_0;
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