classdef Model < handle 
    properties
        grid
        c
        H
        Nz
        Kv_iso
        Kv_cva = 1e5
        dt
        state
        kpp
        rad
        f
    end
    methods

        function m = Model(H, Nz, Kv_iso, dt, f)

            grid = makeGrid(H, Nz);
            c = Const();

            m.c = c;
            m.grid = grid;
            m.H = H;
            m.Nz = Nz;
            m.Kv_iso = Kv_iso;
            m.dt = dt;
            m.state = State(Nz);
            m.kpp = KPP();
            m.rad = Radiation(grid);
            m.f = f;
        end
        
        function showModelInfo(m)
            fprintf("===== Model Info =====\n");
            fprintf("Nz= %d\n", m.Nz);
            fprintf("H = %f\n", m.H);
            fprintf("K_iso = %f\n", m.Kv_iso);
            fprintf("K_cva = %f\n", m.Kv_cva);
            fprintf("dt = %f\n", m.dt);
            fprintf("f = %f\n", m.f);
        end
        
        function [ diag_kpp ] = stepModel(m)
            % calculate surface forcing for temperature and salinity
            % m.calSurfaceForcing();
            T_sfc = m.state.T(1);
            S_sfc = m.state.S(1);
            m.state.wb_sfc = m.c.g * ( TS2alpha(T_sfc, S_sfc) * m.state.wT_sfc - TS2beta(T_sfc, S_sfc) * m.state.wS_sfc );
            
            %m.stepModel_radiation();

            m.stepModel_surface_heatflux();
%            m.stepModel_surface_hydrology();

           diag_kpp = m.stepModel_KPP();
        end
        
        function stepModel_radiation(m)
            [ ~, Q ] = m.rad.calRadiation(m.state.I0 / (m.c.cp_sw * m.c.rho_sw));
            m.state.T = m.state.T + m.dt * Q;
        end
        
        function stepModel_surface_heatflux(m)
            flux = zeros(m.grid.W_pts, 1);
            flux(1) = m.state.wT_sfc;
            m.state.T = m.state.T - m.dt * m.grid.sop.T_ddz_W * flux;
        end
 
        function stepModel_diffusion(m)
            
            % update buoyancy
            m.update_b();
            
            % Backward Euler diffusion
            %
            % x_t+1 - x_t = dt * OP * x_t+1
            % (I - dt * OP) x_t+1 = x_t
            % x_t+1 = (I - dt*OP) \ x_t
            op_diffz = mkOp_diffz(m.grid, m.state.b, m.Kv_iso,  m.Kv_cva);
            m.state.S = ( (m.grid.T_I_T - m.dt * op_diffz) \ m.state.S );
            m.state.T = (m.grid.T_I_T - m.dt * op_diffz) \ m.state.T;
            
            m.update_b();
            
        end

        function diag_kpp = stepModel_KPP(m)
            
            d0 = @(v) spdiags(v(:),0,length(v(:)),length(v(:)));
            
            m.update_b();

            % determine mixed-layer depth
            m.update_ML();
            

            % step model - nonlocal transport
            nloc_flux_T = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_sfc, m.state.wT_sfc);
            nloc_flux_S = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_sfc, m.state.wS_sfc);
            diag_kpp.nloc_flux_T = nloc_flux_T;
            diag_kpp.nloc_flux_S = nloc_flux_S;
    
            m.state.T = m.state.T - m.dt * m.grid.sop.T_ddz_W * nloc_flux_T;
            m.state.S = m.state.S - m.dt * m.grid.sop.T_ddz_W * nloc_flux_S;
            
            %disp(size(nloc_flux_T));
            %disp(size(m.state.b));
            %disp(size(m.state.T));
            %disp(size(m.state.S));
            %disp('=====still in top=====')
            
            m.update_b();


            
            % step model - diffusion, including convective adjustment
            % calculate K_s for temperature and salinity
            [ K_s_ML, K_s_INT ]  = m.kpp.calK_s(m.grid, m.state.h_k, m.state.tau0, m.state.wb_sfc, m.state.b, m.state.u, m.state.v);
            %K_cva = m.calConvectiveAdjustmentK(m.grid, m.state.b, m.Kv_cva);

            %disp('dimension before add');
            %disp(size(K_s_ML));
            %disp(size(K_s_INT));
            %disp(size(K_cva));

            K_s_total = K_s_ML + K_s_INT ;%+ K_cva;
            diag_kpp.K_s_ML = K_s_ML;
            diag_kpp.K_s_INT = K_s_INT;
            
            %disp('dimension');
            %disp(size(K_s_ML));
            %%disp(size(K_s_INT));
            %disp(size(K_cva));
            %disp(size(K_s_total));
            op_diffusion = m.grid.sop.T_ddz_W * d0(K_s_total) * m.grid.sop.W_ddz_T;
            M = m.grid.T_I_T - m.dt * op_diffusion;
            
            %disp('size of S, T');
            %disp(size(m.state.S));
            %disp(size(m.state.T));
            m.state.S = M \ m.state.S;
            m.state.T = M \ m.state.T;
                        
            %disp('size of S, T 2');
            %disp(size(m.state.S));
            %disp(size(m.state.T));
            
            m.update_b();
                    
        end
        
        function update_b(m)
            m.state.b = TS2b(m.state.T, m.state.S);
        end
        
        function [ h, Ri, db, du_sqr, Vt_sqr ] = update_ML(m)
            [ u_star, ~ ] = m.kpp.calMOSTscales(m.state.tau0, m.state.wb_sfc);
            [Ri, db, du_sqr, Vt_sqr ]  = m.kpp.calBulkRichardsonNumber(m.grid, m.state.wb_sfc, m.state.b, m.state.u, m.state.v);
            
            m.state.Ri(:) = Ri;
            [m.state.h, m.state.h_k] = m.kpp.calMixedLayerDepth(m.grid, m.state.Ri, u_star, m.f);
            
            h = m.state.h;
        end
        
        function [ wT_sfc, wS_sfc, wb_sfc, tau0 ] = calSurfaceForcing(m)
            
            %, I0, U, T_air, T_sea, RH, P, E)
            % calculate wT_sfc
            
            % calculate wS_sfc
            
            % calculate wb_sfc through wT_sfc and wS_sfc
            
            % calculate tau0
            
        end
        
        function [ K ] = calConvectiveAdjustmentK(m, grid, b, K_cva)
            db = (grid.W_DN_T - grid.W_UP_T) * b;
            K = zeros(grid.W_pts, 1);
            K( db < 0 ) = K_cva;
        end
    end
end