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
            T_0 = m.state.T(1);
            S_0 = m.state.S(1);
            m.state.tau0 = sqrt(m.state.taux0^2 + m.state.tauy0^2);
            m.state.wb_0 = m.c.g * ( TS2alpha(T_0, S_0) * m.state.wT_0 - TS2beta(T_0, S_0) * m.state.wS_0 );
            m.state.wu_0 = - m.state.taux0 / m.c.rho_sw;
            m.state.wv_0 = - m.state.tauy0 / m.c.rho_sw;
            
            % Step T and S
            %m.stepModel_radiation();

            
            m.stepModel_momentum();
%            m.stepModel_surface_hydrology();

%            m.stepModel_surface_heat_flux();
            diag_kpp = m.stepModel_KPP();
        end
        
        function stepModel_radiation(m)
            [ ~, Q ] = m.rad.calRadiation(m.state.I0 / (m.c.cp_sw * m.c.rho_sw));
            m.state.T = m.state.T + m.dt * Q;
            m.update_b();
        end

        % For Coriolis terms. Use Runge-Kutta 4th order instead of
        % Adam-Bashforth 3rd scheme to avoid storing extra values.
        function stepModel_momentum(m)
            
            [ m.state.u(:), m.state.v(:) ] = m.calRK4Coriolis(m.f, m.dt, m.state.u, m.state.v);

        end
        
        function stepModel_surface_heat_flux(m)
            flux = zeros(m.grid.W_pts, 1);
            flux(1) = m.state.wT_0;
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

            % 1. Determine mixed-layer depth
            m.update_ML();
            
            % 2. Calculate incoming surface flux
            sfc_flux_T = zeros(m.grid.W_pts, 1);
            sfc_flux_T(1) = m.state.wT_0;
            
            sfc_flux_S = zeros(m.grid.W_pts, 1);
            sfc_flux_S(1) = m.state.wS_0;
            
            sfc_flux_u    = zeros(m.grid.W_pts, 1);
            sfc_flux_v    = zeros(m.grid.W_pts, 1);
            sfc_flux_u(1) = m.state.wu_0;
            sfc_flux_v(1) = m.state.wv_0;
            
            % 3. Calculate nonlocal transport
            nloc_flux_T = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_0, m.state.wT_0);
            nloc_flux_S = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_0, m.state.wS_0);
            
            diag_kpp.nloc_flux_T = nloc_flux_T;
            diag_kpp.nloc_flux_S = nloc_flux_S;
    
            %m.state.T = m.state.T - m.dt * m.grid.sop.T_ddz_W * nloc_flux_T;
            %m.state.S = m.state.S - m.dt * m.grid.sop.T_ddz_W * nloc_flux_S;
            %m.update_b();


            % 4. Build diffusion operator for Euler backward method
            % calculate K_s for temperature and salinity
            [ K_s_ML, K_s_INT ]  = m.kpp.calK_x(m.kpp.c.SCALAR, m.grid, m.state.h_k, m.state.tau0, m.state.wb_0, m.state.b, m.state.u, m.state.v);
            % calculate K_m for u, v
            [ K_m_ML, K_m_INT ]  = m.kpp.calK_x(m.kpp.c.MOMENTUM, m.grid, m.state.h_k, m.state.tau0, m.state.wb_0, m.state.b, m.state.u, m.state.v);
            %K_cva = m.calConvectiveAdjustmentK(m.grid, m.state.b, m.Kv_cva);

            K_s_total = K_s_ML + K_s_INT ;%+ K_cva;
            diag_kpp.K_s_ML = K_s_ML;
            diag_kpp.K_s_INT = K_s_INT;
            
            K_m_total = K_m_ML + K_m_INT;
            
            op_diffusion_s = m.grid.sop.T_ddz_W * d0(K_s_total) * m.grid.sop.W_ddz_T;
            op_diffusion_m = m.grid.sop.T_ddz_W * d0(K_m_total) * m.grid.sop.W_ddz_T;
            
            M_s = m.grid.T_I_T - m.dt * op_diffusion_s;
            M_m = m.grid.T_I_T - m.dt * op_diffusion_m;
            
            % 5. Step the model states
            % nonlocal and surface flux
            m.state.S = M_s \ ( m.state.S - m.dt * m.grid.sop.T_ddz_W * ( nloc_flux_S + sfc_flux_S));
            m.state.T = M_s \ ( m.state.T - m.dt * m.grid.sop.T_ddz_W * ( nloc_flux_T + sfc_flux_T));
            m.state.u = M_m \ ( m.state.u - m.dt * m.grid.sop.T_ddz_W * sfc_flux_u);
            m.state.v = M_m \ ( m.state.v - m.dt * m.grid.sop.T_ddz_W * sfc_flux_v);
            
            m.update_b();
                    
        end
        
        function update_b(m)
            m.state.b = TS2b(m.state.T, m.state.S);
        end
        
        function [ h, Ri, db, du_sqr, Vt_sqr ] = update_ML(m)
            [ u_star, L_star ] = m.kpp.calMOSTscales(m.state.tau0, m.state.wb_0);
            [Ri, db, du_sqr, Vt_sqr ]  = m.kpp.calBulkRichardsonNumber(m.grid, m.state.wb_0, m.state.b, m.state.u, m.state.v);
            
            m.state.Ri(:) = Ri;
            [m.state.h, m.state.h_k] = m.kpp.calMixedLayerDepth(m.grid, m.state.Ri, u_star, L_star, m.f);
            
            h = m.state.h;
        end
        
        function [ wT_0, wS_0, wb_0, tau0 ] = calSurfaceForcing(m)
            
            %, I0, U, T_air, T_sea, RH, P, E)
            % calculate wT_0
            
            % calculate wS_0
            
            % calculate wb_0 through wT_0 and wS_0
            
            % calculate tau0
            
        end
        
        function [ K ] = calConvectiveAdjustmentK(m, grid, b, K_cva)
            db = (grid.W_DN_T - grid.W_UP_T) * b;
            K = zeros(grid.W_pts, 1);
            K( db < 0 ) = K_cva;
        end
        
        function [ G_u, G_v ] = calCoriolisTerms(m, f, u, v)
            G_u = f * v;
            G_v = - f * u;
        end
        
        function [u_new, v_new] = calRK4Coriolis(m, dt, f, u, v)
            
            [ G_u_1, G_v_1 ] = m.calCoriolisTerms(f, u, v);
            
            u_tmp = u + dt/2 * G_u_1;
            v_tmp = v + dt/2 * G_v_1;
            [ G_u_2, G_v_2 ] = m.calCoriolisTerms(f, u_tmp, v_tmp);
            
            u_tmp = u + dt/2 * G_u_2;
            v_tmp = v + dt/2 * G_v_2;
            [ G_u_3, G_v_3 ] = m.calCoriolisTerms(f, u_tmp, v_tmp);
            
            u_tmp = u + dt * G_u_3;
            v_tmp = v + dt * G_v_3;
            [ G_u_4, G_v_4 ] = m.calCoriolisTerms(f, u_tmp, v_tmp);
            
            u_new = u + dt/6 * (G_u_1 + 2*G_u_2 + 2*G_u_3 + G_u_4);
            v_new = v + dt/6 * (G_v_1 + 2*G_v_2 + 2*G_v_3 + G_v_4);
        end
    end
end