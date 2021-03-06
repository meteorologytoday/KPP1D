classdef Model < handle 
    properties
        grid
        c
        sf
        dt
        state
        kpp
        rad
        f
        momentum_nudging_timescale = 30 * 60
        scalar_nudging_timescale = 30 * 60
        SURFFLUX_REALISTIC = 1
        SURFFLUX_SIMPLE = 2
    end
    methods

        function m = Model(z_W, dt, f)

            grid = makeGrid(z_W);
            c = Const();
            
            m.c = c;
            m.sf = SurfaceFlux();
            m.grid = grid;
            m.dt = dt;
            m.state = State(grid.Nz);
            m.kpp = KPP();
            m.rad = Radiation(grid);
            m.f = f;
            
            m.state.h_k = 1;
            m.state.h = m.grid.d_W(m.state.h_k + 1);
        end
        
        function showModelInfo(m)
            fprintf("===== Model Info =====\n");
            fprintf("Nz= %d\n", m.grid.Nz);
            fprintf("dt = %f\n", m.dt);
            fprintf("f = %f\n", m.f);
        end
        
        function [ diag_kpp ] = stepModel(m, surf_flux_calculation_type)
            
            if (surf_flux_calculation_type == m.SURFFLUX_REALISTIC)
                
                % Need to input: U10, V10, T_a, q_a, T_o, I_0
                
                % Diagnose wu, wT, wq fluxes using the method
                % "calSurfaceFluxes" in the class "SurfaceFlux" which is based
                % on Monin-Obhukov similarity theory estimation. See that method
                % for more paper details.
                [ wu, wv, wT_sen, wT_lw, wq ] = m.sf.calSurfaceFluxes(m.state.U10, m.state.V10, m.state.T_a, m.state.q_a, m.state.T(1));

                m.state.taux0  = - wu * m.c.rho_a;
                m.state.tauy0  = - wv * m.c.rho_a;
                m.state.Hf_sen = wT_sen * m.c.cprho_a;
                m.state.Hf_lw  = wT_lw;
                m.state.Hf_lat = wq * m.c.Lv_w;
                m.state.evap = wq / m.c.rho_fw;   % kg / m^2 / s / (kg/m^3)

            elseif (surf_flux_calculation_type == m.SURFFLUX_SIMPLE)
                %m.state.taux0  = 0.01;
                % do nothing
            else
                error('Unknown surf_flux_calculation_type: %d', surf_flux_calculation_type);
            end
            
            
            T_0 = m.state.T(1);
            S_0 = m.state.S(1);
            alpha_0 = TS2alpha(T_0, S_0);
            beta_0 = TS2beta(T_0, S_0);
            
            % The flux calculation needs the following variables specified:
            % taux0, tauy0, Hf_sen, Hf_lat, Hf_lw, I_0, albedo, h_k,
            % precip, evap
            m.state.tau0 = sqrt(m.state.taux0^2 + m.state.tauy0^2);
            m.state.wT_0 = (m.state.Hf_sen + m.state.Hf_lat + m.state.Hf_lw) / m.c.cprho_sw;
            m.state.wT_R = m.rad.coe_turbulent_flux_T(m.state.h_k) * (1 - m.state.albedo) * m.state.I_0 / m.c.cprho_sw;
            m.state.wS_0 = (m.state.precip - m.state.evap) * S_0;
            
            m.state.wb_R = m.c.g * alpha_0 * m.state.wT_R;
            m.state.wb_0 = m.c.g * ( alpha_0 * m.state.wT_0 - beta_0 * m.state.wS_0 );
            m.state.wu_0 = - m.state.taux0 / m.c.rho_sw;
            m.state.wv_0 = - m.state.tauy0 / m.c.rho_sw;
            m.state.B_f = m.state.wb_R + m.state.wb_0;
            
            
            m.stepModel_momentum();
            
            m.stepModel_momentumNudging();
            m.stepModel_scalarNudging();
            m.stepModel_radiation();
            
            diag_kpp = m.stepModel_KPP();
        end
                
        function stepModel_radiation(m)
            [ ~, Q ] = m.rad.calRadiation( (1 - m.state.albedo) * m.state.I_0 / m.c.cprho_sw);
            m.state.T = m.state.T + m.dt * Q;
        end

        % For Coriolis terms. Use Runge-Kutta 4th order instead of
        % Adam-Bashforth 3rd scheme to avoid storing extra values.
        function stepModel_momentum(m)
            
            [ m.state.u(:), m.state.v(:) ] = m.calRK4Coriolis(m.f, m.dt, m.state.u, m.state.v);
            
        end
 
        function stepModel_momentumNudging(m)
            
            nudge_idx = isfinite(m.state.u_nudging);
            % (u_t+1 - u_t) / dt = - 1 / tau * ( u_t+1 - u_nudging )
            % u_t+1 - u_t = - dt / tau * ( u_t+1 - u_nudging )
            % u_t+1 (1 + dt / tau) = u_t + (dt / tau) * u_nudging
            % u_t+1 = 1 / (1+dt/tau) * u_t + (dt/tau) * u_nudging
            r = m.dt / m.momentum_nudging_timescale;
            m.state.u(nudge_idx) = (m.state.u(nudge_idx) + r * m.state.u_nudging(nudge_idx)) / (1+r);
            m.state.v(nudge_idx) = (m.state.v(nudge_idx) + r * m.state.v_nudging(nudge_idx)) / (1+r);

        end
        

        
        function stepModel_scalarNudging(m)
            % x_t+1 = (x_t + r * X) / (1+r)
            % (x_t+1 - x_t) / dt = (x_t + r*X - x_t - r * x_t)/( (1+r)*dt )
            % fix_rate_(t+0.5) = (r*X - r*x_t) / ( ( 1 + r ) * dt )
            
            [ m.state.T(:), m.state.T_nudging_fix(:) ] = nudgingProfileHelper(m.state.T, m.state.T_nudging, m.dt, m.scalar_nudging_timescale);
            [ m.state.S(:), m.state.S_nudging_fix ] = nudgingProfileHelper(m.state.S, m.state.S_nudging, m.dt, m.scalar_nudging_timescale);
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
            nloc_flux_T = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_0, m.state.wT_0 + m.state.wT_R);
            nloc_flux_S = m.kpp.calNonLocalFlux_s(m.grid, m.state.h_k, m.state.wb_0, m.state.wS_0);
            
            
            %m.state.T = m.state.T - m.dt * m.grid.sop.T_ddz_W * nloc_flux_T;
            %m.state.S = m.state.S - m.dt * m.grid.sop.T_ddz_W * nloc_flux_S;
            %m.update_b();


            % 4. Build diffusion operator for Euler backward method
            % calculate K_s for temperature and salinity
            [ K_s_ML, K_s_INT ]  = m.kpp.calK_x(m.kpp.c.SCALAR, m.grid, m.state.h_k, m.state.tau0, m.state.wb_0, m.state.b, m.state.u, m.state.v);

            % calculate K_m for u, v
            [ K_m_ML, K_m_INT ]  = m.kpp.calK_x(m.kpp.c.MOMENTUM, m.grid, m.state.h_k, m.state.tau0, m.state.wb_0, m.state.b, m.state.u, m.state.v);

            K_s_total = K_s_ML + K_s_INT ;
            K_m_total = K_m_ML + K_m_INT;
            
            op_local_flux_s = - d0(K_s_total) * m.grid.sop.W_ddz_T;
            op_diffusion_s = - m.grid.sop.T_ddz_W * op_local_flux_s;
            M_s = m.grid.T_I_T - m.dt * op_diffusion_s;
            
            op_local_flux_m = - d0(K_m_total) * m.grid.sop.W_ddz_T;
            op_diffusion_m = - m.grid.sop.T_ddz_W * op_local_flux_m;
            
            
            M_m = m.grid.T_I_T - m.dt * op_diffusion_m;
            
            diag_kpp.loc_flux_S = op_local_flux_s * m.state.S;
            diag_kpp.loc_flux_T = op_local_flux_s * m.state.T;
            
            % 5. Step the model states
            % nonlocal and surface flux
            m.state.S = M_s \ ( m.state.S - m.dt * m.grid.sop.T_ddz_W * ( nloc_flux_S + sfc_flux_S ));
            m.state.T = M_s \ ( m.state.T - m.dt * m.grid.sop.T_ddz_W * ( nloc_flux_T + sfc_flux_T ));
    
            m.state.u = M_m \ ( m.state.u - m.dt * m.grid.sop.T_ddz_W * sfc_flux_u);
            m.state.v = M_m \ ( m.state.v - m.dt * m.grid.sop.T_ddz_W * sfc_flux_v);
            
            m.update_b();

            
            diag_kpp.nloc_flux_T = nloc_flux_T;
            diag_kpp.nloc_flux_S = nloc_flux_S;
            diag_kpp.K_s_ML = K_s_ML;
            diag_kpp.K_s_INT = K_s_INT;
            diag_kpp.K_m_ML = K_m_ML;
            diag_kpp.K_m_INT = K_m_INT;
            
        end
        
        function update_b(m)
            m.state.b = TS2b(m.state.T, m.state.S);
        end
        
        function [ h, Ri, db, du_sqr, Vt_sqr ] = update_ML(m)
            [ u_star, L_star ] = m.kpp.calMOSTscales(m.state.tau0, m.state.B_f);
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

function [ nudged_profile, fix_rate ] = nudgingProfileHelper(profile, target_profile, dt, nudging_timescale)
            % x_t+1 = (x_t + r * X) / (1+r)
            % (x_t+1 - x_t) / dt = (x_t + r*X - x_t - r * x_t)/( (1+r)*dt )
            % fix_rate_(t+0.5) = (r*X - r*x_t) / ( ( 1 + r ) * dt )
            
            nudged_profile = profile * 1;
            fix_rate = profile * 0;

            r = dt / nudging_timescale;
            
            nudge_idx = isfinite(target_profile);
            nudged = (profile(nudge_idx) + r * target_profile(nudge_idx)) / (1+r);
            
            fix_rate(nudge_idx) = (nudged - profile(nudge_idx)) / dt;
            nudged_profile(nudge_idx) = nudged;
end