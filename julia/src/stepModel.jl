
function stepModel!(
    m :: Model,
    dt :: Float64,
)

    if m.fo.surf_flux_calculation_type == :ATM_FORCING
                
        # Need to input: U10, V10, T_a, q_a, T_o, I_0
                
        # Diagnose wu, wT, wq fluxes using the method
        # "calSurfaceFluxes" in the class "SurfaceFlux" which is based
        # on Monin-Obhukov similarity theory estimation. See that method
        #  for more paper details.
        wu, wv, wT_sen, wT_lw, wq = KPP.calSurfaceFluxesUsingAtmForcing(
            m.st.U10,
            m.st.V10,
            m.st.T_a,
            m.st.q_a,
            m.st.T[1],
        )

        m.st.τx0    = - wu * ρ_a
        m.st.τy0    = - wv * ρ_a
        m.st.Hf_sen = wT_sen * cp_ρ_a
        m.st.Hf_lw  = wT_lw
        m.st.Hf_lat = wq * m.c.Lv_w
        m.st.evap = wq / m.c.ρ_fw      # kg / m^2 / s / (kg/m^3)
    
    elseif m.fo.surf_flux_calculation_type == :NOCALCULATION
        # do nothing
    else
        throw(ErrorException("Unknown surf_flux_calculation_type: $(m.st.surf_flux_calculation_type)"))
    end
            
    T_0 = m.st.T[1]
    S_0 = m.st.S[1]
    α_T_0 = Buoyancy.TS2α_T(T_0, S_0)
    α_S_0 = Buoyancy.TS2α_S(T_0, S_0)
    
    # The flux calculation needs the following variables specified:
    # taux0, tauy0, Hf_sen, Hf_lat, Hf_lw, I_0, albedo, h_k,
    # precip, evap
    m.fo.τ0   = sqrt(m.fo.τx0^2 + m.fo.τy0^2)
    m.fo.wT_0 = (m.fo.Hf_sen + m.fo.Hf_lat + m.fo.Hf_lw) / cp_ρ_sw
    m.fo.wT_R = m.co.rad.coe_turbulent_flux_T[m.st.h_k] * (1.0 - m.fo.albedo) * m.fo.I_0 / cp_ρ_sw
    m.fo.wS_0 = (m.fo.precip - m.fo.evap) * S_0
    
    m.fo.wb_R = g * α_T_0 * m.fo.wT_R;
    m.fo.wb_0 = g * ( α_T_0 * m.fo.wT_0 - α_S_0 * m.fo.wS_0 )
    m.fo.wu_0 = - m.fo.τx0 / ρ_sw
    m.fo.wv_0 = - m.fo.τy0 / ρ_sw
    m.fo.B_f = m.fo.wb_R + m.fo.wb_0
    
    #stepModel_momentum(m)
    stepModel_radiation!(m, dt) 
    #m.stepModel_momentumNudging();
    #m.stepModel_scalarNudging();
    #m.stepModel_radiation();
    
    #diag_kpp = stepModel_KPP!(m, dt)
    updateBuoyancy!(m)
end

function stepModel_radiation!(
    m :: Model,
    dt :: Float64,
)
    _, Q = calRadiation(m.co.rad, (1.0 - m.fo.albedo) * m.fo.I_0 / cp_ρ_sw )
    #println(Q * dt)
    @. m.st.T += dt * Q
end

function stepModel_KPP!(
    m :: Model,
    dt :: Float64,
)
    
    amo = m.co.amo
    bmo = amo.bmo
    gd  = m.ev.gd
    st  = m.st
    # 1. Determine mixed-layer depth
    update_ML!(m)
    
    # 2. Calculate incoming surface flux
    sfc_flux_T    = zeros(Float64, bmo.W_pts, 1)
    sfc_flux_T[1] = m.fo.wT_0
    
    sfc_flux_S    = zeros(Float64, bmo.W_pts, 1)
    sfc_flux_S[1] = m.fo.wS_0
    
    sfc_flux_u    = zeros(Float64, bmo.W_pts, 1)
    sfc_flux_v    = zeros(Float64, bmo.W_pts, 1)
    sfc_flux_u[1] = m.fo.wu_0
    sfc_flux_v[1] = m.fo.wv_0
    
    # 3. Calculate nonlocal transport
    nloc_flux_T = calNonLocalFlux_s(m.st.h_k, m.fo.wb_0, m.fo.wT_0 + m.fo.wT_R, amo, gd)
    nloc_flux_S = calNonLocalFlux_s(m.st.h_k, m.fo.wb_0, m.fo.wS_0, amo, gd)
    
    
    #m.state.T = m.state.T - m.dt * m.grid.sop.T_ddz_W * nloc_flux_T;
    #m.state.S = m.state.S - m.dt * m.grid.sop.T_ddz_W * nloc_flux_S;
    #m.update_b();


    # 4. Build diffusion operator for Euler backward method
    # calculate K_s of temperature and salinity
    K_s_ML, K_s_INT = calK_x(
        :SCALAR, 
        m.st.h_k,
        m.fo.τ0,
        m.fo.wb_0,
        m.st.b,
        m.st.u,
        m.st.v,
        amo,
        gd,
    )

    # calculate K_m for u, v
    K_m_ML, K_m_INT = calK_x(
        :MOMENTUM,
        m.st.h_k,
        m.fo.τ0,
        m.fo.wb_0,
        m.st.b,
        m.st.u,
        m.st.v,
        amo,
        gd,
    )

    # Total diffusivity is the mixed-layer one plus the internal one
    K_s_total = K_s_ML + K_s_INT
    K_m_total = K_m_ML + K_m_INT
    
    op_local_flux_s = - d0(K_s_total) * amo.W_∂z_T
    op_diffusion_s  = - amo.T_DIVz_W * op_local_flux_s
    M_s = bmo.T_I_T - dt * op_diffusion_s
    
    op_local_flux_m = - d0(K_m_total) * amo.W_∂z_T
    op_diffusion_m = - amo.T_DIVz_W * op_local_flux_m
    
    
    M_m = bmo.T_I_T - dt * op_diffusion_m
    
    loc_flux_S = op_local_flux_s * m.st.S
    loc_flux_T = op_local_flux_s * m.st.T
    
    # 5. Step the model states
    # nonlocal and surface flux
    m.st.S[:] = M_s \ ( m.st.S - dt * amo.T_DIVz_W * ( nloc_flux_S + sfc_flux_S ))
    m.st.T[:] = M_s \ ( m.st.T - dt * amo.T_DIVz_W * ( nloc_flux_T + sfc_flux_T ))

    m.st.u[:] = M_m \ ( m.st.u - dt * amo.T_DIVz_W * sfc_flux_u )
    m.st.v[:] = M_m \ ( m.st.v - dt * amo.T_DIVz_W * sfc_flux_v )
    
    updateBuoyancy!(m)
 
    
    diag_kpp = (
        loc_flux_T  = loc_flux_T,
        loc_flux_S  = loc_flux_S,
        nloc_flux_T = nloc_flux_T,
        nloc_flux_S = nloc_flux_S,
        K_s_ML      = K_s_ML,
        K_s_INT     = K_s_INT,
        K_m_ML      = K_m_ML,
        K_m_INT     = K_m_INT,
    )

    return diag_kpp
end

function update_ML!(
    m :: Model,
)
    u_star, L_star = KPP.calMOSTscales(m.fo.τ0, m.fo.B_f)
    m.st.Ri, db, du_sqr, Vt_sqr = KPP.calBulkRichardsonNumber(m.fo.wb_0, m.st.b, m.st.u, m.st.v, m.co.amo, m.ev.gd)
    
    m.st.h, m.st.h_k = KPP.calMixedLayerDepth(m.st.Ri, u_star, L_star, m.ev.f, m.co.amo, m.ev.gd)

    return m.st.h, m.st.Ri, db, du_sqr, Vt_sqr
end


#=
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

=#
