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
    end
    methods

        function m = Model(H, Nz, Kv_iso, dt)

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
            
        end
        
        function showModelInfo(m)
            fprintf("===== Model Info =====\n");
            fprintf("Nz= %d\n", m.Nz);
            fprintf("H = %f\n", m.H);
            fprintf("K_iso = %f\n", m.Kv_iso);
            fprintf("K_cva = %f\n", m.Kv_cva);
            fprintf("dt = %f\n", m.dt);
        end
        
        function stepModel(m)
            m.stepModel_radiation();
            m.stepModel_diffusion();
        end
        
        function stepModel_radiation(m)
            [ ~, Q ] = m.rad.calRadiation(m.state.I0 / (m.c.cp_sw * m.c.rho_sw));
            m.state.T = m.state.T + m.dt * Q;
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
        
        function update_b(m)
            m.state.b = TS2b(m.state.T, m.state.S);
        end
        
        function [ h, Ri, db, du_sqr, Vt_sqr ] = update_ML(m)
            [Ri, db, du_sqr, Vt_sqr ]  = m.kpp.calBulkRichardsonNumber(m.grid, m.state.wb_sfc, m.state.b, m.state.u, m.state.v);
            m.state.Ri(:) = Ri;
            
            [m.state.h, m.state.h_k] = m.kpp.calMixedLayerDepth(m.grid, m.state.Ri);
            
            h = m.state.h;
        end
    end
end