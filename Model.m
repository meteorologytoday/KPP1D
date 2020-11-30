classdef Model < handle 
    properties
        grid
        sop
        H
        Nz
        Kv_iso
        Kv_cva = 1e5
        dt
        state
        kpp
    end
    methods

        function m = Model(H, Nz, Kv_iso, dt)

            grid = makeGrid(H, Nz);
            sop = makeSpatialOperators(grid);

            m.grid = grid;
            m.sop = sop;
            m.H = H;
            m.Nz = Nz;
            m.Kv_iso = Kv_iso;
            m.dt = dt;
            m.state = State(Nz);
            m.kpp = KPPConstants();
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
            
            % update buoyancy
            m.update_b();
            
            % Backward Euler diffusion
            %
            % x_t+1 - x_t = dt * OP * x_t+1
            % (I - dt * OP) x_t+1 = x_t
            % x_t+1 = (I - dt*OP) \ x_t
            op_diffz = mkOp_diffz(m.grid, m.sop, m.state.b, m.Kv_iso,  m.Kv_cva);
            m.state.S = ( (m.grid.T_I_T - m.dt * op_diffz) \ m.state.S );
            m.state.T = (m.grid.T_I_T - m.dt * op_diffz) \ m.state.T;
            
            m.update_b();
            
        end
        
        function update_b(m)
            m.state.b = TS2b(m.state.T, m.state.S);
        end
        
    end
end