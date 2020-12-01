classdef State < handle
    properties
        u
        v
        T
        S
        b
        Ri
        wb_sfc
        I0
        h
        h_k
    end
    methods

        function s = State(Nz)
            s.T = zeros(Nz,1);
            s.S = s.T * 0;
            s.b = s.T * 0;
            s.Ri = s.T * 0;
            s.u = s.T * 0;
            s.v = s.T * 0;
            
            s.wb_sfc = 0.0;
            s.I0 = 0.0;
        end
        
    end
end