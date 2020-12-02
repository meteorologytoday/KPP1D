classdef State < handle
    properties
        u
        v
        T
        S
        b
        Ri
        wb_sfc = 0
        wT_sfc = 0
        wS_sfc = 0
        U      = 0  % surface wind speed on x direction
        V      = 0  % surface wind speed on y direction
        I0     = 0  % radiation           positive upward
        Hf_sen = 0  % sensible heat flux  positive upward
        Hf_lat = 0  % latent heat flux    positive upward
        h      = 0
        h_k    = 0
        Precip = 0  % m / s
        Evap   = 0 % m / s
        tau0   = 0
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