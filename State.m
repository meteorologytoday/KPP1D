classdef State < handle
    properties
        u
        v
        T
        S
        b
        Ri
        u_nudging
        v_nudging
        B_f = 0
        wT_0 = 0
        wS_0 = 0
        wb_0 = 0
        wT_R = 0
        wb_R = 0
        wu_0 = 0
        wv_0 = 0
        U10    = 0  % surface 10m wind speed on x direction
        V10    = 0  % surface 10m wind speed on y direction
        I_0    = 0  % radiation           positive upward
        Hf_sen = 0  % sensible heat flux  positive upward
        Hf_lat = 0  % latent heat flux    positive upward
        Hf_lw  = 0  % longwave radiation heat flux. This is the radiation from the blackbody radiation rather than solar irradiance I0
        h      = 0
        h_k    = 0
        precip = 0  % m / s
        evap   = 0 % m / s
        tau0   = 0 % Pa
        taux0  = 0 % Pa
        tauy0  = 0 % Pa
        albedo = 0.06;
        T_a = 0;
        q_a = 0;
    end
    methods

        function s = State(Nz)
            s.T = zeros(Nz,1);
            s.S = s.T * 0;
            s.b = s.T * 0;
            s.Ri = s.T * 0;
            s.u = s.T * 0;
            s.v = s.T * 0;
            s.u_nudging = s.T * 0;
            s.v_nudging = s.T * 0;
        end
        
    end
end