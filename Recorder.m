classdef Recorder < handle 
    
    properties
        Nz
        Nt
        T
        S
        u
        v
    end
    
    methods

        function rec = Recorder(Nz, Nt)
            rec.Nz = Nz;
            rec.Nt = Nt;
            rec.T = zeros(Nz, Nt);
            rec.S = zeros(Nz, Nt);
            rec.u = zeros(Nz, Nt);
            rec.v = zeros(Nz, Nt);
        end
        
        function [ ] = record(rec, t_idx, m)
            rec.T(:, t_idx) = m.state.T;
            rec.S(:, t_idx) = m.state.S;
            rec.u(:, t_idx) = m.state.u;
            rec.v(:, t_idx) = m.state.v;
        end
    end
end