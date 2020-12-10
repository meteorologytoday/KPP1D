classdef Recorder < handle 
    
    properties
        var2d
        var1d
        Nz
        Nt
        T
        S
        u
        v
        Hf_lat
        Hf_sen
        Hf_lw
        I_0
        T_a
        h
        precip
        evap
        taux0
        tauy0
        tau0
        U10
        V10
        u_nudging_fix
        v_nudging_fix
        T_nudging_fix
        S_nudging_fix
    end
    
    methods

        function rec = Recorder(Nz, Nt)
            rec.var2d = {'T', 'S', 'u', 'v'};
            rec.var1d = {'Hf_sen', 'Hf_lat', 'Hf_lw', 'T_a', 'I_0', 'h', ...
                'precip', 'evap', 'taux0', 'tauy0', 'tau0', 'U10', 'V10', ...
                %'u_nudging_fix', 'v_nudging_fix', 'T_nudging_fix', 'S_nudging_fix', ...
            };
            
            rec.Nz = Nz;
            rec.Nt = Nt;
            
            for i = 1:length(rec.var2d)
                varname = rec.var2d{i};
                rec.(varname) = zeros(Nz, Nt);
            end
            
            for i = 1:length(rec.var1d)
                varname = rec.var1d{i};
                rec.(varname) = zeros(1, Nt);
            end
        end
        
        function [ ] = record(rec, t_idx, m)
            for i = 1:length(rec.var2d)
                varname = rec.var2d{i};
                rec.(varname)(:, t_idx) = m.state.(varname);
            end
            
            for i = 1:length(rec.var1d)
                varname = rec.var1d{i};
                rec.(varname)(t_idx) = m.state.(varname);
            end
            
%            rec.S(:, t_idx) = m.state.S;
%            rec.u(:, t_idx) = m.state.u;
%            rec.v(:, t_idx) = m.state.v;
%            rec.Hf_sen(t_idx) = m.state.Hf_sen;
%            rec.Hf_lat(t_idx) = m.state.Hf_lat;
%            rec.Hf_lw(t_idx) = m.state.Hf_lw;
        end
    end
end