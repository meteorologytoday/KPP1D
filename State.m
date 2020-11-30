classdef State
    properties
        T
        S
        b
    end
    methods

        function s = State(Nz)
            s.T = zeros(Nz,1);
            s.S = s.T*0;
            s.b = s.T*0;
        end
    end
end