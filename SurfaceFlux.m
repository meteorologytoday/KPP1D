classdef SurfaceFlux < handle 
    
    properties
        c
    end
    
    methods

        function sf = SurfaceFlux()
            sf.c = Const();
        end
        
        % Calculation of heat fluxes follows Large & Pond (1982)
        % bulk formula equation (8)
        % Parameterization of coefficient follows Large et al. (1994)
        % equation (A1). 
        function [ wu, wv, wT_sen, wT_R, wq ] = calSurfaceFluxes(sf, u10, v10, T_a, q_a, T_o)
            U10 = sqrt(u10^2 + v10^2);
            stable = T_a > T_o;
            [ C_D, C_T, C_q ] = sf.calCoefficients(U10, stable);

            % unit of q = kg / m^3
            q_o = sf.calSeaSurfaceSpecificHumidity(T_o);
            
            % always assume wind at the ocean surface is zero
            wU = - C_D * U10.^2.0;
            wu = wU * u10/U10; 
            wv = wU * v10/U10; 
            
            wT_sen =   C_T * U10 * (T_o - T_a);
            
            % The longwave radiation. Remember that absolute
            % temperature is used to do the Taylor expansion.
            wT_R = 4 * sf.c.sigma_sb * (T_o + 273.15)^3 * (T_o - T_a);
            
            
            wq =   C_q * U10 * (q_o - q_a);
        end

        % C_T = C_\theta in Large et al. (1994)
        % C_q = C_E      in Large et al. (1994)
        function [ C_D, C_T, C_q ] = calCoefficients(sf, U10, stable)
            
            if U10 == 0
                error('U10 cannot be zero.');
            end
            
            C_D = (2.7 ./ U10 + 0.142 + 0.0764 .* U10) * 1e-3;
            C_D_sqrt = sqrt(C_D);
            
            if (stable)
                C_T = 18.0 * C_D_sqrt * 1e-3;
            else
                C_T = 32.7 * C_D_sqrt * 1e-3;
            end
            
            C_q = 34.6 * C_D_sqrt; 
        end

        function q_o = calSeaSurfaceSpecificHumidity(sf, T)
            
            % Claus-Clapyreon relation.
            % Directly given by Large and Pond (1982)
            % in the paragraph between equation (9) and (10)
            q_o = 0.98 * sf.calSaturatedSpecificHumidity(T);
            
        end
        
        function q_s = calSaturatedSpecificHumidity(sf, T)
            
            % Claus-Clapyreon relation.
            % Directly given by Large and Pond (1982)
            % in the paragraph between equation (9) and (10)
            % Notice that the original formulat gives the unit of g / m^3
            % so we convert it to kg / m^3 by dividing it by 1000.
            q_s =  64038e4 * exp( -5107.4 ./ (T + sf.c.Cel_Kel_offset) ) / 1000;
            
        end
        
    end
end