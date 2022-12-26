module Buoyancy_linear
    
    const T_ref =  298.15 - 273.15
    const S_ref =  35.0
    const α_T   = 2e-4
    const α_S   = 8e-4
    const g     = 9.80616     # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
    
    function TS2b(
        T :: Float64,
        S :: Float64,
    )
        return - g * (1.0 - α_T * (T - T_ref) + α_S * (S - S_ref)) 
    end

    function TS2b!(
        T :: AbstractArray{Float64},
        S :: AbstractArray{Float64},
        b :: AbstractArray{Float64},
    )

        for i = 1:length(T)
            b[i] = TS2b(T[i], S[i])
        end
    end

    
    # Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
    # Here we use Table 3, and arbitrary pick Z=0m as the formula.

    # Thermal expansion coefficient  := (dρ/dT)_S
    function TS2α_T(T::Float64, S::Float64)
        return α_T
    end
    
    # Salinity expansion coefficient := - (dρ/dS)_T  (notice the negative sign)
    function TS2α_S(T :: Float64, S :: Float64)
        return α_S
    end

end
