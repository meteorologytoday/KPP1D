module Buoyancy_nonlinear
 
    # Follow Bryan and Cox (1972) : An Approximate Equation of State for Numerical Models of Ocean Circulation
    # Here we use Table 3, and arbitrary pick Z=0m as the formula.
   
    const g               = 9.80616     # m / s^2      copied from models/csm_share/shr/shr_const_mod.F90
 
    const T_ref =   13.5
    const S_ref =   32.6
    const ρ_ref = 1024.458

    const ρ1 = -.20134    / ρ_ref
    const ρ2 =  .77096    / ρ_ref
    const ρ3 = -.49261e-2 / ρ_ref
    const ρ4 =  .46092e-3 / ρ_ref
    const ρ5 = -.20105e-2 / ρ_ref
    const ρ6 =  .36597e-4 / ρ_ref
    const ρ7 =  .47371e-5 / ρ_ref
    const ρ8 =  .37735e-4 / ρ_ref
    const ρ9 =  .65493e-5 / ρ_ref
   
    function TS2b(
        T :: Float64,
        S :: Float64,
    )

        dT = T - T_ref
        dS = S - S_ref
        b  = - g * ( ρ1 * dT + ρ2 * dS +  ρ3*(dT^2) + ρ4*(dS^2) 
            + ρ5 * (dT * dS) + ρ6 * (dT^3) + ρ7 * (dS^2)*dT + ρ8*(dT^2)*dS + ρ9*(dS^3) 
        )
        return b
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
        dT = T - T_ref
        dS = S - S_ref
        α_T = - ( ρ1 + 2*ρ3*dT + ρ5*dS + 3*ρ6*(dT^2) + ρ7*(dS^2) + 2*ρ8*dS*dT )
        return α_T
    end
    
    # Salinity expansion coefficient := - (dρ/dS)_T  (notice the negative sign)
    function TS2α_S(T :: Float64, S :: Float64)
        dT = T - T_ref
        dS = S - S_ref
        α_S = ρ2 + 2*ρ4*dS + ρ5*dT + 2*ρ7*dS*dT + ρ8*(dT^2) + 3*ρ9*(dS^2)
        return α_S
    end

end
