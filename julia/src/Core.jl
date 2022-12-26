mutable struct Core

    amo :: AdvancedMatrixOperators
    rad :: Radiation

    function Core(
        ev :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd, mask_T=ev.mask_T)
        rad = Radiation(amo=amo)

        return new(
            amo,
            rad,
        )    
    end
end
