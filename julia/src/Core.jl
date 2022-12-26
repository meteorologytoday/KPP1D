mutable struct Core

    amo :: AdvancedMatrixOperators
    rad :: Radiation

    function Core(
        env :: Env,
    )

        amo = AdvancedMatrixOperators(gd=ev.gd, mask=ev.mask_T)
        rad = Radiation(amo=amo)

        return new(
            amo,
            rad,
        )    
    end
end
