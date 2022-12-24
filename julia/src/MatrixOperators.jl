if !(:Grid in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "PolelikeCoordinate.jl")))
end

module MatrixOperators

    using ..PolelikeCoordinate
    
    export BasicMatrixOperators, AdvancedMatrixOperators



end
