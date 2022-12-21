if !(:PolelikeCoordinate in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "PolelikeCoordinate.jl")))
end

if !(:MatrixOperators in names(Main))
    include(normpath(joinpath(dirname(@__FILE__), "MatrixOperators.jl")))
end


module KPP

    using SparseArrays
    using Formatting
    using ..PolelikeCoordinate
    using ..MatrixOperators: BasicMatrixOperators, AdvancedMatrixOperators


    include("constants.jl")
    include("KPP_constants.jl")
    include("KPP_functions.jl")




end
