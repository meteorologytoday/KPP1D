#if !(:MatrixOperators in names(Main))
#    include(normpath(joinpath(dirname(@__FILE__), "MatrixOperators.jl")))
#end


module KPP

    using SparseArrays
    using Formatting


    include("Grid.jl")
    include("BasicMatrixOperators.jl")
    include("AdvancedMatrixOperators.jl")

    include("MatrixOperators.jl")
    include("constants.jl")
    include("KPP_constants.jl")
    include("KPP_functions.jl")


    include("Env.jl") 
    include("State.jl") 
    include("Core.jl") 
    include("Forcing.jl") 
    include("Model.jl") 

end
