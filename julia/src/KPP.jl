#if !(:MatrixOperators in names(Main))
#    include(normpath(joinpath(dirname(@__FILE__), "MatrixOperators.jl")))
#end


module KPP

    using SparseArrays
    using Formatting

    include("constants.jl")

    if BUOYANCY_TYPE == :LINEAR
        include("Buoyancy_linear.jl")
        const Buoyancy = Buoyancy_linear
    elseif BUOYANCY_TYPE == :NONLINEAR
        include("Buoyancy_nonlinear.jl")
        const Buoyancy = Buoyancy_nonlinear
    end

    include("Grid.jl")
    include("BasicMatrixOperators.jl")
    include("AdvancedMatrixOperators.jl")

    include("MatrixOperators.jl")


    include("Env.jl") 
    include("State.jl") 
    include("Core.jl") 
    include("Forcing.jl") 
    include("Model.jl") 

    include("KPP_constants.jl")
    include("KPP_functions.jl")
    include("other_functions.jl")

end
