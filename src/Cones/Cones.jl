module Cones

using PolyhedralOmegaMono.Polynomials
include("Ray.jl")
include("Cone.jl")
include("CombinationOfCones.jl")
include("EnumerateFundamentalParallelepiped.jl")
include("ComputeRationalFunction.jl")
include("Utils.jl")

export
    # Types
    Ray,
    Cone,
    CombinationOfCones,
    # Functions
    vrep_matrix,
    flip,
    isforward,
    primitive,
    substitute_ray,
    substitute_cone,
    enumerate_fundamental_parallelepiped,
    compute_rational_function,
    compute_rational_function_str

end
