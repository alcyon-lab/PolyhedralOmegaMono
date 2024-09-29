module Polynomials
using Symbolics
include("Polynomial.jl")
include("CombinationOfRationalFunctions.jl")

export
    # Types
    CombinationOfRationalFunctions, Polynomial,
    # Functions
    get_ratfun,
    create_monomial_from_exponents,
    get_variable_array,
    evaluate,
    evaluate_all_with,
    simplify
end
