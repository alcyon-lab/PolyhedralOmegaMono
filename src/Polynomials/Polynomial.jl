using Symbolics

const Polynomial = Num

function create_monomial_from_exponents(z::Vector{<:Number}, x::Union{Vector})::Polynomial
    assert(length(z) <= length(x))
    tmp = 1
    for i in eachindex(z)
        tmp = tmp * (x[i]^Int(z[i]))
    end
    return tmp
end
