using PolyhedralOmegaMono.Polynomials

function compute_rational_function(cone::Cone{T}, parallelepipeds, x::Union{Vector}) where {T<:Number}
    num = sum([prod(x .^ Int32.(p)) for p in parallelepipeds])
    den = prod([(1 - prod(x .^ Int32.(ray.direction))) for ray in cone.rays])
    return CombinationOfRationalFunctions(Pair(num * (-1)^(!cone.sign), den))
end


function compute_rational_function_str(cone::Cone{T}, parallelepipeds; counting::Bool=false) where {T<:Number}
    function create_monomial_str_from_exponents(z::Vector{<:Number})
        tmp = ""
        for i in eachindex(z)
            if z[i] == 0
                t = "1"
            else
                if counting
                    t = "q^($(Int32(z[i])))"
                else
                    t = "x$(i)^($(Int32(z[i])))"
                end
            end
            if length(tmp) > 0
                if tmp == "1"
                    tmp = t
                elseif t != "1"
                    tmp = tmp * " * " * t
                end
            else
                tmp = t
            end
        end
        return "(" * tmp * ")"
    end
    num = ""
    for p in parallelepipeds
        if (length(num) > 0)
            num = num * " + " * create_monomial_str_from_exponents(p)
        else
            num = create_monomial_str_from_exponents(p)
        end
    end
    den = ""
    for ray in cone.rays
        if (length(den) > 0)
            den = den * " * " * "(1 - $(create_monomial_str_from_exponents(ray.direction)))"
        else
            den = "(1 - $(create_monomial_str_from_exponents(ray.direction)))"
        end
    end
    res = "(" * num * ")/(" * den * ")"
    if (!cone.sign)
        res = "-(" * res * ")"
    end
    return res
end
