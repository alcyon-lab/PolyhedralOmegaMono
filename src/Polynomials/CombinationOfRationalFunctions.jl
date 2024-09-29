struct CombinationOfRationalFunctions
    ratfuns::Vector{Pair{Polynomial,Polynomial}}

    function CombinationOfRationalFunctions(ratfuns::Vector{Pair{Polynomial,Polynomial}})
        new(ratfuns)
    end
    function CombinationOfRationalFunctions(ratfun::Pair{Polynomial,Polynomial})
        new([ratfun])
    end

    function CombinationOfRationalFunctions()
        new(Vector{Pair{Polynomial,Polynomial}}([]))
    end
end

function Base.:(+)(combination::CombinationOfRationalFunctions, ratfun::Pair{Polynomial,Polynomial})
    # TODO: DO in place instead
    ## Why we need to do this in place?
    newpairs = Vector{Pair{Polynomial,Polynomial}}([])
    found = false
    for (numerator, denominator) in combination.ratfuns
        # This is potentially very expensive computationally
        if isequal(denominator, ratfun[2]) # check if this expands. if so, consider alternstive representation.
            push!(newpairs, Pair(numerator + ratfun[1], denominator))
            found = true
        else
            push!(newpairs, Pair(numerator, denominator))
        end
    end
    if !found
        push!(newpairs, ratfun)
        # push!(combination.ratfuns, ratfun)
    end
    return CombinationOfRationalFunctions(newpairs)
end


function Base.:(+)(combination1::CombinationOfRationalFunctions, combination2::CombinationOfRationalFunctions)
    combination = CombinationOfRationalFunctions(combination1.ratfuns)
    for pair in combination2.ratfuns
        combination += pair
    end
    return combination
end

function Base.:(*)(combination::CombinationOfRationalFunctions, c::Integer)
    newpairs = Vector{Pair{Polynomial,Polynomial}}([])
    for ratfun in combination.ratfuns
        push!(newpairs, Pair(ratfun[1] * c, ratfun[2]))
    end
    return CombinationOfRationalFunctions(newpairs)
end

function get_ratfun(c::CombinationOfRationalFunctions)
    println(parent(c.ratfuns[1]))
    FF = fraction_field(parent(c.ratfuns[1]))
    sum([n / d for (n, d) in c.ratfuns])
end

function simplify(c::CombinationOfRationalFunctions)
    return Symbolics.simplify(sum([n / d for (n, d) in c.ratfuns]))
end

function evaluate_all_with(c::CombinationOfRationalFunctions, variables::Vector, val::Any)
    return evaluate(c, Dict([variables[i] => val for i in eachindex(variables)]))
end

function evaluate(c::CombinationOfRationalFunctions)
    R = [n / d for (n, d) in c.ratfuns]
    return Symbolics.simplify(sum([r for r in R]))
end

function evaluate(c::CombinationOfRationalFunctions, values::Union{Dict,Pair})
    R = [n / d for (n, d) in c.ratfuns]
    return Symbolics.substitute(Symbolics.simplify(sum([r for r in R])), values)
end

function subs(c::CombinationOfRationalFunctions, v::Vector)
    R = [n / d for (n, d) in c.ratfuns]
    return sum([r(v) for r in R])
end

function Base.show(io::IO, c::CombinationOfRationalFunctions)
    print(io, "CombinationOfRationalFunctions: ")
    for (i, ratfun) in enumerate(c.ratfuns)
        if (i == length(c.ratfuns))
            print(io, "(($(ratfun[1]))/($(ratfun[2])))")
        else
            print(io, "(($(ratfun[1]))/($(ratfun[2]))) + ")
        end
    end
end
