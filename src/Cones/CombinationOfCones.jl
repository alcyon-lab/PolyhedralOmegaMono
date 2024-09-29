struct CombinationOfCones{T}
    cones::Vector{Pair{Cone{T},Int}}

    function CombinationOfCones(cones::Vector{Pair{Cone{T},Int}}) where {T}
        new{T}(cones)
    end

    function CombinationOfCones{T}() where {T}
        new{T}([])
    end
end

function Base.:(+)(combination::CombinationOfCones{T}, cone::Cone{T}) where {T}
    newcones::Vector{Pair{Cone{T},Int}} = []
    found = false
    for (c, count) in combination.cones
        if cone == c
            found = true
            if c.sign == cone.sign
                push!(newcones, Pair(c, count + 1))
            else
                if count - 1 != 0
                    push!(newcones, Pair(c, count - 1))
                end
            end
        else
            push!(newcones, Pair(c, count))
        end
    end
    if !found
        push!(newcones, Pair(cone, 1))
    end
    return CombinationOfCones(newcones)
end

function Base.:(+)(combination::CombinationOfCones{T}, cone::Pair{Cone{T},Int}) where {T}
    c, count = cone
    for i in 1:count
        combination += c
    end
    return combination
end

function Base.:(+)(combination1::CombinationOfCones{T}, combination2::CombinationOfCones{T}) where {T}
    combination = CombinationOfCones(combination1.cones)
    for pair in combination2.cones
        combination += pair
    end
    return combination
end

function Base.:(+)(combination::CombinationOfCones{T}, cones::Vector{Cone{T}}) where {T}
    for cone in cones
        combination += cone
    end
    return combination
end

function Base.show(io::IO, c::CombinationOfCones{T}) where {T}
    println(io, "CombinationOfCone{::$(T)}: ")
    for (cone, count) in c.cones
        println(io, "\t$(summary(cone)) x $(count)")
    end
end
