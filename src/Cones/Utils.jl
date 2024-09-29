using PolyhedralOmegaMono.Values

function substitute_vector(vector::Vector{T}, subs::Dict) where {T<:Union{Number,Value}}
    return Vector{Number}([isa(x, Value) ? eval(substitute(x, subs)) : x for x in vector])
end

function substitute_ray(ray::Ray{T}, subs::Dict) where {T<:Union{Number,Value}}
    return Ray{Number}(
        substitute_vector(ray.direction, subs),
        substitute_vector(ray.apex, subs),
    )
end

function substitute_cone(cone::Cone{T}, subs::Dict) where {T<:Union{Number,Value}}
    return Cone{Number}(
        [substitute_ray(r, subs) for r in cone.rays],
        substitute_vector(cone.apex, subs),
        cone.openness,
        cone.sign,
    )
end
