using PolyhedralOmegaMono.Values

struct Ray{T}
    direction::Vector{T}
    apex::Vector{T}

    function Ray(direction::Vector{T}) where {T}
        new{T}(direction, zeros(T, size(direction, 1)))
    end

    function Ray(direction::Vector{T}, apex::Vector{T}) where {T}
        new{T}(direction, apex)
    end

    function Ray{T}(direction::Vector) where {T}
        new{T}(convert(Vector{T}, direction), zeros(T, size(direction, 1)))
    end

    function Ray{T}(direction::Vector, apex::Vector) where {T}
        new{T}(convert(Vector{T}, direction), convert(Vector{T}, apex))
    end

    function Ray{T}(ray::Ray) where {T}
        new{T}(ray.direction, ray.apex)
    end

    function Ray(ray::Ray{T}) where {T}
        new{T}(ray.direction, ray.apex)
    end
end

function Base.:+(ray1::Ray{T}, ray2::Ray) where {T}
    return Ray{T}(ray1.direction + ray2.direction, ray1.apex + ray2.apex)
end

function Base.:-(ray1::Ray{T}, ray2::Ray) where {T}
    return Ray{T}(ray1.direction - ray2.direction, ray1.apex - ray2.apex)
end

function Base.:*(ray::Ray{T}, scalar::Number) where {T}
    return Ray{T}(ray.direction .* scalar, ray.apex)
end

function Base.:*(scalar::Number, ray::Ray{T}) where {T}
    return ray * scalar
end

function Base.:*(ray::Ray{T}, scalar::Value) where {T}
    return Ray{T}(ray.direction .* scalar.val, ray.apex)
end

function Base.:*(scalar::Value, ray::Ray{T}) where {T}
    return ray * scalar.val
end

function Base.:-(ray::Ray{T}) where {T}
    return ray * -1
end

function Base.:(==)(r1::Ray{T1}, r2::Ray{T2}) where {T1,T2}
    return r1.direction == r2.direction && r1.apex == r2.apex
end

function Base.show(io::IO, r::Ray{T}) where {T}
    direction = join(repr.(r.direction), ", ")
    apex = join(repr.(r.apex), ", ")
    print(io, "Ray{::$(T)} ($(direction)) apex: ($(apex))")
end

function flip(ray::Ray{T}) where {T}
    return Ray{T}([-e for e in ray.direction], ray.apex)
end

function isforward(ray::Ray)
    for e in ray.direction
        if e == 0
            continue
        elseif e > 0
            return true
        else
            return false
        end
    end
    # if ray is 0
    @warn "isforward used on a ray with no direction"
    return true
end

function primitive(ray::Ray{T}) where {T}
    @assert !(T <: AbstractFloat) "primitive used on a float ray (Ray{$(T)})"
    filtered_coordinates = collect(map(x -> Number(x), filter(x -> isinteger(x) || x isa Number, ray.direction)))
    if isempty(filtered_coordinates)
        return ray
    end
    g = gcd(filtered_coordinates)
    if g == 0
        g = 1
    end
    # Throw exception if convert(T,e/g) should not happen
    return Ray{T}(convert(Vector{T}, [e // g for e in ray.direction]), ray.apex)
end
