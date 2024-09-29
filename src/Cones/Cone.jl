struct Cone{T}
    rays::Vector{Ray{T}}
    apex::Vector{T}
    openness::Vector{Bool}

    # Extra
    sign::Bool

    function Cone(rays::Vector{Vector{T}}, apex::Vector{T}, openness::Vector{Bool}, sign::Bool=true) where {T}
        ray_objects = [Ray{T}(ray) for ray in rays]
        sorted_rays, sorted_openness = _sort_rays(ray_objects, openness)
        new{T}(sorted_rays, apex, sorted_openness, sign)
    end

    function Cone(rays::Vector{Ray{T}}, apex::Vector{T}, openness::Vector{Bool}, sign::Bool=true) where {T}
        sorted_rays, sorted_openness = _sort_rays(rays, openness)
        new{T}(sorted_rays, apex, sorted_openness, sign)
    end

    function Cone{T}(rays::Vector{Ray{R}}, apex::Vector{A}, openness::Vector{Bool}, sign::Bool=true) where {T,R,A}
        ray_objects = [Ray{T}(ray) for ray in rays]
        sorted_rays, sorted_openness = _sort_rays(ray_objects, openness)
        new{T}(sorted_rays, convert(Vector{T}, apex), sorted_openness, sign)
    end

    function Cone{T}(rays::Vector{Vector{R}}, apex::Vector{A}, openness::Vector{Bool}, sign::Bool=true) where {T,R,A}
        ray_objects = [Ray{T}(ray) for ray in rays]
        sorted_rays, sorted_openness = _sort_rays(ray_objects, openness)
        new{T}(sorted_rays, convert(Vector{T}, apex), sorted_openness, sign)
    end

    function _sort_rays(rays::Vector{Ray{T}}, openness::Vector{Bool})::Tuple{Vector{Ray{T}},Vector{Bool}} where {T}
        paired = sort!([(primitive(rays[i]), openness[i]) for i in 1:length(rays)], lt=(x, y) -> isless(x[1].direction, y[1].direction))
        sorted_rays = [p[1] for p in paired]
        sorted_openness = [p[2] for p in paired]
        return sorted_rays, sorted_openness
    end
end

function Base.:(==)(c1::Cone{T1}, c2::Cone{T2}) where {T1,T2}
    # TODO: fix ray check
    return c1.apex == c2.apex && c1.rays == c2.rays && c1.openness == c2.openness
end


function Base.summary(io::IO, c::Cone{T}) where {T}
    print(io, "Cone: ")
    print(io, "{ R: ")
    for ray in c.rays
        print(io, "$(ray.direction)")
    end
    print(io, " ")
    print(io, "A: $(c.apex) ")
    print(io, "O: $(c.openness) ")
    print(io, "($(c.sign ? '+' : '-')) }")
end

function Base.show(io::IO, c::Cone{T}) where {T}
    println(io, "Cone{::$(T)}")
    println(io, "Rays:")
    for ray in c.rays
        println(io, "\t$(ray)")
    end
    println(io, "Apex:")
    println(io, "\t$(c.apex)")
    println(io, "Openness:")
    println(io, "\t$(c.openness)")
    println(io, "Sign:")
    print(io, "\t$(c.sign ? '+' : '-')")
end

function vrep_matrix(cone::Cone{T}) where {T}
    if isempty(cone.rays) || isempty(cone.apex)
        throw(ArgumentError("Cone must have rays and apex defined for V-rep."))
    end

    dim = length(cone.rays[1].direction)
    vrep = Matrix{T}(undef, dim, length(cone.rays))

    for (i, ray) in enumerate(cone.rays)
        vrep[:, i] = ray.direction
    end

    return vrep
end


function flip(cone::Cone{T}) where {T}
    sign::Bool = false
    openness = deepcopy(cone.openness)
    rays = deepcopy(cone.rays)
    for i in 1:size(cone.rays, 1)
        if !isforward(cone.rays[i])
            # change direction 
            rays[i] = flip(cone.rays[i])
            # change openness
            openness[i] = !cone.openness[i]
            # change sign   
            sign = !sign
        end
    end
    return Cone(rays, cone.apex, openness, xor(sign, cone.sign))
end
