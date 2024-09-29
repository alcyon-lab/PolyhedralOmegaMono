using PolyhedralOmegaMono.Cones

function eliminate_last_coordinate(cone::Cone{T}; equality::Bool=false)::CombinationOfCones{T} where {T<:Union{Value,Number}}
    result = CombinationOfCones{T}()
    last_apex = cone.apex[end]
    indices = []
    if last_apex < 0
        for j in eachindex(cone.rays)
            if cone.rays[j].direction[end] > 0
                push!(indices, j)
            end
        end
    else
        for j in eachindex(cone.rays)
            if cone.rays[j].direction[end] < 0
                push!(indices, j)
            end
        end
        if !equality
            rays = collect(map(r -> Ray{T}(r.direction[1:end-1], r.apex[1:end-1]), cone.rays)) # CPrime
            result += flip(Cone(rays, cone.apex[1:end-1], cone.openness[:], cone.sign))
        end
    end
    if equality
        if ((cone.apex[end] < 0) || (cone.apex[end] == 0 && sum([r.direction[end] > 0 for r in cone.rays]) > 0) # exists(range(k),lambda j: V[j][n] > 0)[0]):
            || ((cone.apex[end] > 0) || (cone.apex[end] == 0 && sum([r.direction[end] < 0 for r in cone.rays]) > 0))
        )
            for j in indices
                ray_j = cone.rays[j]
                rays = Vector{Ray{T}}([])
                for (i, ray_i) in enumerate(cone.rays)
                    if i != j
                        sign = last_apex >= 0 ? 1 : -1
                        if i == j
                            ray = -sign * ray_i
                        else
                            ray = sign * (ray_j * ray_i.direction[end] - ray_i * ray_j.direction[end])
                        end
                        reduced_ray = primitive(Ray{T}(ray.direction[1:end-1], ray.apex[1:end-1]))
                        push!(rays, reduced_ray)
                    end
                end

                apex = Vector{T}(cone.apex - (last_apex // ray_j.direction[end]) .* ray_j.direction)
                reduced_apex = apex[1:end-1]

                openness = deepcopy(cone.openness)
                deleteat!(openness, j)

                new_cone = Cone(rays, reduced_apex, openness, cone.sign)
                new_cone = flip(new_cone)
                result += new_cone
            end
        else # q[n] == 0 and forall(range(k), lambda j: V[j][n] == 0)[0]
            rays = collect(map(r -> Ray{T}(r.direction[1:end-1], r.apex[1:end-1]), cone.rays)) # CPrime
            result += flip(Cone(rays, cone.apex[1:end-1], cone.openness[:], cone.sign))
        end
    else
        for j in indices
            ray_j = cone.rays[j]
            rays = Vector{Ray{T}}([])
            for (i, ray_i) in enumerate(cone.rays)
                sign = last_apex >= 0 ? 1 : -1
                if i == j
                    ray = -sign * ray_i
                else
                    ray = sign * (ray_j * ray_i.direction[end] - ray_i * ray_j.direction[end])
                end
                reduced_ray = primitive(Ray{T}(ray.direction[1:end-1], ray.apex[1:end-1]))
                push!(rays, reduced_ray)
            end

            apex = Vector{T}(cone.apex - (last_apex // ray_j.direction[end]) .* ray_j.direction)
            reduced_apex = apex[1:end-1]

            openness = deepcopy(cone.openness)
            openness[j] = false

            new_cone = Cone(rays, reduced_apex, openness, cone.sign)
            new_cone = flip(new_cone)
            result += new_cone
        end
    end
    return result
end

function eliminate_coordinates(cone::Cone{T}, k::Int, equality::Union{Vector{Bool},Nothing})::CombinationOfCones{T} where {T<:Union{Value,Number}}
    combination = CombinationOfCones{T}()
    combination += cone
    if isnothing(equality)
        equality = fill(false, (k,))
    end
    for i = 1:k
        innercombination = CombinationOfCones{T}()
        for (cone, count) in combination.cones
            res = eliminate_last_coordinate(cone, equality=equality[k-i+1])
            for j in 1:count
                innercombination += res
            end
        end
        combination = innercombination
    end
    return combination
end
