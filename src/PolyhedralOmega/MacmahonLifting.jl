using LinearAlgebra # for I
using PolyhedralOmegaMono.Cones

function macmahon_lifting(A::Matrix{T}, b::Vector{T})::Cone{T} where {T<:Union{Number,Value}}
    size_a = size(A)
    Id = Matrix{T}(Matrix(1I, size_a[2], size_a[2]))
    new_matrix = vcat(Id, A)
    apex = append!(zeros(T, size_a[2]), -b)
    openness = zeros(Bool, size_a[2])
    rays = convert(Vector{Vector{T}}, collect(eachcol(new_matrix)))
    return Cone(rays, apex, openness)
end

function macmahon_lifting(A::Matrix{T}, b::Vector{T}, f::Vector{T}; symbol::Symbol=Symbol('a'))::Cone{T} where {T<:Value}
    size_a = size(A)
    Id = Matrix{T}(Matrix(1I, size_a[2], size_a[2]))
    new_matrix = vcat(Id, transpose(f), A)
    apex = append!(zeros(T, size_a[2]), [-Expr(symbol)], -b)
    openness = zeros(Bool, size_a[2])
    rays = convert(Vector{Vector{T}}, collect(eachcol(new_matrix)))
    return Cone(rays, apex, openness)
end
