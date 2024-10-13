using SmithNormalForm
using IterTools
using LinearAlgebra

function enumerate_fundamental_parallelepiped(cone::Cone{T}) where {T<:Number}
    vrep = Matrix{Int}(vrep_matrix(cone))
    SNFRes = SmithNormalForm.smith(vrep)
    S = SmithNormalForm.diagm(SNFRes)
    Uinv = SNFRes.Sinv
    Winv = SNFRes.Tinv
    dimension = size(vrep, 2) # num of rows
    ambientDimension = size(vrep, 1) # num of cols
    diagonals = Int64[]
    apex = cone.apex
    p = zeros(T, (ambientDimension,))
    if dimension < ambientDimension
        @info "vrep $(vrep)"
        @info "dim $(dimension) ambDim $(ambientDimension)"
        A_ = Matrix{Int}(round.(det(vrep' * vrep) * pinv(vrep)))

        Ap = A_ # [end-(ambientDimension-dimension)+1:end, :]
        b = Ap * apex
        A, b = Ap, b

        SNFRes_ = SmithNormalForm.smith(A)
        S_, Uinv_, Vinv_ = SmithNormalForm.diagm(SNFRes_), SNFRes_.Sinv, SNFRes_.Tinv
        @info "A $(A)"
        @info "b $(b)"
        @info "S_ $(S_)"
        @info "Vinv_ $(Vinv_)"
        @info "Uinv_ $(Uinv_)"
        m, n = size(A)
        p_ = Uinv_ * b
        q = vcat([p_[i] // S_[i, i] for i in 1:m], zeros(Int, n - m))
        p = Vinv_ * q
        if !all(isinteger.(p))
            return Int[]
        end
    end
    for i in 1:dimension
        if (i <= size(S, 1) && i <= size(S, 2) && S[i, i] != 0)
            push!(diagonals, S[i, i])
        end
    end
    lastDiagonal = diagonals[end]
    # sprime = [Integer(sk / si) for si in s]
    sprimeDiagonals = Int64[]
    for d in diagonals
        push!(sprimeDiagonals, Int64(lastDiagonal // d))
    end
    for _ in 1:(dimension-length(diagonals))
        push!(sprimeDiagonals, 1)
    end
    sprime = diagm(sprimeDiagonals)
    qhat = Uinv * (apex - p)
    Wprime = Winv * sprime
    qtrans = [sum([-Wprime[j, i] * qhat[i] for i = 1:dimension]) for j = 1:dimension]
    qfrac = [qtrans[i] - floor(Int, qtrans[i]) for i = 1:dimension]
    qint = [floor(Int, qi) for qi in qtrans]
    qsummand = [Int64(qi) for qi in (lastDiagonal * apex + vrep * qfrac)]
    openness = [(qfrac[j] == 0 ? cone.openness[j] : false) for j in 1:dimension]
    #bigP
    #res1 = [[1:1:diagonals[i];] for i= 1:dimension]
    # CartesianProduct( *[xrange(s[i]) for i in 1:k] )
    f = let qint = qint, Wprime = Wprime, openness = openness, lastDiagonal = lastDiagonal, ambientDimension = ambientDimension, vrep = vrep
        function (v)
            innerRes = mod.(Wprime * v .+ qint, lastDiagonal)
            innerRes = [inner == 0 && openness[j] ? lastDiagonal : inner for (inner, j) in zip(innerRes, 1:length(openness))]
            outerRes = vrep * innerRes
            finalRes = collect(Int64(ai + bi) // lastDiagonal for (ai, bi) in zip(outerRes, qsummand))
            return finalRes
        end
    end
    return (f(collect(x)) for x in IterTools.product([(0:diagonals[i]-1) for i in 1:dimension]...))
end
