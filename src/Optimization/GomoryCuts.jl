function gomory_cuts(A::Matrix{<:Number}, b::Vector{<:Number}, c::Vector{<:Number})
    return gomory_cuts(Matrix{Rational}(A), Vector{Rational}(b), Vector{Rational}(c))
end
function gomory_cuts(A::Matrix{Rational}, b::Vector{Rational}, c::Vector{Rational})
    m, n = size(A)
    slack_var_count = m
    Tableau = Matrix{Rational}(undef, m + 1, n + slack_var_count + 1)

    Tableau[1:m, 1:n] = A
    Tableau[1:m, n+slack_var_count+1] = b
    for i in 1:m
        for j in n+1:n+slack_var_count
            if i == j - n
                Tableau[i, j] = 1
            else
                Tableau[i, j] = 0
            end
        end
    end
    # objective row
    Tableau[m+1, 1:n] = c
    Tableau[m+1, n+1:n+slack_var_count+1] .= 0
    OriginalTableau = deepcopy(Tableau)
    Tableau = _simplex_method_row_ops(Tableau)
    slack_cuts = _compute_gomory_cuts(Tableau)
    cuts = Vector{Rational}[]
    for cut in slack_cuts
        res = cut
        for j in 1:slack_var_count
            slack_coeff = cut[j+n]
            res = res - (slack_coeff * OriginalTableau[j, :])
        end
        if sum(res) != 0
            println(sum(res))
            push!(cuts, -res)
        end

    end
    newA = Matrix{Rational}(undef, m + length(cuts), n)
    newB = Vector{Rational}(undef, m + length(cuts))
    newA[1:m, :] = A
    newB[1:m] = b
    for (i, cut) in enumerate(cuts)
        newA[m+i, 1:n] = cut[1:n]
        newB[m+i] = cut[n+slack_var_count+1]
    end
    return newA, newB, cuts
end

function _simplex_method_row_ops(Tableau::Matrix{Rational})
    m, n = size(Tableau)
    while true
        obj_row = Tableau[m, :]
        pivot_col = argmin(obj_row)
        if obj_row[pivot_col] >= 0
            return Tableau
        end
        if any(e -> e < 0, Tableau[:, pivot_col])
            pivot_row = argmin(Tableau[begin:(m-1), n] ./ Tableau[begin:(m-1), pivot_col])
            Tableau[pivot_row, :] = Tableau[pivot_row, :] / Tableau[pivot_row, pivot_col]
            for i in 1:m
                if i != pivot_row
                    Tableau[i, :] = Tableau[i, :] - Tableau[pivot_row, :] * Tableau[i, pivot_col]
                end
            end
        else
            throw("No optimal solution")
        end
    end
end

function _compute_gomory_cuts(Tableau::Matrix{Rational})
    m, n = size(Tableau)
    res = Vector{Rational}[]
    for i in 1:m-1
        nums = Int[]
        denoms = Int[]

        for j in 1:n
            denom = denominator(Tableau[i, j])
            push!(nums, mod(numerator(Tableau[i, j]), denom))
            push!(denoms, denom)
        end
        push!(res, nums .// denoms)
    end
    return res
end

""" visualization
using Plots
using ImplicitEquations

function generate_equality_functions(A, b)
    f_array = Array{Function}(undef, length(b))
    for j in eachindex(b)
        f_array[j] = x -> (b[j] - A[j, 1] * x) / A[j, 2]
    end
    return f_array
end

function plot_gomory(A, b, c)
    nA, nb, cuts = gomory_cuts(A, b, c)
    p = plot(framestyle=:origin)
    plot!(p, generate_equality_functions(nA, nb), xlims=(-1, 10), ylims=(-1, 10), lc=:red)
    plot!(p, generate_equality_functions(A, b), xlims=(-1, 10), ylims=(-1, 10), lc=:blue)
    display(p)
end
"""
