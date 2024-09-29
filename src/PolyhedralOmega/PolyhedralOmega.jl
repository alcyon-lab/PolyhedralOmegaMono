module PolyhedralOmega

using PolyhedralOmegaMono.Polynomials
using PolyhedralOmegaMono.Cones
using PolyhedralOmegaMono.Values
using Symbolics

include("EliminateCoordinates.jl")
include("MacmahonLifting.jl")

export solve, optimize


function solve(A::Matrix{T}, b::Vector{T}, equality::Union{Vector{Bool},Nothing}=nothing; write_rf_to_out::Bool=false, out::IO=stdout, counting::Bool=false, enable_multi_thread::Bool=false) where {T<:Number}
    return solve(Matrix{Rational}(A), Vector{Rational}(b), equality, write_rf_to_out=write_rf_to_out, out=out, counting=counting)
end

function solve(A::Matrix{T}, b::Vector{T}, equality::Union{Vector{Bool},Nothing}=nothing; write_rf_to_out::Bool=false, out::IO=stdout, counting::Bool=false, enable_multi_thread::Bool=false) where {T<:Union{Value,Rational}}
    A = -A
    b = -b
    macmahon_cone = macmahon_lifting(A, b)
    list_of_cones = eliminate_coordinates(macmahon_cone, size(b, 1), equality)
    fpps = Dict()
    global r = CombinationOfRationalFunctions()
    @variables x[1:size(A, 2)]
    @variables q
    first_loop = true
    if enable_multi_thread
        lk = ReentrantLock()
        Threads.@threads for (cone, count) in list_of_cones.cones
            cone = Cone{Number}(cone.rays, cone.apex, cone.openness, cone.sign)
            fpp = enumerate_fundamental_parallelepiped(cone)
            lock(lk) do
                fpps[cone] = fpp
                if write_rf_to_out
                    cone_rf_str = compute_rational_function_str(cone, fpp, counting=counting)
                    if count != 1
                        cone_rf_str = "($(count)*($(cone_rf_str)))"
                    end
                    if first_loop
                        write(out, cone_rf_str)
                        first_loop = false
                    else
                        write(out, " + ")
                        write(out, cone_rf_str)
                    end
                else
                    if counting
                        cone_rf_s = compute_rational_function(cone, fpp, [q for i in 1:size(A, 2)]) * count
                    else
                        cone_rf_s = compute_rational_function(cone, fpp, collect(x[1:size(A, 2)])) * count
                    end
                    r = cone_rf_s + r
                end
            end
        end
    else
        for (cone, count) in list_of_cones.cones
            cone = Cone{Number}(cone.rays, cone.apex, cone.openness, cone.sign)
            fpp = enumerate_fundamental_parallelepiped(cone)
            fpps[cone] = fpp
            if write_rf_to_out
                cone_rf_str = compute_rational_function_str(cone, fpp, counting=counting)
                if count != 1
                    cone_rf_str = "($(count)*($(cone_rf_str)))"
                end
                if first_loop
                    write(out, cone_rf_str)
                    first_loop = false
                else
                    write(out, " + ")
                    write(out, cone_rf_str)
                end
            else
                if counting
                    cone_rf_s = compute_rational_function(cone, fpp, [q for i in 1:size(A, 2)]) * count
                else
                    cone_rf_s = compute_rational_function(cone, fpp, collect(x[1:size(A, 2)])) * count
                end
                r = cone_rf_s + r
            end
        end
    end

    if write_rf_to_out
        return list_of_cones, fpps, nothing
    else
        return list_of_cones, fpps, r
    end
end

function optimize(A::Matrix{T}, b::Vector{T}, f::Union{Vector{T},Nothing}=nothing; max_value::Number=32, equality::Union{Vector{Bool},Nothing}=nothing) where {T<:Number}
    newf::Union{Vector{Rational},Nothing} = nothing
    if isnothing(f)
        newf = f
    else
        newf = Vector{Rational}(f)
    end
    return optimize(Matrix{Rational}(A), Vector{Rational}(b), newf, max_value=max_value, equality=equality)
end

function optimize(A::Matrix{T}, b::Vector{T}, f::Union{Vector{T},Nothing}=nothing; max_value::Number=32, equality::Union{Vector{Bool},Nothing}=nothing) where {T<:Rational}
    newf::Union{Vector{Value},Nothing} = nothing
    if isnothing(f)
        newf = f
    else
        newf = Vector{Value}(f)
    end
    return optimize(Matrix{Value}(A), Vector{Value}(b), newf, max_value=max_value, equality=equality)
end

function optimize(A::Matrix{T}, b::Vector{T}, f::Union{Vector{T},Nothing}=nothing; max_value::Number=32, equality::Union{Vector{Bool},Nothing}=nothing) where {T<:Value}
    A = -A
    b = -b
    α = Symbol("α")
    if isnothing(f)
        f = ones(T, (size(A, 2)))
    end
    macmahon_cone = macmahon_lifting(A, b, f, symbol=α)
    list_of_cones = eliminate_coordinates(macmahon_cone, size(b, 1), equality)
    value = max_value // 2
    min_value = 0
    optimal_rf = (-1 => Value)
    capacity = -1
    @variables x[1:size(A, 2)]
    while true
        print("Value: ", value)
        rf = CombinationOfRationalFunctions()
        for (cone, count) in list_of_cones.cones
            s_cone = substitute_cone(cone, Dict(α => value))
            s_cone_eliminated = eliminate_last_coordinate(s_cone)
            for (s, s_count) in s_cone_eliminated.cones
                fpp = enumerate_fundamental_parallelepiped(substitute_cone(s, Dict()))
                rf += (compute_rational_function(s, fpp, collect(x)) * count * s_count)
            end
        end
        eval_res = evaluate_all_with(rf, collect(x), 1)
        res = floor(eval_res)
        println(" | Res: ", res)
        if isequal(res, 1)
            return Polynomials.simplify(rf)
        elseif isequal(res, 0)
            tmp_value = value
            value = floor(min_value + (value - min_value) / 2)
            max_value = tmp_value
            capacity = tmp_value
            if isequal(value, min_value)
                if optimal_rf[1] != -1
                    return Polynomials.simplify(optimal_rf[2])
                end
                return Polynomials.simplify(rf)
            end
        else
            tmp_value = value
            if capacity == -1
                value *= 2
            else
                value = floor(value + (capacity - value) / 2)
            end
            min_value = tmp_value
            if optimal_rf[1] == -1 || optimal_rf[1] > res
                optimal_rf = res => rf
            end
            if isequal(value, min_value)
                return Polynomials.simplify(rf)
            end
        end
    end
    assert("Logic error")
end

end
