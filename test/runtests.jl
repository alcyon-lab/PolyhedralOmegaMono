using PolyhedralOmegaMono
using Test

@testset "PolyhedralOmegaMono.jl" begin
    @testset "PolyhedralOmega.jl" begin
        @testset "Square with segment" begin
            A = [1 0; 0 1; -1 0; 0 -1; 1 -1]
            b = [5, 5, -2, -2, 0]
            e = Vector{Bool}([0, 0, 0, 0, 1])
            r1, r2, r3 = PolyhedralOmegaMono.PolyhedralOmega.solve(A, b, e)
            expected_cone_count = 2
            expected_rat_fun_str = "CombinationOfRationalFunctions: (((x[1]^2)*(x[2]^2) - (x[1]^6)*(x[2]^6))/(1 - x[1]*x[2]))"
            @test length(r1.cones) == expected_cone_count
            @test string(r3) == expected_rat_fun_str
        end
        @testset "hyp7.1" begin
            A = -[
                -1 -1 -1 -1 -1 -1 -1;
                -1 0 0 0 0 0 0;
                0 -1 0 0 0 0 0;
                0 0 -1 0 0 0 0;
                0 0 0 -1 0 0 0;
                0 0 0 0 -1 0 0;
                0 0 0 0 0 -1 0;
                0 0 0 0 0 0 -1;
                1 0 0 0 0 0 0;
                0 1 0 0 0 0 0;
                0 0 1 0 0 0 0;
                0 0 0 1 0 0 0;
                0 0 0 0 1 0 0;
                0 0 0 0 0 1 0;
                0 0 0 0 0 0 1
            ]
            b = -[-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
            e = Vector{Bool}([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            r1, r2, r3 = PolyhedralOmegaMono.PolyhedralOmega.solve(A, b, e)
            expected_cone_count = 42
            expected_rat_fun_str = ""
            open("$(@__DIR__)/expected-outputs/hyp7.1", "r") do file
                expected_rat_fun_str = read(file, String)
            end
            @test length(r1.cones) == expected_cone_count
            @test string(r3) == expected_rat_fun_str
        end
        @testset "Magic4x4" begin
            A = -[[-1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0];
                [0 0 0 0 -1 -1 -1 -1 0 0 0 0 0 0 0 0];
                [0 0 0 0 0 0 0 0 -1 -1 -1 -1 0 0 0 0];
                [0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1];
                [-1 0 0 0 -1 0 0 0 -1 0 0 0 -1 0 0 0];
                [0 -1 0 0 0 -1 0 0 0 -1 0 0 0 -1 0 0];
                [0 0 -1 0 0 0 -1 0 0 0 -1 0 0 0 -1 0];
                [0 0 0 -1 0 0 0 -1 0 0 0 -1 0 0 0 -1];
                [-1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1];
                [0 0 0 -1 0 0 -1 0 0 -1 0 0 -1 0 0 0]]
            b = -[-100, -100, -100, -100, -100, -100, -100, -100, -100, -100]
            e = Vector{Bool}([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
            r1, r2, r3 = PolyhedralOmegaMono.PolyhedralOmega.solve(A, b, e)
            expected_cone_count = 196
            expected_rat_fun_str = ""
            open("$(@__DIR__)/expected-outputs/magic4x4", "r") do file
                expected_rat_fun_str = read(file, String)
            end
            @test length(r1.cones) == expected_cone_count
            @test string(r3) == expected_rat_fun_str
        end
    end
end
