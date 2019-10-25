# Testing for the MPSolve.jl package. 

using MPSolve
using Test

unity_roots(n) = [exp(j*2im*BigFloat(pi)/n) for j = 1:n]

function solve_test(p, rts)
    (app, rad) = mps_roots(p, 54)
    I = sortperm(app, by=angle)
    dist = map(abs, app[I] - sort(rts, by=angle))
    @test all(dist <= rad[I])
end

function test_roots_of_unity(n)
    p = Int64[0 for i = 1:n + 1]
    p[1] = 1
    p[end] = -1
    solve_test(p, unity_roots(n))
end

function test_roots_of_unity_fp(n)
    p = zeros(n+1)
    p[1] = 1.0
    p[end] = -1.0
    solve_test(p, unity_roots(n))
end

function test_roots_of_unity_bigint(n)
    p = [BigInt(0) for i = 1 : n + 1]
    p[1] = BigInt(1)
    p[end] = BigInt(-1)
    solve_test(p, unity_roots(n))
end

"""
Test if solving a polynomial with complex integer coefficients
works. 
"""
function roots2coeffs(roots)
    coeffs = [0im for n = 1:length(roots) + 1]
    coeffs[1] = 1
    for i = 1:length(roots)
        c = -roots[i]*coeffs[1:i + 1]
        c[2:i + 1] += coeffs[1:i]
        coeffs[1:i + 1] = c
    end
    coeffs
end

function test_complex_int()
    sols = [ 2 ; 3+1im ; 5-1im ; -2 ]
    p = roots2coeffs(sols)
    solve_test(p, sols)
end

function test_complex_bigint()
    sols = Complex{BigInt}[ 2 ; 3+1im ; 5-1im ; -2 ]
    p = roots2coeffs(sols)
    solve_test(p, sols)
end

function test_roots_of_unity_bigfloat(n)
    p = [ BigFloat(0) for i = 1 : n + 1 ]
    p[1] = 1
    p[end] = -1
    solve_test(p, unity_roots(n))
end

if VERSION < v"1.3-rc4"
    error("This module only works with Julia version 1.3-rc4 or greater")
end

N = 100

test_roots_of_unity(N)
test_roots_of_unity_fp(N)
test_roots_of_unity_bigint(N)
test_roots_of_unity_bigfloat(N)
test_complex_int()
test_complex_bigint()
