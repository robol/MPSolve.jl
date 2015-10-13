using Polynomials
using MPSolve
using Base.Test

function test_roots_of_unity(n)
    nroots = [ exp(j * 1im * 2 * pi / n) for j = 1 : n ]
    
    p = [ Int64(0) for i = 1 : n + 1 ]
    p[1] = 1
    p[end] = -1

    p = Poly(p)

    (app, rad) = mps_roots(p)

    for i = 1 : n
        (err, ind) = findmin(abs(app - nroots[i]))
        @test err <= rad[ind]
    end
end

function test_roots_of_unity_fp(n)
    nroots = [ exp(j * 1im * 2 * pi / n) for j = 1 : n ]
    
    p = zeros(n+1)
    p[1] = 1.0
    p[end] = -1.0

    p = Poly(p)

    (app, rad) = mps_roots(p)

    for i = 1 : n
        (err, ind) = findmin(abs(app - nroots[i]))
        @test err <= rad[ind]
    end
end

function test_roots_of_unity_bigint(n)
    nroots = [ exp(j * 1im * 2 * pi / n) for j = 1 : n ]
    
    p = [ BigInt(0) for i = 1 : n + 1 ]
    p[1] = BigInt(1)
    p[end] = BigInt(-1)

    p = Poly(p)

    (app, rad) = mps_roots(p)

    for i = 1 : n
        (err, ind) = findmin(abs(app - nroots[i]))
        @test err <= rad[ind]
    end    
end

"""
Test if solving a polynomial with complex integer coefficienta
works. 
"""
function test_complex_int()
    sols = [ 2 ; 3+1im ; 5-1im ; -2 ]

    p = Poly([ 1 ])
    for i = 1 : length(sols)
        p = p * Poly([ - sols[i] ; 1 ])
    end

    (app,rad) = mps_roots(p)
    
    for i = 1 : length(sols)
        (err, ind) = findmin(abs(app - sols))
        @test err <= rad[ind]
    end        
end

function test_complex_bigint()
    sols = Complex{BigInt}[ 2 ; 3+1im ; 5-1im ; -2 ]

    p = Poly([ BigInt(1) ])
    for i = 1 : length(sols)
        p = p * Poly([ - sols[i] ; BigInt(1) ])
    end

    (app,rad) = mps_roots(p)
    
    for i = 1 : length(sols)
        (err, ind) = findmin(abs(app - sols))
        @test err <= rad[ind]
    end
end

test_roots_of_unity(100)
test_roots_of_unity_fp(100)
test_roots_of_unity_bigint(100)
test_complex_int()
test_complex_bigint()
