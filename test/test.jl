using MPSolve
using Base.Test

function test_roots_of_unity(n)
    nroots = exp(1im * 2 * pi / n)
    
    p = zeros(n+1)
    p[1] = 1
    p[end] = -1

    (app, rad) = roots(p)

    for i = 1 : n
        (err, ind) = findmin(abs(app - nroots))
        @test err <= rad[ind]
    end
end

function test_roots_of_unity_fp(n)
    nroots = exp(1im * 2 * pi / n)
    
    p = zeros(n+1)
    p[1] = 1.0
    p[end] = -1.0

    (app, rad) = roots(p)

    for i = 1 : n
        (err, ind) = findmin(abs(app - nroots))
        @test err <= rad[ind]
    end
end

test_roots_of_unity(100)
test_roots_of_unity_fp(100)
