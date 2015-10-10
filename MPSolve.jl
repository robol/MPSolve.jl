module MPSolve

export polyval, roots

function setupMPSolve(n :: Integer)
    ctx = ccall((:mps_context_new, "libmps"), Ptr{Void}, ())
    mp  = ccall((:mps_monomial_poly_new, "libmps"), 
                Ptr{Void}, (Ptr{Void}, Int), ctx, n)
    (ctx, mp)
end

function solvePoly(ctx, mp)
    ccall((:mps_context_set_input_poly, "libmps"), Void, 
          (Ptr{Void}, Ptr{Void}), 
          ctx, mp)
    ccall((:mps_mpsolve, "libmps"), Void, (Ptr{Void},), ctx)
end

function getFloatingPointRoots(ctx, mp, n)
    # Here we prepare some containers that will be used to pass a 
    # double ** to MPSolve, using the function mps_context_get_roots_d(). 
    # 
    # That function will allocate memory to represent the results
    # whose ownership will be transferred to Julia through the use 
    # of pointer_to_array. 
    app_container = Array(Ptr{Complex128}, 1)
    radius_container = Array(Ptr{Float64}, 1)
    app_container[1] = 0
    radius_container[1] = 0

    ccall ((:mps_context_get_roots_d, "libmps"), Void, 
           (Ptr{Void}, Ptr{Ptr{Complex128}}, Ptr{Ptr{Float64}}),
           ctx, app_container, radius_container)

    approximations = pointer_to_array(app_container[1], n, true)
    radius = pointer_to_array(radius_container[1], n, true)

    (approximations, radius)
end

function releaseMPSolveContext(ctx, mp)
    ccall((:mps_monomial_poly_free, "libmps"), Void, (Ptr{Void}, Ptr{Void}), 
          ctx, mp)
    ccall((:mps_context_free, "libmps"), Void, (Ptr{Void},), ctx)
end

"""
r = polyval (coefficients, x) evaluates a polynomial at a point.

 - coefficients is a vector containing the coefficients
 of the polynomial in decreasing degree order, so the
 leading coefficient is the first element of the vector. 
 - x is the point where the polynomial is evaluted. 
"""
function polyval(coefficients :: AbstractVector, x)
    n = length(coefficients)
    r = 0    
    for i = 1 : n - 1
        r = (r + coefficients[i]) * x
    end
    r = r + coefficients[n]
end

"""
(approximations, radii) = roots(coefficients) approximates
the roots of the polynomial p(x) defined by the coefficients
in the input vector. 

 - coefficients is a vector containing the coefficients, 
 with the leading coefficient first. 
"""
function roots(coefficients :: AbstractVector{Float64})

    n = length(coefficients) - 1

    (ctx, mp) = setupMPSolve(n)
    
    for i = 1 : n + 1
        ccall((:mps_monomial_poly_set_coefficient_d, "libmps"), Void, 
              (Ptr{Void}, Ptr{Void}, Int, Float64, Float64), 
              ctx, mp, n-i+1, real(coefficients[i]), imag(coefficients[i]))
    end

    solvePoly(ctx, mp)
    (approximations, radius) = getFloatingPointRoots(ctx, mp, n)
    releaseMPSolveContext(ctx, mp)

    (approximations, radius)
end

function roots(coefficients :: AbstractVector{Int64})
    n = length(coefficients) - 1

    (ctx, mp) = setupMPSolve(n)
    
    for i = 1 : n + 1
        ccall((:mps_monomial_poly_set_coefficient_int, "libmps"), Void, 
              (Ptr{Void}, Ptr{Void}, Int, Int64, Int64), 
              ctx, mp, n-i+1, real(coefficients[i]), imag(coefficients[i]))
    end

    solvePoly(ctx, mp)
    (approximations, radius) = getFloatingPointRoots(ctx, mp, n)
    releaseMPSolveContext(ctx, mp)

    (approximations, radius)
end

function roots(real_coefficients :: AbstractVector{BigFloat},
               imag_coefficients :: AbstractVector{BigFloat})
    n = length(real_coefficients) - 1
    (ctx, mp) = setupMPSolve(n)

    for i = 1 : n + 1
        println(string(real_coefficients[i]))
        ccall ((:mps_monomial_poly_set_coefficient_s, "libmps"), Void,
               (Ptr{Void}, Ptr{Void}, Int, Ptr{UInt8}, Ptr{UInt8}), 
               ctx, mp, n-i+1, string(real_coefficients[i]), 
               string(imag_coefficients[i]))
    end

    solvePoly(ctx, mp)

    # TODO: We should get the output as BigFloat numbers, instead
    # of truncating it to floating point. 
    (approximations, radius) = getFloatingPointRoots(ctx, mp, n)

    releaseMPSolveContext(ctx, mp)

    (approximations, radius)
end

end
