module MPSolve

using Polynomials

export mps_roots

# 
# The following is taken from the discussion that can be found
# at https://groups.google.com/forum/#!msg/julia-dev/uqp7LziUEfY/9RkiymZY5twJ. 
# 
immutable mpz_struct
   alloc::Cint
   size::Cint
   d::Ptr{Void}
end

Base.convert(::Type{mpz_struct}, i::BigInt) = mpz_struct(i.alloc, i.size, i.d)

type mpq_struct
     num::mpz_struct
     den::mpz_struct
     mpq_struct(num, den) = new(num, den)
     mpq_struct(q::Rational{BigInt}) = mpq_struct(q.num, q.den)
     mpq_struct() = new() # uninitialized
end

immutable mpf_struct
    mp_prec::Cint
    mp_size::Cint
    mp_exp::Cint
    mp_d::Ptr{Void}
end

type mpc_struct
    r::mpf_struct
    i::mpf_struct
    mpc_struct() = new() # uninitialized
end

type rdpe_struct
    r::Cdouble
    e::Clong
end

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

    ccall((:mps_context_get_roots_d, "libmps"), Void, 
          (Ptr{Void}, Ptr{Ptr{Complex128}}, Ptr{Ptr{Float64}}),
          ctx, app_container, radius_container)
    display(app_container[1])

    approximations = pointer_to_array(app_container[1], n, true)
    radius = pointer_to_array(radius_container[1], n, true)

    (approximations, radius)
end

function getGMPRoots(ctx, mp, n)
    # Obtain a copy of the approximations and radii represented
    # as Complex{BigFloat} and BigFloat, respectively. 

    app_container = Array(Ptr{mpc_struct}, 1)
    app_container[1] = 0

    radius_container = Array(Ptr{rdpe_struct}, 1)
    radius_container[1] = 0

    ccall((:mps_context_get_roots_m, "libmps"), Void, 
          (Ptr{Void}, Ptr{Ptr{mpc_struct}}, Ptr{Ptr{rdpe_struct}}),
          ctx, app_container, radius_container)

    display(app_container[1])
    # mpf_approximations = pointer_to_array(app_container[1], n, true)
    
    # Convert mpf_t approximatino to the internal mpfr type of 
    # Julia, so we can map them back to BigFloats
end

function releaseMPSolveContext(ctx, mp)
    ccall((:mps_monomial_poly_free, "libmps"), Void, (Ptr{Void}, Ptr{Void}), 
          ctx, mp)
    ccall((:mps_context_free, "libmps"), Void, (Ptr{Void},), ctx)
end

"""
(approximations, radii) = mps_roots(p) approximates
the roots of the polynomial p(x). 
"""
function mps_roots(p::Poly{Complex{Float64}})

    coefficients = p.a
    n = length(coefficients) - 1

    (ctx, mp) = setupMPSolve(n)
    
    for i = 1 : n + 1
        ccall((:mps_monomial_poly_set_coefficient_d, "libmps"), Void, 
              (Ptr{Void}, Ptr{Void}, Int, Float64, Float64), 
              ctx, mp, i-1, real(coefficients[i]), imag(coefficients[i]))
    end

    solvePoly(ctx, mp)
    (approximations, radius) = getFloatingPointRoots(ctx, mp, n)
    releaseMPSolveContext(ctx, mp)

    (approximations, radius)
end

mps_roots(p::Poly{Float64}) = mps_roots(convert(Poly{Complex{Float64}}, p))

function mps_roots(p::Poly{Complex{Rational{BigInt}}})

    real_coefficients = real(p.a)
    imag_coefficients = imag(p.a)

    n = length(real_coefficients) - 1
    (ctx, mp) = setupMPSolve(n)

    for i = 1 : n + 1
        x = mpq_struct(real_coefficients[i])
        y = mpq_struct(imag_coefficients[i])

        ccall ((:mps_monomial_poly_set_coefficient_q, "libmps"), Void,
               (Ptr{Void}, Ptr{Void}, Int, 
                Ptr{mpq_struct}, Ptr{mpq_struct}), 
               ctx, mp, i-1, &x, &y)
    end    

    solvePoly(ctx, mp)

    # TODO: We should get the output as BigFloat numbers, instead
    # of truncating it to floating point. 
    (approximations, radius) = getFloatingPointRoots(ctx, mp, n)
    # (approximations, radius) = getGMPRoots(ctx, mp, n)
    releaseMPSolveContext(ctx, mp)
    (approximations, radius)
end

# Generic methods defined using the conversion of Integer types to 
# BigInts and Rational{BigInt}s. 
mps_roots{T<:Integer}(p::Poly{Complex{T}}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))
mps_roots{T<:Integer}(p::Poly{T}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))
mps_roots{T<:Integer}(p::Poly{Rational{T}}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))

function mps_roots(p::Poly{Complex{BigFloat}})
    real_coefficients = real(p.a)
    imag_coefficients = imag(p.a)

    n = length(real_coefficients) - 1
    (ctx, mp) = setupMPSolve(n)

    # TODO: Finish to implement this function
end

mps_roots(p::Poly{BigFloat}) = mps_roots(convert(Poly{Complex{BigFloat}}, p))

end ## End of Module MPSolve
