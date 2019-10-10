module MPSolve

using Polynomials

import Base: convert,show
import Base.GMP: BigInt
using Base.GMP: Limb, MPZ, BigInt

#export mps_roots
export MpzStruct,MpqStruct,prnt,clear,BigInt # ,mpz_t

struct MpzStruct
    alloc::Cint
    size::Cint
    d::Ptr{Limb}

    MpzStruct() = MpzStruct(0)
    MpzStruct(b::T) where T<:Signed = MpzStruct(BigInt(b))
    function MpzStruct(b::BigInt)
        m = Ref{MpzStruct}()
        GC.@preserve m begin
            ccall((:__gmpz_init, :libgmp), Cvoid, (Ref{MpzStruct},), m)
            ccall((:__gmpz_set, :libgmp), Cvoid, (Ref{MpzStruct},Ref{BigInt}),
                  m, b)
            return m[]
        end
    end
end

function BigInt(m::MpzStruct)
    b = BigInt()
    ccall((:__gmpz_set, :libgmp), Cvoid, (Ref{BigInt},Ref{MpzStruct}), b, m)
    return b
end

function prnt(m::MpzStruct)
    ccall((:__gmp_printf,:libgmp), Cvoid, (Cstring,Ref{MpzStruct},), "%Zd",m)
end

function show(io::IO, m::MpzStruct)
    b = BigInt(m)
    print(io,"MpzStruct(", b, ")")
end

function prnt(m::MpzStruct)
    ccall((:__gmp_printf,:libgmp), Cvoid, (Cstring,Ref{MpzStruct},), "%Zd\n",m)
end

clear(m::MpzStruct) = begin
    ccall((:__gmpz_clear, :libgmp), Cvoid, (Ref{MpzStruct},), m)
end

# const mpz_t = Ref{MpzStruct}
             
mutable struct MpqStruct
    num::MpzStruct
    den::MpzStruct

    MpqStruct() = MpqStruct(MpzStruct(), MpzStruct())
    MpqStruct(q::Rational{T}) where T<: Signed = MpqStruct(q.num, q.den)
    function MpqStruct(num, den)
        q = new(MpzStruct(num), MpzStruct(den))
        finalizer(x-> begin
                  clear(x.num)
                  clear(x.den)
                  end, q)
    end
end

show(io::IO, q::MpqStruct) = begin
    show(io, q.num)
    print(io, "//")
    show(io, q.den)
end

# Arbitrary precision floating point type from gmp.h
# is different from Julia's BigFloat
struct MpfStruct
    mp_prec::Cint
    mp_size::Cint
    mp_exp::Clong
    mp_d::Ptr{Cvoid}
end

mutable struct MpscStruct
    r::MpfStruct
    i::MpfStruct
    mpsc_struct() = new()
end

# struct RdpeStruct
#     r::Cdouble
#     e::Clong
# end

# mutable struct Cplx_struct
#      r::Cdouble
#      i::Cdouble
# end

function setupMPSolve(degree :: Integer)
    ctx = ccall((:mps_context_new, "libmps"), Ptr{Cvoid}, ())
    mp  = ccall((:mps_monomial_poly_new, "libmps"), 
                            Ptr{Cvoid}, (Ptr{Cvoid}, Int), ctx, degree)
    (ctx, mp)
end

function solve_poly(cf::Vector{Float64})
    degree = length(cf) - 1
    
end
    

# function solve_poly(ctx, mp)
#     ccall((:mps_context_set_input_poly, "libmps"), Cvoid, 
#           (Ptr{Cvoid}, Ptr{Cvoid}),
#           ctx, mp)
#     print("Enter mps_mpsolve()\n")
#     ccall((:mps_mpsolve, "libmps"), Cvoid, (Ptr{Cvoid},), ctx)
#     print("Leave mps_mpsolve()\n")
# end

# function get_floating_point_roots(ctx, mp, n)
#     # Here we prepare some containers that will be used to pass a 
#     # double ** to MPSolve, using the function mps_context_get_roots_d(). 
#     # 
#     # That function will allocate memory to represent the results
#     # whose ownership will be transferred to Julia through the use 
#     # of pointer_to_array. 
#     app_container = Array{Ptr{Nothing}}(undef, 1)
#     radius_container = Array{Ptr{Nothing}}(undef, 1)
#     app_container[1] = C_NULL
#     radius_container[1] = C_NULL
#     #app_container = Ref(C_NULL)
#     #radius_container = Ref(C_NULL)
    
#     ccall((:mps_context_get_roots_d, "libmps"), Cvoid, 
#           (Ptr{Cvoid}, Ref{Ptr{Nothing}}, Ref{Ptr{Nothing}}),
#           ctx, app_container, radius_container)
#     # display(app_container[1])

#     # approximations = unsafe_wrap(app_container[1], n, true)
#     # radius = unsafe_wrap(radius_container[1], n, true)

#     # (approximations, radius)
# end

# # function getGMPRoots(ctx, mp, n)
# #     # Obtain a copy of the approximations and radii represented
# #     # as Complex{BigFloat} and BigFloat, respectively. 

# #     app_container = Array(Ptr{MpscStruct}, 1)
# #     app_container[1] = 0

# #     radius_container = Array(Ptr{RdpeStruct}, 1)
# #     radius_container[1] = 0

# #     ccall((:mps_context_get_roots_m, "libmps"), Cvoid, 
# #           (Ptr{Cvoid}, Ptr{Ptr{MpscStruct}}, Ptr{Ptr{RdpeStruct}}),
# #           ctx, app_container, radius_container)

# #     display(app_container[1])
# #     # mpf_approximations = pointer_to_array(app_container[1], n, true)
    
# #     # Convert mpf_t approximatino to the internal mpfr type of 
# #     # Julia, so we can map them back to BigFloats
# # end

# function release_mpsolve_context(ctx, mp)
#     ccall((:mps_monomial_poly_free, "libmps"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), 
#           ctx, mp)
#     ccall((:mps_context_free, "libmps"), Cvoid, (Ptr{Cvoid},), ctx)
# end

# # """
# # (approximations, radii) = mps_roots(p) approximates
# # the roots of the polynomial p(x). 
# # """
# function mps_roots(p::Poly{Complex{Float64}})

#     coefficients = p.a
#     n = length(coefficients) - 1

#     (ctx, mp) = setupMPSolve(n)
    
#     for i = 1 : n + 1
#         ccall((:mps_monomial_poly_set_coefficient_d, "libmps"), Cvoid, 
#               (Ptr{Cvoid}, Ptr{Cvoid}, Int, Float64, Float64), 
#               ctx, mp, i-1, real(coefficients[i]), imag(coefficients[i]))
#     end

#     solve_poly(ctx, mp)
#     (approximations, radius) = get_floating_point_roots(ctx, mp, n)
#     release_mpsolve_context(ctx, mp)

#     (approximations, radius)
# end

# mps_roots(p::Poly{Float64}) = mps_roots(convert(Poly{Complex{Float64}}, p))

# function mps_roots(p::Poly{Complex{Rational{BigInt}}})

#     real_coefficients = real(p.a)
#     imag_coefficients = imag(p.a)

#     n = length(real_coefficients) - 1
#     (ctx, mp) = setupMPSolve(n)

#     for i = 1 : n + 1
#         x = MpqStruct(real_coefficients[i])
#         y = MpqStruct(imag_coefficients[i])

#         ccall ((:mps_monomial_poly_set_coefficient_q, "libmps"), Cvoid,
#                (Ptr{Cvoid}, Ptr{Cvoid}, Int, 
#                 Ptr{MpqStruct}, Ptr{MpqStruct}), 
#                ctx, mp, i-1, &x, &y)
#     end    

#     solvePoly(ctx, mp)

#     # TODO: We should get the output as BigFloat numbers, instead
#     # of truncating it to floating point. 
#     (approximations, radius) = getFloatingPointRoots(ctx, mp, n)
#     # (approximations, radius) = getGMPRoots(ctx, mp, n)
#     releaseMPSolveContext(ctx, mp)
#     (approximations, radius)
# end

# # Generic methods defined using the conversion of Integer types to 
# # BigInts and Rational{BigInt}s. 
# mps_roots{T<:Integer}(p::Poly{Complex{T}}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))
# mps_roots{T<:Integer}(p::Poly{T}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))
# mps_roots{T<:Integer}(p::Poly{Rational{T}}) = mps_roots(convert(Poly{Complex{Rational{BigInt}}}, p))

# function mps_roots(p::Poly{Complex{BigFloat}})
#     real_coefficients = real(p.a)
#     imag_coefficients = imag(p.a)

#     n = length(real_coefficients) - 1
#     (ctx, mp) = setupMPSolve(n)

#     # TODO: Finish to implement this function
# end

# mps_roots(p::Poly{BigFloat}) = mps_roots(convert(Poly{Complex{BigFloat}}, p))

end
