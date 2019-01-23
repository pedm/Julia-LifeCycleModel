
## Bisection method for root finding (for FLoat64 values)

## From Jason Merrill https://gist.github.com/jwmerrill/9012954
## cf. http://squishythinking.com/2014/02/22/bisecting-floats/
# Alternative "mean" definition that operates on the binary representation
# of a float. Using this definition, bisection will never take more than
# 64 steps.

function _middle(x::Float64, y::Float64)
  # Use the usual float rules for combining non-finite numbers
  if !isfinite(x) || !isfinite(y)
    return x + y
  end

  # Always return 0.0 when inputs have opposite sign
  if sign(x) != sign(y) && x != 0.0 && y != 0.0
    return 0.0
  end

  negate = x < 0.0 || y < 0.0

  xint = reinterpret(UInt64, abs(x))
  yint = reinterpret(UInt64, abs(y))
  unsigned = reinterpret(Float64, (xint + yint) >> 1)

  return negate ? -unsigned : unsigned
end

"""

    `bisection64_custom(f, a, b)`

* `f`: a callable object, like a function

* `a`, `b`: Real values specifying a *bracketing* interval (one with
`f(a) * f(b) < 0`). These will be converted to `Float64` values.

Runs the bisection method using midpoints determined by a trick
leveraging 64-bit floating point numbers. After ensuring the
intermediate bracketing interval does not straddle 0, the "midpoint"
is half way between the two values onces converted to unsigned 64-bit
integers. This means no more than 64 steps will be taken, fewer if `a`
and `b` already share some bits.

The process is guaranteed to return a value `c` with `f(c)` one of
`0`, `Inf`, or `NaN`; *or* one of `f(prevfloat(c))*f(c) < 0` or
`f(c)*f(nextfloat(c)) > 0` holding.

This function is a bit faster than the slightly more general
`find_zero(f, [a,b], Bisection())` call.

Due to Jason Merrill. Customized by Patrick Moran to allow for a tolerance

"""
function bisection64_custom(f, a::Float64, b::Float64, tol::Float64)

    if a > b
        b,a = a, b
    end

    m = _middle(a,b)
    fa, fb = sign(f(a)), sign(f(b))

    fa * fb > 0 && throw(ArgumentError(bracketing_error))
    (iszero(fa) || isnan(fa) || isinf(fa)) && return a
    (iszero(fb) || isnan(fb) || isinf(fb)) && return b

    while a < m < b
        f_val = f(m)
        fm = sign(f_val)
        # println("m = $m and f_val = $f_val")

        if (abs(f_val) < tol) || iszero(fm) || isnan(fm) || isinf(fm)
            return m
        elseif fa * fm < 0
            b,fb=m,fm
        else
            a,fa=m,fm
        end
        m = _middle(a,b)
    end
    return m
end

# I pass in the args to f for speed, as suggested here:
# https://mmas.github.io/bisection-method-julia
function bisection64_with_args(f, a::Float64, b::Float64, tol::Float64, args=()::Tuple)

    if a > b
        b,a = a, b
    end

    m = _middle(a,b)
    fa, fb = sign(f(a, args...)), sign(f(b, args...))

    fa * fb > 0 && throw(ArgumentError(bracketing_error))
    (iszero(fa) || isnan(fa) || isinf(fa)) && return a
    (iszero(fb) || isnan(fb) || isinf(fb)) && return b

    while a < m < b
        f_val = f(m, args...)
        fm = sign(f_val)
        # println("m = $m and f_val = $f_val")

        if (abs(f_val) < tol) || iszero(fm) || isnan(fm) || isinf(fm)
            return m
        elseif fa * fm < 0
            b,fb=m,fm
        else
            a,fa=m,fm
        end
        m = _middle(a,b)
    end
    return m
end

################################################################################
#
# # https://mmas.github.io/bisection-method-julia
# function bisection_test(f::Function, a::Number, b::Number;
#                    tol::AbstractFloat=1e-5, maxiter::Integer=100)
#     fa = f(a)
#     fa*f(b) <= 0 || error("No real root in [a,b]")
#     i = 0
#     local c
#     while b-a > tol
#         i += 1
#         i != maxiter || error("Max iteration exceeded")
#         c = (a+b)/2
#         fc = f(c)
#         if fc == 0
#             break
#         elseif fa*fc > 0
#             a = c  # Root is in the right half of [a,b].
#             fa = fc
#         else
#             b = c  # Root is in the left half of [a,b].
#         end
#     end
#     return c
# end
