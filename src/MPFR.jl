module MPFR

export
    MPFRFloat,
    prec,
    get_default_precision,
    set_default_precision,
    with_precision
    
import
    Base.(*),
    Base.+,
    Base.-,
    Base./,
    Base.<,
    Base.<<,
    Base.>>,
    Base.<=,
    Base.==,
    Base.>,
    Base.>=,
    Base.^,
    Base.(~),
    Base.(&),
    Base.(|),
    Base.($),
    Base.binomial,
    Base.ceil,
    Base.cmp,
    Base.complex,
    Base.convert,
    Base.copysign,
    Base.div,
    Base.exponent,
    Base.factorial,
    Base.fld,
    Base.floor,
    Base.gcd,
    Base.gcdx,
    Base.integer_valued,
    Base.isfinite,
    Base.isinf,
    Base.isnan,
    Base.lcm,
    Base.max,
    Base.min,
    Base.mod,
    Base.modf,
    Base.ndigits,
    Base.nextfloat,
    Base.prevfloat,
    Base.promote_rule,
    Base.rem,
    Base.round,
    Base.show,
    Base.showcompact,
    Base.sqrt,
    Base.string,
    Base.trunc

const ROUNDING_MODE = 0
const DEFAULT_PRECISION = [53]

# Basic type and initialization definitions

immutable MPFRFloat{N} <: FloatingPoint
    mpfr::Vector{Int32}
    # function MPFRFloat()
    #     z = Array(Int32, 5)
    #     ccall((:mpfr_init,:libmpfr), Void, (Ptr{Void}, ), z)
    #     b = new(z)
    #     finalizer(b.mpfr, MPFR_clear)
    #     return b
    # end
    function MPFRFloat()
        if N < 2
            error("Invalid precision")
        end
        z = Array(Int32, 5)
        ccall((:mpfr_init2,:libmpfr), Void, (Ptr{Void}, Int), z, N)
        b = new(z)
        finalizer(b.mpfr, MPFR_clear)
        return b
    end
end

MPFR_clear(mpfr::Vector{Int32}) = ccall((:mpfr_clear, :libmpfr), Void, (Ptr{Void},), mpfr)

function MPFRFloat{N}(x::MPFRFloat{N})
    z = MPFRFloat{N}()
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    return z
end

for (fJ, fC) in ((:si,:Int), (:ui,:Uint), (:d,:Float64))
    @eval begin
        function MPFRFloat(x::($fC))
            z = MPFRFloat{DEFAULT_PRECISION[end]}()
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{Void}, ($fC), Int32), z.mpfr, x, ROUNDING_MODE)   
            return z
        end
    end
end

function MPFRFloat(x::BigInt)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpz, ROUNDING_MODE)   
    return z
end

function MPFRFloat(x::BigFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_set_f, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpf, ROUNDING_MODE)   
    return z
end

function MPFRFloat(x::String, base::Int)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{Void}, Ptr{Uint8}, Int32, Int32), z.mpfr, x, base, ROUNDING_MODE)
    if err != 0; error("Invalid input"); end
    return z
end
MPFRFloat(x::String) = MPFRFloat(x, 10)


MPFRFloat(x::Bool) = MPFRFloat(uint(x))
MPFRFloat(x::Signed) = MPFRFloat(int(x))
MPFRFloat(x::Unsigned) = MPFRFloat(uint(x))
if WORD_SIZE == 32
    MPFRFloat(x::Int64) = MPFRFloat(string(x))
    MPFRFloat(x::Uint64) = MPFRFloat(BigInt(x))
end
MPFRFloat(x::Float32) = MPFRFloat(float64(x))
MPFRFloat(x::Rational) = MPFRFloat(num(x)) / MPFRFloat(den(x))

# TODO: fix the precision support here
convert{N}(::Type{MPFRFloat{N}}, x::Rational) = MPFRFloat(x) # to resolve ambiguity
convert{N}(::Type{MPFRFloat{N}}, x::Real) = MPFRFloat(x)
convert(::Type{MPFRFloat}, x::Real) = MPFRFloat(x)

convert(::Type{Float64}, x::MPFRFloat) = ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{Void},), x.mpfr)
convert(::Type{Float32}, x::MPFRFloat) = ccall((:mpfr_get_flt,:libmpfr), Float32, (Ptr{Void},), x.mpfr)

promote_rule{T<:Real,N}(::Type{MPFRFloat{N}}, ::Type{T}) = MPFRFloat{N}
promote_rule{T<:Real}(::Type{MPFRFloat}, ::Type{T}) = MPFRFloat

# TODO: Decide if overwriting the default BigFloat rule is good
#promote_rule{T<:FloatingPoint}(::Type{BigInt},::Type{T}) = MPFRFloat{DEFAULT_PRECISION[end]}
#promote_rule{T<:FloatingPoint}(::Type{BigFloat},::Type{T}) = MPFRFloat{DEFAULT_PRECISION[end]}

# Basic operations

for (fJ, fC) in ((:+,:add), (:-,:sub), (:*,:mul), (:/,:div), (:^, :pow))
    @eval begin 
        function ($fJ)(x::MPFRFloat, y::MPFRFloat)
            z = MPFRFloat{DEFAULT_PRECISION[end]}()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
            return z
        end
    end
end

function -(x::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_neg, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    return z
end

function cmp(x::MPFRFloat, y::MPFRFloat)
    ccall((:mpfr_cmp, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr)
end

function sqrt(x::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_sqrt, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    if isnan(z)
        throw(DomainError())
    end
    return z
end

for f in (:ceil, :floor, :trunc)
    @eval begin
        function ($f)(x::MPFRFloat)
            z = MPFRFloat{DEFAULT_PRECISION[end]}()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), z.mpfr, x.mpfr)
            return z
        end
    end
end

function ^(x::MPFRFloat, y::Uint)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Uint, Int32), z.mpfr, x.mpfr, y, ROUNDING_MODE)
    return z
end

function ^(x::MPFRFloat, y::Int)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int, Int32), z.mpfr, x.mpfr, y, ROUNDING_MODE)
    return z
end

function ^(x::MPFRFloat, y::BigInt)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpz, ROUNDING_MODE)
    return z
end

function max(x::MPFRFloat, y::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_max, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
    return z
end

function min(x::MPFRFloat, y::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_min, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
    return z
end

function modf(x::MPFRFloat)
    if isinf(x)
        return (MPFRFloat(NaN), x)
    end
    zint = MPFRFloat{DEFAULT_PRECISION[end]}()
    zfloat = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_modf, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), zint.mpfr, zfloat.mpfr, x.mpfr, ROUNDING_MODE)
    return (zfloat, zint)
end

function rem(x::MPFRFloat, y::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_remainder, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
    return z
end

# Utility functions
==(x::MPFRFloat, y::MPFRFloat) = ccall((:mpfr_equal_p, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr) != 0
<=(x::MPFRFloat, y::MPFRFloat) = ccall((:mpfr_lessequal_p, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr) != 0
>=(x::MPFRFloat, y::MPFRFloat) = ccall((:mpfr_greaterequal_p, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr) != 0
<(x::MPFRFloat, y::MPFRFloat) = ccall((:mpfr_less_p, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr) != 0
>(x::MPFRFloat, y::MPFRFloat) = ccall((:mpfr_greater_p, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr) != 0

function prec(x::MPFRFloat)
    return ccall((:mpfr_get_prec, :libmpfr), Int, (Ptr{Void},), x.mpfr)
end

get_default_precision() = DEFAULT_PRECISION[end]
function set_default_precision(x::Int)
    if x < 2
        throw(DomainError())
    end
    DEFAULT_PRECISION[end] = x
end

function copysign(x::MPFRFloat, y::MPFRFloat)
    z = MPFRFloat{DEFAULT_PRECISION[end]}()
    ccall((:mpfr_copysign, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
    return z
end

function exponent(x::MPFRFloat)
    if x == 0 || !isfinite(x)
        throw(DomainError())
    end
    # The '- 1' is to make it work as Base.exponent
    return ccall((:mpfr_get_exp, :libmpfr), Int, (Ptr{Void},), x.mpfr) - 1
end

function integer_valued(x::MPFRFloat)
    return ccall((:mpfr_integer_p, :libmpfr), Int32, (Ptr{Void},), x.mpfr) != 0
end

isfinite(x::MPFRFloat) = !isinf(x)
function isinf(x::MPFRFloat)
    return ccall((:mpfr_inf_p, :libmpfr), Int32, (Ptr{Void},), x.mpfr) != 0
end

function isnan(x::MPFRFloat)
    return ccall((:mpfr_nan_p, :libmpfr), Int32, (Ptr{Void},), x.mpfr) != 0
end

function nextfloat(x::MPFRFloat)
    z = MPFRFloat(x)
   ccall((:mpfr_nextabove, :libmpfr), Int32, (Ptr{Void},), z.mpfr) != 0
   return z
end

function prevfloat(x::MPFRFloat)
    z = MPFRFloat(x)
   ccall((:mpfr_nextbelow, :libmpfr), Int32, (Ptr{Void},), z.mpfr) != 0
   return z
end

function with_precision(f::Function, precision::Integer)
    old_precision = get_default_precision()
    set_default_precision(precision)
    ret = f()
    set_default_precision(old_precision)
    return ret
end

# WARNING: it rounds to prec bits, and not to prec digits.
function round(x::MPFRFloat, prec::Int)
    if prec < 2
        throw(DomainError())
    end
    z = MPFRFloat(x)
    ccall((:mpfr_prec_round, :libmpfr), Int32, (Ptr{Void}, Int, Int32), z.mpfr, prec, ROUNDING_MODE)
    return z
end

function string(x::MPFRFloat)
    lng = 128
    for i = 1:2
        z = Array(Uint8, lng)
        lng = ccall((:mpfr_snprintf,:libmpfr), Int32, (Ptr{Uint8}, Uint, Ptr{Uint8}, Ptr{Void}...), z, lng, "%.Re", x.mpfr)
        if lng < 128 || i == 2
            return bytestring(convert(Ptr{Uint8}, z[1:lng]))
        end
    end
end

show(io::IO, b::MPFRFloat) = print(io, string(b) * " with $(prec(b)) bits of precision")
showcompact(io::IO, b::MPFRFloat) = print(io, string(b))

end #module