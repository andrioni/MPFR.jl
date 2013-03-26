module MPFR

export
    MPFRFloat,
    prec,
    set_prec
    
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
    Base.div,
    Base.factorial,
    Base.fld,
    Base.floor,
    Base.gcd,
    Base.gcdx,
    Base.isinf,
    Base.isnan,
    Base.lcm,
    Base.mod,
    Base.ndigits,
    Base.promote_rule,
    Base.rem,
    Base.round,
    Base.show,
    Base.showcompact,
    Base.sqrt,
    Base.string,
    Base.trunc

const ROUNDING_MODE = 0
const DEFAULT_PRECISION = 53

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

function MPFRFloat(x::MPFRFloat)
    z = MPFRFloat(prec(x))
    ccall((:mpfr_set, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    return z
end
MPFRFloat(x::MPFRFloat, N::Int) = round(x, N)

for (fJ, fC) in ((:si,:Int), (:ui,:Uint), (:d,:Float64))
    @eval begin
        function MPFRFloat(x::($fC), N::Int)
            z = MPFRFloat{N}()
            ccall(($(string(:mpfr_set_,fJ)), :libmpfr), Int32, (Ptr{Void}, ($fC), Int32), z.mpfr, x, ROUNDING_MODE)   
            return z
        end
        MPFRFloat(x::($fC)) = MPFRFloat(x, DEFAULT_PRECISION)
    end
end

function MPFRFloat(x::BigInt, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_set_z, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpz, ROUNDING_MODE)   
    return z
end
MPFRFloat(x::BigInt) = MPFRFloat(x, DEFAULT_PRECISION)

function MPFRFloat(x::BigFloat, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_set_f, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpf, ROUNDING_MODE)   
    return z
end
MPFRFloat(x::BigFloat) = MPFRFloat(x, DEFAULT_PRECISION)

function MPFRFloat(x::String, base::Int, N::Int)
    z = MPFRFloat{N}()
    err = ccall((:mpfr_set_str, :libmpfr), Int32, (Ptr{Void}, Ptr{Uint8}, Int32, Int32), z.mpfr, x, base, ROUNDING_MODE)
    if err != 0; error("Invalid input"); end
    return z
end
MPFRFloat(x::String, base::Int) = MPFRFloat(x, base, DEFAULT_PRECISION)
MPFRFloat(x::String) = MPFRFloat(x, 10, DEFAULT_PRECISION)


MPFRFloat(x::Bool) = MPFRFloat(uint(x))
MPFRFloat(x::Signed) = MPFRFloat(int(x))
MPFRFloat(x::Unsigned) = MPFRFloat(uint(x))
if WORD_SIZE == 32
    MPFRFloat(x::Int64) = MPFRFloat(string(x))
    MPFRFloat(x::Uint64) = MPFRFloat(BigInt(x))
end
MPFRFloat(x::Float32) = MPFRFloat(float64(x))
MPFRFloat(x::Rational) = MPFRFloat(num(x)) / MPFRFloat(den(x))

convert(::Type{MPFRFloat}, x::Rational) = MPFRFloat(x) # to resolve ambiguity
convert(::Type{MPFRFloat}, x::Real) = MPFRFloat(x)

convert(::Type{Float64}, x::MPFRFloat) = ccall((:mpfr_get_d,:libmpfr), Float64, (Ptr{Void},), x.mpfr)
#convert(::Type{FloatingPoint}, x::BigInt) = MPFRFloat(x)

promote_rule{T<:Union(Integer,FloatingPoint)}(::Type{MPFRFloat}, ::Type{T}) = MPFRFloat
#promote_rule{T<:FloatingPoint}(::Type{BigInt},::Type{T}) = MPFRFloat
promote_rule{T<:FloatingPoint}(::Type{BigFloat},::Type{T}) = MPFRFloat

# Basic operations

for (fJ, fC) in ((:+,:add), (:-,:sub), (:*,:mul), (:/,:div), (:^, :pow))
    @eval begin 
        function ($fJ)(x::MPFRFloat, y::MPFRFloat, N::Int)
            z = MPFRFloat{N}()
            ccall(($(string(:mpfr_,fC)),:libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpfr, ROUNDING_MODE)
            return z
        end
        ($fJ)(x::MPFRFloat, y::MPFRFloat) = ($fJ)(x, y, max(prec(x),prec(y)))
    end
end

function -(x::MPFRFloat, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_neg, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    return z
end
-(x::MPFRFloat) = -(x::MPFRFloat, prec(x))

function cmp(x::MPFRFloat, y::MPFRFloat)
    ccall((:mpfr_cmp, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), x.mpfr, y.mpfr)
end

function sqrt(x::MPFRFloat, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_sqrt, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, ROUNDING_MODE)
    return z
end
sqrt(x::MPFRFloat) = sqrt(x, prec(x))

for f in (:ceil, :floor, :trunc)
    @eval begin
        function ($f)(x::MPFRFloat, N::Int)
            z = MPFRFloat{Int}()
            ccall(($(string(:mpfr_,f)), :libmpfr), Int32, (Ptr{Void}, Ptr{Void}), z.mpfr, x.mpfr)
            return z
        end
        ($f)(x::MPFRFloat) = ($f)(x, prec(x))
    end
end

function ^(x::MPFRFloat, y::Uint, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_pow_ui, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Uint, Int32), z.mpfr, x.mpfr, y, ROUNDING_MODE)
    return z
end
^(x::MPFRFloat, y::Uint) = ^(x, y, prec(x))

function ^(x::MPFRFloat, y::Int, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_pow_si, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Int, Int32), z.mpfr, x.mpfr, y, ROUNDING_MODE)
    return z
end
^(x::MPFRFloat, y::Int) = ^(x, y, prec(x))

function ^(x::MPFRFloat, y::BigInt, N::Int)
    z = MPFRFloat{N}()
    ccall((:mpfr_pow_z, :libmpfr), Int32, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Int32), z.mpfr, x.mpfr, y.mpz, ROUNDING_MODE)
    return z
end
^(x::MPFRFloat, y::BigInt) = ^(x, y, prec(x))

==(x::MPFRFloat, y::MPFRFloat) = cmp(x,y) == 0
<=(x::MPFRFloat, y::MPFRFloat) = cmp(x,y) <= 0
>=(x::MPFRFloat, y::MPFRFloat) = cmp(x,y) >= 0
<(x::MPFRFloat, y::MPFRFloat) = cmp(x,y) < 0
>(x::MPFRFloat, y::MPFRFloat) = cmp(x,y) > 0

# Utility functions

function prec(x::MPFRFloat)
    return ccall((:mpfr_get_prec, :libmpfr), Int, (Ptr{Void},), x.mpfr)
end

# WARNING: it rounds to prec bits, and not to prec digits.
function round(x::MPFRFloat, prec::Int)
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