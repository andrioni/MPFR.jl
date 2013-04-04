using Test
using MPFR

# constructors
x = MPFRFloat{53}()
x = MPFRFloat(12)
y = MPFRFloat(x)
@test x == y
y = MPFRFloat(0xc)
@test x == y
y = MPFRFloat(12.)
@test x == y
y = MPFRFloat(BigInt(12))
@test x == y
y = MPFRFloat(BigFloat(12))
@test x == y
y = MPFRFloat("12")
@test x == y
y = MPFRFloat(float32(12.))
@test x == y
y = MPFRFloat(12//1)
@test x == y

# basic constructor precision
x = MPFRFloat(12)
y = MPFRFloat(x; precision = 42)
@test prec(y) == 42
y = MPFRFloat(0xc; precision = 42)
@test prec(y) == 42
y = MPFRFloat(12.; precision = 42)
@test prec(y) == 42
y = MPFRFloat(BigInt(12); precision = 42)
@test prec(y) == 42
y = MPFRFloat(BigFloat(12); precision = 42)
@test prec(y) == 42
y = MPFRFloat("12"; precision = 42)
@test prec(y) == 42
y = MPFRFloat(float32(12.); precision = 42)
@test prec(y) == 42
y = MPFRFloat(12//1; precision = 42)
@test prec(y) == 42

# +
x = MPFRFloat(12)
y = MPFRFloat(30)
@test x + y == MPFRFloat(42)

# -
x = MPFRFloat(12)
y = MPFRFloat(-30)
@test x - y == MPFRFloat(42)

# *
x = MPFRFloat(6)
y = MPFRFloat(9)
@test x * y != MPFRFloat(42)
@test x * y == MPFRFloat(54)

# /
x = MPFRFloat(9)
y = MPFRFloat(6)
@test x / y == MPFRFloat(9/6)

# < / > / <= / >=
x = MPFRFloat(12)
y = MPFRFloat(42)
z = MPFRFloat(30)
@test y > x
@test y >= x
@test y > z
@test y >= z
@test x < y
@test x <= y
@test z < y
@test z <= y
@test y - x >= z
@test y - x <= z
@test !(x >= z)
@test !(y <= z)

# ^
x = MPFRFloat(12)
y = MPFRFloat(4)
@test x^y == MPFRFloat(20736)

# ceil
x = MPFRFloat(12.042)
@test MPFRFloat(13) == ceil(x)

# copysign
x = MPFRFloat(1)
y = MPFRFloat(-1)
@test copysign(x, y) == y
@test copysign(y, x) == x

# isfinite / isinf
x = MPFRFloat(Inf)
y = MPFRFloat(1)
@test isinf(x) == true
@test isinf(y) == false
@test isfinite(x) == false
@test isinf(x) == true

# isnan
x = MPFRFloat(NaN)
y = MPFRFloat(1)
@test isnan(x) == true
@test isnan(y) == false

# convert to
@test convert(MPFRFloat{53}, 1//2) == MPFRFloat("0.5")
@test convert(MPFRFloat{53}, 0.5) == MPFRFloat("0.5")
@test convert(MPFRFloat{53}, 40) == MPFRFloat("40")
@test convert(MPFRFloat{53}, float32(0.5)) == MPFRFloat("0.5")

# convert from
@test convert(Float64, MPFRFloat(0.5)) == 0.5
@test convert(Float32, MPFRFloat(0.5)) == float32(0.5)

# exponent
x = MPFRFloat(0)
@test_fails exponent(x)
x = MPFRFloat(Inf)
@test_fails exponent(x)
x = MPFRFloat(15.674)
@test exponent(x) == exponent(15.674)

# sqrt DomainError
@test_fails sqrt(MPFRFloat(-1))

# precision
old_precision = get_default_precision()
x = MPFRFloat(0)
@test prec(x) == old_precision
set_default_precision(256)
x = MPFRFloat(0)
@test prec(x) == 256
set_default_precision(old_precision)
z = with_precision(240) do
    z = x + 20
    return z
end
@test float(z) == 20.
@test prec(z) == 240
x = MPFRFloat(12)
@test prec(x) == old_precision
@test_fails set_default_precision(1)

# integer_valued
@test integer_valued(MPFRFloat(12))
@test !integer_valued(MPFRFloat(12.12))

# nextfloat / prevfloat
with_precision(53) do
    x = MPFRFloat(12.12)
    @test MPFRFloat(nextfloat(12.12)) == nextfloat(x)
    @test MPFRFloat(prevfloat(12.12)) == prevfloat(x)
end
@test isnan(nextfloat(MPFRFloat(NaN)))
@test isnan(prevfloat(MPFRFloat(NaN)))

# comparisons
x = MPFRFloat(1)
y = MPFRFloat(-1)
z = MPFRFloat(NaN)
ipl = MPFRFloat(Inf)
imi = MPFRFloat(-Inf)
@test x > y
@test x >= y
@test x >= x
@test y < x
@test y <= x
@test y <= y
@test x < ipl
@test x <= ipl
@test x > imi
@test x >= imi
@test imi == imi
@test ipl == ipl
@test imi < ipl
@test z != z
@test !(z == z)
@test !(z <= z)
@test !(z < z)
@test !(z >= z)
@test !(z > z)

# modf
x = MPFRFloat(12)
y = MPFRFloat(0.5)
@test modf(x+y) == (y, x)
x = MPFRFloat(NaN)
@test map(isnan, modf(x)) == (true, true)
x = MPFRFloat(Inf)
y = modf(x)
@test (isnan(y[1]), isinf(y[2])) == (true, true)

# rem
with_precision(53) do
    x = MPFRFloat(2)
    y = MPFRFloat(1.67)
    @test rem(x,y) == rem(2, 1.67)
    y = MPFRFloat(NaN)
    @test isnan(rem(x,y))
    @test isnan(rem(y,x))
    y = MPFRFloat(Inf)
    @test rem(x,y) == x
    @test isnan(rem(y,x))
end

# min/max
x = MPFRFloat(4)
y = MPFRFloat(2)
@test max(x,y) == x
@test min(x,y) == y
y = MPFRFloat(NaN)
@test max(x,y) == x
@test min(x,y) == x
@test isnan(max(y,y))
@test isnan(min(y,y))

