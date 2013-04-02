using Test
using MPFR

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

