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
