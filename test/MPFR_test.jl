using Test
using MPFR

# copysign
x = MPFRFloat(1)
y = MPFRFloat(-1)
@test copysign(x, y) == y
@test copysign(y, x) == x
