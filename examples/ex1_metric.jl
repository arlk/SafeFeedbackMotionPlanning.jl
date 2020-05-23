using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra

λ=1.74
βupp=5.0
βlow=3.0

### Setup
@polyvar x[1:2]
B = [0.5; -2]
Bp = [2; 0.5]

# Singh ICRA Example
dx1 =   -x[1] + 2*x[2]
dx2 = -3*x[1] + 4*x[2] -0.25*x[2]^3
f = [dx1, dx2]
fx = differentiate(f, x)

##### SOS
model = SOSModel(Mosek.Optimizer)
@variable(model, W[1:2, 1:2])
@constraint(model, βupp*I-W in PSDCone())
@constraint(model, W-βlow*I in PSDCone())
ccm=-Bp'*(fx*W+W*fx'+2*λ*W)*Bp;
@constraint(model, ccm in SOSCone())

optimize!(model)
Wsol = round.(value.(W); digits=7)
