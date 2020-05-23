using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra

X=0.1
λ=1.0
βupp=0.26
βlow=0.17

### Setup
@polyvar x[1:3]
B = [0; 0; 1]
Bp = [1 0; 0 1; 0 0]

# Slotine TAC, Leung ACC example
dx1 = -x[1] + x[3]
dx2 =  x[1]^2 - x[2] - 2*x[1]*x[3] + x[3]
dx3 = -x[2]
f = [dx1, dx2, dx3]
fx = differentiate(f, x)

model = SOSModel(optimizer_with_attributes(Mosek.Optimizer, "LOG"=>0))
@variable(model, W[1:3,1:3], Poly(monomials(x[1], 0:2)))
W = W + W' # enforce symmetricity
dWf = differentiate.(W, x[1])*dx1
ccm=-Bp'*(-dWf + fx*W + W*fx' + 2*λ*W)*Bp;

Wset = @set x[1]^2 ≤ X^2
upper = βupp*I - W
lower = W - βlow*I
@constraint(model, upper in PSDCone(), domain = Wset)
@constraint(model, lower in PSDCone(), domain = Wset)
CCMset = @set Wset && x[2]^2 ≤ X^2 && x[3]^2 ≤ X^2
@constraint(model, ccm in PSDCone(), domain = CCMset)
optimize!(model)

Wsol = round.(value.(W); digits=5)

Wsolx = differentiate.(Wsol, x[1])
Wsolxx = differentiate.(Wsolx, x[1])

subsMat(wij, val) = wij(x[1]=>val)

Wc = subsMat.(Wsol, 0)
Wl = subsMat.(Wsolx, 0)
Wq = subsMat.(Wsolxx, 0)/2
