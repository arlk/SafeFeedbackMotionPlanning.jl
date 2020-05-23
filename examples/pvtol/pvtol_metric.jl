using SumOfSquares
using DynamicPolynomials
using MosekTools
using LinearAlgebra
using JLD

λ=0.8
βupp=160
βlow=1
@polyvar x[1:6] # px, pz, ϕ, vx, vz, dϕ̇
m = 0.486 # kg
J = 0.00383 # kg m²
l = 0.25 # m
g = 9.81 # m/s²
B = [zeros(4,2); 1/m 1/m; l/J -l/J]
Bp = [I; zeros(2,4)]

# some sanity checks
@assert norm(Bp'*B) == 0
@assert λ > 0

# trignometric functions as polynomials: from Stanford ASL's RobustMP code
# replace later if works
sinp(x) = 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
cosp(x) = 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

# Singh IJRR example: planar quadrotor
dx1 =  x[4]*cosp(x[3]) - x[5]*sinp(x[3])
dx2 =  x[4]*sinp(x[3]) + x[5]*cosp(x[3])
dx3 =  x[6]
dx4 =  x[5]*x[6] - g*sinp(x[3])
dx5 = -x[4]*x[6] - g*cosp(x[3])
dx6 =  0
f = [dx1, dx2, dx3, dx4, dx5, dx6]
fx = differentiate(f, x)

##### SOS
model = SOSModel(Mosek.Optimizer)
@variable(model, W[1:6,1:6], Poly(monomials(x[3:4], 0:3)))
#  @variable(model, Whot[1:4,1:4], Poly(monomials(x[3:4], 3:4)))
#  W = W + [Whot zeros(4,2); zeros(2,6)]
W = W + W' # enforce symmetricity

Wset = @set x[3]^2 ≤ (π/2)^2 && x[4]^2 ≤ 2^2
upper = βupp*I - W
lower = W - βlow*I
@constraint(model, upper in PSDCone(), domain = Wset)
@constraint(model, lower in PSDCone(), domain = Wset)

CCMset = @set Wset && x[6]^2 ≤ (π/3)^2 && x[5]^2 ≤ 1
dWf = differentiate.(W, x[3])*dx3 + differentiate.(W, x[4])*dx4
ccm=-Bp'*(-dWf + fx*W + W*fx' + 2*λ*W)*Bp;
@constraint(model, ccm in PSDCone(), domain = CCMset)

optimize!(model)
@show termination_status(model)
@show βupp βlow λ
Wsol = round.(value.(W); digits=7)

Wsolx3 = differentiate.(Wsol, x[3])
Wsolx32 = differentiate.(Wsolx3, x[3])/2
Wsolx33 = differentiate.(Wsolx32, x[3])/3
Wsolx32x4 = differentiate.(Wsolx32, x[4])
Wsolx4 = differentiate.(Wsol, x[4])
Wsolx42 = differentiate.(Wsolx4, x[4])/2
Wsolx43 = differentiate.(Wsolx42, x[4])/3
Wsolx3x4 = differentiate.(Wsolx3, x[4])
Wsolx3x42 = differentiate.(Wsolx3x4, x[4])/2

subsMat(wij, val) = wij(x[3]=>val,x[4]=>val)

Wx33 = subsMat.(Wsolx33,0)
Wx32x4 = subsMat.(Wsolx32x4,0)
Wx3x42 = subsMat.(Wsolx3x42,0)
Wx43 = subsMat.(Wsolx43,0)
Wx32 = subsMat.(Wsolx32, 0)
Wx3x4 = subsMat.(Wsolx3x4, 0)
Wx42 = subsMat.(Wsolx42, 0)
Wx3 = subsMat.(Wsolx3, 0)
Wx4 = subsMat.(Wsolx4, 0)
W = subsMat.(Wsol, 0)

@save "pvtol_W.jld" Wx33 Wx32x4 Wx3x42 Wx43 Wx32 Wx3x4 Wx42 Wx3 Wx4 W
