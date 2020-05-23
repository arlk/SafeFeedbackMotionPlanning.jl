using SafeFeedbackMotionPlanning
using DifferentialEquations
using Plots
using LinearAlgebra

# System
f(x) = [-x[1] + 2*x[2];
        -0.25*x[2]^3 - 3*x[1] + 4*x[2]]
B = [0.5, -2.0]

# Sim time
tspan = (0., 5.)

# Initial condition
x0 = [3.4, -2.4]

# CCM
W = [4.25828 -0.93423; -0.93423 3.76692]
λ = 1.74

# Uncertainty and upper bound
h(t, x) = -2*sin(2*t) -0.01*x'*x
Δh = 2.1

# L1 filter bandwidth and adaptation rate
ω = 90
Γ = 4e7

############### Generating Feasible Desired Trajectory
include("traj.jl")
function TO.dynamics(model::SinghICRA, x, u)
    f(x) + B*u[1]
end
Base.size(::SinghICRA) = 2,1
xs0 = SA[3.4, -2.4]
xsf = SA[0.0, 0.0]
model = SinghICRA()
plan = simple_traj(model, xs0, xsf, tspan[2], 0.05)
xs, us = desired_traj(plan)
xy = Array(hcat(states(plan)...)')
######################################################

sys_p = sys_params(f, B)
ccm_p = ccm_params(xs, us, λ, W)
l1_p = l1_params(ω, Γ, Δh)

# CCM only system with perturbation
ptb_sys = nominal_system(sys_p, ccm_p, h)
ptb_sol = solve(ptb_sys, x0, tspan, Tsit5(), progress = true, progress_steps = 1)

# L1 + CCM system to handle disturbances
l1_sys = l1_system(sys_p, ccm_p, l1_p, h)
l1_sol = solve(l1_sys, x0, tspan, Rosenbrock23(), progress = true, progress_steps = 1, saveat = 0.01)

plot(xy[:,1], xy[:,2], linecolor=:black, linestyle=:dash, label="Desired")
plot!(ptb_sol, vars=[(1,2)], linecolor=:red, label="CCM")
plot!(l1_sol, vars=[(1,2)], linecolor=:blue, linewidth=2.0, label="L1+CCM")
