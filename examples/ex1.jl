using SafeFeedbackMotionPlanning
using DifferentialEquations
using Plots
using LinearAlgebra

# System
f(x) = [ -x[1] + x[3];
          x[1]^2 - 2*x[1]*x[3] - x[2] + x[3];
         -x[2]  ]
B = [0; 0; 1]

# Sim time
tspan = (0., 10.)

# Initial condition
x0 = [0.1, -0.1, 0.1]

# CCM
const Wc = [0.20365   0.00015  -0.00566;
            0.00015   0.21671  -0.00566;
           -0.00566  -0.00566   0.21607]
const Wl = [0.0      -0.4073   -0.00133;
           -0.4073   -0.00061   0.00998;
           -0.00133   0.00998   0.06785]
const Wq = [0.0      0.0      0.00053;
            0.0      0.81459  0.00308;
            0.00053  0.00308  0.00343]
Wf(x) = Wc .+ Wl.*x[1] .+ Wq.*(x[1]^2)

# CCM convergence rate
λ = 1.0

# Desired regulation point
xs = [0.0, 0.0, 0.0]
us = 0.0

# Uncertainty and upper bound
h(t, x) = 0.1*sin(2*t)
Δh = 0.1

# L1 filter bandwidth and adaptation rate
ω = 50
Γ = 5e6


sys_p = sys_params(f, B)
ccm_p = ccm_params(xs, us, λ, Wf)
l1_p = l1_params(ω, Γ, Δh)

# CCM only system with perturbation
ptb_sys = nominal_system(sys_p, ccm_p, h)
ptb_sol = solve(ptb_sys, x0, tspan, Tsit5(), progress = true, progress_steps = 1)

# L1 + CCM system to handle disturbances
l1_sys = l1_system(sys_p, ccm_p, l1_p, h)
l1_sol = solve(l1_sys, x0, tspan, Rosenbrock23(), progress = true, progress_steps = 1, saveat = 0.1)

plot(ptb_sol, vars=[(0,1)], linecolor=:red, linestyle=:dash, label="x1(t): CCM")
plot!(ptb_sol, vars=[(0,2)], linecolor=:blue, linestyle=:dash, label="x2(t): CCM")
plot!(ptb_sol, vars=[(0,3)], linecolor=:green, linestyle=:dash, label="x3(t): CCM")
plot!(l1_sol, vars=[(0,1)], linecolor=:red, label="x1(t): CCM+L1")
plot!(l1_sol, vars=[(0,2)], linecolor=:blue, label="x2(t): CCM+L1")
plot!(l1_sol, vars=[(0,3)], linecolor=:green, label="x3(t): CCM+L1")
