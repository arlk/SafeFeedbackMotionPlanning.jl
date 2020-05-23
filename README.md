# SafeFeedbackMotionPlanning

*This package is under active development*

Design safe controllers for nonlinear systems in the presence of disturbances. The controllers ensure that the system state remains in a tunable *safe tube* around the desired state. Find out more about it in our [paper](https://github.com/arlk/SafeFeedbackMotionPlanning/blob/master/arxiv.pdf):
```
@article{lakshmanan2020safe,
  title={Safe Feedback Motion Planning: A Contraction Theory and $\mathcal{L}_1$-Adaptive Control Based Approach},
  author={Lakshmanan, Arun and Gahlawat, Aditya and Hovakimyan, Naira},
  journal={arXiv preprint arXiv:2004.01142},
  year={2020}
}
```
## Problem Formulation
Consider a system of the form
```julia
ẋ = f(x) + B(x)u
```
where `x(t), f(x) ∈ ℝⁿ`, `u(t) ∈ ℝᵐ`, and `B(x) ∈ ℝⁿˣᵐ`, that is facing disturbances `h(t,x)  ∈ ℝᵐ` (either/both state and time dependent) that is matched with the control channel:
```julia
ẋ = f(x) + B(x)(u + h(t,x))
```
Given a desired state-control trajectory pair `(x*,u*)` (feasible under the nominal system) that the actual system is required to follow, design a controller `u(t)` such that the states remain inside a tube `Ω(ρ, x*(t)) = {ℝⁿ : ||y - x*(t)|| ≤ ρ}` of some user defined width `ρ > 0`.

In this package you can design controllers to do exactly this! A more in-depth review of the control architecture can be found in the paper.

## Usage

### Installation

```julia
julia> ] add https://github.com/arlk/SafeFeedbackMotionPlanning.jl
```

### Example

```julia
using SafeFeedbackMotionPlanning
using DifferentialEquations
using Plots
using LinearAlgebra

# Define system matrices or functions
f(x) = [-x[1] + 2*x[2];
        -0.25*x[2]^3 - 3*x[1] + 4*x[2]]
B = [0.5, -2.0]

# Simulation time span
tspan = (0., 5.)

# Intial condition of the system
x0 = [1.0, -1.0]

# Define desired regulation point or trajectories (as functions)
xs = [0.0, 0.0]
us = 0.0

# Control Contraction Metric
W = [4.25828 -0.93423; -0.93423 3.76692]
λ = 1.74

# Uncertainty unknown to the controller
h(t, x) = -2*sin(2*t) -0.01*x'*x

# Upper bound on the norm of the uncertainty
# for any x ∈ Ω(ρ) and t ≥ 0 (see paper for more details)
Δh = 2.1
ω = 90
Γ = 4e7

# Construct objects using the parameters you just defined
sys_p = sys_params(f, B)
ccm_p = ccm_params(xs, us, λ, W)
l1_p = l1_params(ω, Γ, Δh)

# CCM only system without perturbations
nom_sys = nominal_system(sys_p, ccm_p)
nom_sol = solve(nom_sys, x0, tspan, Tsit5(), progress = true, progress_steps = 1)

# CCM only system with perturbations
ptb_sys = nominal_system(sys_p, ccm_p, h)
ptb_sol = solve(ptb_sys, x0, tspan, Tsit5(), progress = true, progress_steps = 1)

# L1 Reference sysem with perturbations AKA (non-implementable) ideal
# performance of the L1 system with full knowledge of the uncertainty
ref_sys = reference_system(sys_p, ccm_p, ω, h)
ref_sol = solve(ref_sys, x0, tspan, Tsit5(), progress = true, progress_steps = 1)

# L1 + CCM system with perturbations
l1_sys = l1_system(sys_p, ccm_p, l1_p, h)
l1_sol = solve(l1_sys, x0, tspan, Rosenbrock23(), progress = true, progress_steps = 1, saveat = 0.01)

l = @layout [a b; c d]
p1 = plot(nom_sol, vars=[(0,1),(0,2)], title="Ideal CCM performance", legend=:none)
p2 = plot(ptb_sol, vars=[(0,1),(0,2)], title="CCM with perturbations", legend=:none)
p3 = plot(ref_sol, vars=[(0,1),(0,2)], title="Ideal L1+CCM with perturbations\n (non-implementable)", legend=:none)
p4 = plot(l1_sol, vars=[(0,1),(0,2)], title="L1+CCM with perturbations", legend=:none)
plot(p1, p2, p3, p4, layout=l)
```

![](https://github.com/arlk/SafeFeedbackMotionPlanning.jl/raw/master/examples/readme_ex.png)

## TODO
 - CCM metric search using SumofSquares.jl
 - Automatically computing ω and Γ from the tube bounds
