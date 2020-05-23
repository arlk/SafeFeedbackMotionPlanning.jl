using StaticArrays
using Interpolations

import TrajectoryOptimization
const TO = TrajectoryOptimization

struct SinghICRA <: TO.AbstractModel end

function simple_traj(model, x0, xf, tf, dt; Q = 0.5, R = 1.0)
    n, m = size(model)
    Qm = Q*Diagonal(@SVector ones(n))
    Qfm = Q*Diagonal(@SVector ones(n))
    Rm = R*Diagonal(@SVector ones(1))
    N = length(0:dt:tf)

    conSet = TO.ConstraintSet(n,m,N)
    goal = TO.GoalConstraint(xf)
    TO.add_constraint!(conSet, goal, N:N)
    obj = TO.LQRObjective(Qm, Rm, Qfm, xf, N)
    prob = TO.Problem(model, obj, xf, tf, constraints=conSet, x0=x0, integration=TO.RK3)
    solver = TO.iLQRSolver(prob)
    TO.solve!(solver)

    return solver
end

function desired_traj(solver)
    n, m = size(solver.model)
    t0 = range(0.0, solver.tf, length=solver.N)
    xs = Array(hcat(states(solver)...)')
    xitp = Interpolations.scale(interpolate(xs,
               (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t0, 1:n)
    xf(t) = xitp.(Ref(t),1:n)
    if m == 1
        us = Array(vcat(controls(solver)...))
        us = vcat(us, 0.0)
        uf = Interpolations.scale(interpolate(us,
                   BSpline(Cubic(Natural(OnGrid())))), t0)
    else
        us = Array(hcat(controls(solver)...)')
        us = vcat(us, zeros(1,m))
        uitp = Interpolations.scale(interpolate(us,
                   (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t0, 1:m)
        uf(t) = xitp.(Ref(t),1:m)
    end
    return xf, uf
end
