module SafeFeedbackMotionPlanning

using ComponentArrays
using ForwardDiff
using LinearAlgebra
using Optim
using OrdinaryDiffEq
using UnPack

import DiffEqBase: solve
import LineSearches: BackTracking

export sys_params, ccm_params, l1_params
export nominal_system, reference_system, l1_system

include("chebyshev.jl")
include("types.jl")
include("sys.jl")
include("ccm.jl")
include("l1.jl")

end # module
