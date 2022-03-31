module FlightGNC
using LinearAlgebra, Statistics

using Reexport
@reexport using FSimBase
const FSBase = FSimBase
import FSimBase: State, Params, Dynamics!, Dynamics, Command, apply_inputs, Simulator, solve, reinit!, log!, step!, step_until!

using UnPack, ComponentArrays
using Plots                             # required for fig_print
using Convex, Mosek, MosekTools         # required for IASCG
using Roots, OrdinaryDiffEq, Optim      # required for Speed_Predictor in IASCG
using DiffEqFlux                        # required for CTPG_train

# exports -------------------------------------------------------------
## algorithms
export BPNG, BPNG_cmd, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D, Bias_Quaternion_IACG, Bias_zero     # defined in BPNG.jl
export IASCG, IASCG_Planner, IASCG_cmd, Path_Generator, Speed_Predictor, Speed_Dynamics         # defined in IASCG.jl
export CTPG_train                                                                               # defined in CTPG.jl

## environments
export PointMassMissile, PursuerEvadorMissile
export VerticalPlaneDynamicMissile
export VerticalPlane_RigidBody_Missile

## utils
export fig_print, equal_AR_3D, view_result

# includes -------------------------------------------------------------
## algorithms
include("algorithms/BPNG.jl")
include("algorithms/IASCG.jl")
include("algorithms/CTPG.jl")

## environments
include("environments/missiles.jl")

## utils
include("utils/fig_utils.jl")
end
