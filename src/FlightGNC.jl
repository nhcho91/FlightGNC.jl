module FlightGNC
using LinearAlgebra

using Reexport
@reexport using FSimBase
const FSBase = FSimBase
import FSimBase: State, Params, Dynamics!, Dynamics
import FSimBase: Command
using UnPack, ComponentArrays
using Plots  # required for fig_print
using Convex, Mosek, MosekTools
using Roots, OrdinaryDiffEq, Optim # required for Speed_Predictor in IASCG

# exports -------------------------------------------------------------
## algorithms
export BPNG, BPNG_cmd, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D, Bias_zero
export IASCG, IASCG_Planner, IASCG_cmd, Path_Generator, Speed_Predictor, Speed_Dynamics

## environments
export PointMassMissile, PursuerEvadorMissile
export VerticalPlaneDynamicMissile

## utils
export fig_print, equal_AR_3D

# includes -------------------------------------------------------------
## algorithms
include("algorithms/BPNG.jl")
include("algorithms/IASCG.jl")

## environments
include("environments/missiles.jl")

## utils
include("utils/fig_utils.jl")
end
