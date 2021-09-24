module FlightGNC
using LinearAlgebra

using Reexport
@reexport using FSimBase
const FSBase = FSimBase
import FSimBase: State, Params, Dynamics!, Dynamics
import FSimBase: Command
using UnPack, ComponentArrays
using Plots  # required for fig_print

# exports -------------------------------------------------------------
## algorithms
export BPNG, BPNG_cmd
export Bias_zero, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D

## environments
export PointMassMissile, PursuerEvadorMissile

## utils
export fig_print

# includes -------------------------------------------------------------
## algorithms
include("algorithms/BPNG.jl")

## environments
include("environments/missiles.jl")

## utils
include("utils/fig_print.jl")
end
