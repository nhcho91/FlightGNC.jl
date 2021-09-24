# using Revise
# includet("../main/main_missiles.jl")
include("../main/main_missiles.jl")

##
include("../main/sim_missiles.jl")

##
sim_plot(N, dim, α, χ_M_0, χ_f_d, σ_M_lim, s_Bias)

# Parallel execution
# Threads.nthreads()

# n   = 100 # number of scenarios
# Ns  = 1:n |> Map(i -> rand(1)[1]) |> collect
# Ns= 1:0.1:3
# sim_results = Ns |> Map(N -> main(N)) |> tcollect