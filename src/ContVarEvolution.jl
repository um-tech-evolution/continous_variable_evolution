# Includes all Julia code except for the file  src/run.jl.
using Statistics
using Random
using Distributed
using Dates
module ContVarEvolution
include("types.jl")     
include("propsel.jl")   # Proportional selection
include("cont_var.jl")  # Primary code for simulation
include("run_cont_var.jl")
include("bin_data.jl")
end
using Main.ContVarEvolution
