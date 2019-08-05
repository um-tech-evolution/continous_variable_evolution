# Includes all Julia code except for the file  src/run.jl.
module ContVarEvolution
#=
try 
  using Random
  using Distributed
catch
end
try 
  using Dates
catch
end
=#
include("types.jl")     
include("propsel.jl")   # Proportional selection
include("cont_var.jl")  # Primary code for simulation
include("cont_entropy.jl")  # Continuous entropy code
include("run_cont_var.jl")
include("bin_data.jl")
end
using Main.ContVarEvolution
