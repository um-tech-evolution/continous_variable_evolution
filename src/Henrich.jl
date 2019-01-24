# Includes all Julia code except for the file  src/run.jl.
module ContVarEvolution
#=
try # for julia v7
  using Random
  using Distributed
catch
end
try 
  using Dates
catch
end
=#
include("htypes.jl")     
include("propsel.jl")   # Proportional selection
include("power_sel.jl")   # Proportional selection applied to a power of fitness
include("henrich.jl")  # Primary code for simulation
include("hrun_cont_var.jl")
include("bin_data.jl")
end
using Main.ContVarEvolution
