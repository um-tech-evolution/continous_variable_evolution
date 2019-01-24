# Julia program to be called from the command line to run continous variable simulation
# Remaining code is in the module ContVarEvolution.
# Example run:  julia run.jl examples/example1
# Example run:  julia run.jl examples/example1  1     # runs with random seed = 1
# Example run:  julia -p 4 run.jl examples/example1
using DataStructures
try   # These are needed julia v7, but will fail in julia v6.  The try ... catch will recover from the failure
  using Distributed
  using Random
  using Dates
catch
end
@everywhere include("Henrich.jl")
#include("Henrich.jl")

@doc """function my_isdefined()
  S  should be a symbol that corresponds to a global variable.
  For example:  global A
  The corresponding symbol is :A
  my_isdefined(:A)  returns true if variable A has a value assigned to it, and false otherwise.
  A replacement for the  julia v6 function "isdefined()" which is deprecated in v7
  This definition works in both v6 and v7
"""
function my_isdefined( S::Symbol )
  try
    eval(S)
    return true
  catch
    return false
  end
end

@doc """ function save_params()
  Save the parameters specified in the configuration file into ContVarEvolution.cont_var_result() data structure.
  The parameters are global variables.
"""
function save_params() 
  # sim_record  is a record containing both the parameters and the results for a trial
  println("mutation_stddev_list: ",mutation_stddev_list)
  println("N_mut_list: ",N_mut_list)
  sim_record = ContVarEvolution.cont_var_result( N_list,mutation_stddev_list, N_mut_list,num_trials,
       num_subpops,num_attributes, ngens, burn_in, 
       mutation_stddev, mutation_bias, ideal, fit_power, renormalize, neutral )
  return sim_record
end


@doc """ function check_parameters()
  Check that all necessary parameters have been read from the parameter file.
  With julia v. 6.2 on CentOS, missing parameters can cause a segmentation fault.
"""
function check_parameters()
  param_list = [:num_trials, :N_list, :num_subpops, :num_attributes, :ngens, :burn_in, :ideal, :fit_power, :renormalize, :neutral]
  for p in param_list
    if !my_isdefined(p)
      error("The parameter $(String(p)) is not defined in the parameter file: $(simname).jl.")
    end
  end
end

if length(ARGS) == 0
  simname = "examples/hexample1"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
include("$(simname).jl")
#println("simname: ",simname)
println("simtype: ",simtype)
if simtype != 4
  error("simtype must be 4 for Henrich model.")
end
check_parameters()
sim_record = save_params()
run_trials( simname, sim_record )
