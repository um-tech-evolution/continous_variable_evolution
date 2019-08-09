# Julia program to be called from the command line to run continous variable simulation
# Remaining code is in the module ContVarEvolution.
# Example run:  julia run.jl examples/example1
# Example run:  julia run.jl examples/example1  1     # runs with random seed = 1
# Example run:  julia -p 4 run.jl examples/example1
using DataStructures
try   # These are needed julia v7, but will fail in julia v6.  The try ... catch will recover from the failure
  using Distributed
  using Random
  using Statistics
  using Dates
catch
end
push!(LOAD_PATH,"../src/")
@everywhere include("ContVarEvolution.jl")
#include("ContVarEvolution.jl")

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
  global mutation_stddev_list
  global N_mut_list
  global w
  int_burn_in = 0
  # sim_record  is a record containing both the parameters and the results for a trial
  if simtype == 2
    w = 0.0
  end
  if my_isdefined(:mutation_stddev_list)
    sim_record = ContVarEvolution.cont_var_result( N_list, num_attributes_list, mutation_stddev_list, Float64[], num_trials,N_list[1],
       num_subpops,num_attributes_list[1], ngens, burn_in, int_burn_in,
       mutation_stddev_list[1], ideal, fit_slope, neutral, w )
  elseif my_isdefined( :N_mut_list )
    sim_record = ContVarEvolution.cont_var_result( N_list, num_attributes_list, Float64[], N_mut_list, num_trials,N_list[1],
       num_subpops,num_attributes_list[1], ngens, burn_in, int_burn_in,
       N_mut_list[1]/N_list[1], ideal, fit_slope, neutral, w ) 
  else
    error("Either mutation_stddev_list or N_mut_list must be defined in the configuration file.")
  end
  return sim_record
end


@doc """ function check_parameters()
  Check that all necessary parameters have been read from the parameter file.
  With julia v. 6.2 on CentOS, missing parameters can cause a segmentation fault.
"""
function check_parameters()
  global num_trials, N_list, num_subpops, num_attributes_list, ngens, burn_in, mutation_stddev_list, N_mut_list, ideal, fit_slope, neutral
  param_list = [:num_trials, :N_list, :num_subpops, :num_attributes_list, :ngens, :burn_in, :ideal, :fit_slope, :neutral]
  for p in param_list
    if !my_isdefined(p)
      error("The parameter $(String(p)) is not defined in the parameter file: $(simname).jl.")
    end
  end
end

if length(ARGS) == 0
  simname = "examples/example2"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    Random.seed!(seed)
  end
end
include("$(simname).jl")
#println("simname: ",simname)
println("simtype: ",simtype)
#println("w: ",w)
check_parameters()
sim_record = save_params()
run_trials( simname, sim_record )
