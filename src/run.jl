# Julia program to be called from the command line to run continous variable simulation
# Remaining code is in the module ContVarEvolution.
# Example run:  julia -L ContVarEvolution.jl run.jl examples/example1
# Example run:  julia -L ContVarEvolution.jl run.jl examples/example1  1     # runs with random seed = 1
# Example run:  julia -p :4 -L ContVarEvolution.jl run.jl examples/example1
using ContVarEvolution

@doc """ function run_trials()
  Runs multiple trials of the simulation using the parameter file "\$(simname).jl" where  simname is the first command line argument.
  Trials can be run in parallel using the Julia parallel map (pmap) facility.
"""
function run_trials( simname::AbstractString ) 
  global mutation_stddev_list
  global N_mut_list
  #circular_variation = extreme_variation = false
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  # sim_record  is a record containing both the parameters and the results for a trial
  if isdefined(:mutation_stddev_list)  # the  sim_record  initialized here is used in the calls to writeheader() below.
    sim_record = ContVarEvolution.cont_var_result(num_trials,N_list[1],num_subpops,num_attributes_list[1], ngens, burn_in,
       mutation_stddev_list[1], ideal, fit_slope, neutral )
  elseif isdefined(:N_mut_list)
    sim_record = ContVarEvolution.cont_var_result(num_trials,N_list[1],num_subpops,num_attributes_list[1], ngens, burn_in,
       N_mut_list[1]/N_list[1]/100, ideal, fit_slope, neutral )
  end
  # collect these parameter records into the list   sim_record_list_run 
  sim_record_list_run = ContVarEvolution.cont_var_result_type[]
  trial=1
  if isdefined(:mutation_stddev_list)
    for N in N_list
      for num_attributes in num_attributes_list
        for mutation_stddev in mutation_stddev_list
          for trial = 1:num_trials
              sim_record = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
                 mutation_stddev, ideal, fit_slope, neutral )
              Base.push!(sim_record_list_run, sim_record )
          end
        end
      end
    end
  elseif isdefined(:N_mut_list)
    for N in N_list
      for N_mut in N_mut_list
        for num_attributes in num_attributes_list
          for trial = 1:num_trials
            mutation_stddev = N_mut/N
            sim_record = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
                 mutation_stddev, ideal, fit_slope, neutral )
            int_burn_in = Int(round(burn_in*sim_record.N+50.0)) 
            Base.push!(sim_record_list_run, sim_record )
          end
        end
      end
    end
  end
  println("===================================")
  # Run the simulation function "cont_var_simulation" on each parameter record in parallel
  # For each trial, "cont_var_simulation" adds the trial results to the parameter/result record
  # TODO:  change back to pmap
  sim_record_list_result = map(cont_var_simulation, sim_record_list_run )
  trial = 1
  writeheader( STDOUT, sim_record )
  writeheader( stream, sim_record )
  for sim_record_result in sim_record_list_result    # write 1 record per trial
    writerow(stream,trial,sim_record_result)
    writerow(STDOUT,trial,sim_record_result)
    trial += 1
  end
end    

if length(ARGS) == 0
  simname = "examples/example2"
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
println("num_trials: ",num_trials)
run_trials( simname )
