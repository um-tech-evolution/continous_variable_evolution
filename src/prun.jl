# Example run:  julia -L NeutralEvolution.jl run.jl examples/sn_example1
# Example run:  julia -p 4 -L NeutralEvolution.jl run.jl examples/sn_example1
using NeutralEvolution

function run_trials( simname::AbstractString )
  global num_trials
  global use_population
  global save_populations
  global log_error    # use log normal error:  work in log space
  global wright_fisher_copy  
  #println("log_error: ",log_error)
  if !(simtype == 3 || simtype == 4)
    error("simtype must be 3 or 4 for simple neutral evolution!")
  end
  # set defaults for parameters that may not be included in the config file.
  if !isdefined(:use_population)
    use_population = true
  end
  if !isdefined(:save_populations)
    save_populations = false
  end
  if !isdefined(:log_error)
    println("log_error not defined")
    log_error = false    # default is to not use log error
  end
  if !isdefined(:wright_fisher_copy)
    println("wright_fisher_copy not defined")
    wright_fisher_copy = true     # default is to use wright-fisher copy
  end
  sn = simple_neutral_init( simtype, N, mutstddev, ngens, initial_value, num_trials, record_interval, use_population,
      save_populations, log_error, wright_fisher_copy )
  csn = cummulative_neutral_init( sn )
  print_simple_neutral_params( sn )
  sn_list_run = simple_neutral_type[]
  for t = 1:num_trials
    #=
    if sn.N > 1 || sn.use_population
      sn = simple_neutral_simulation( sn )
    elseif sn.N == 1
      #sn = ces_log( sn )
      sn = ces( sn )
    else
      error(" sn.N must be positive ")
    end
    =#
    Base.push!(sn_list_run,deepcopy(sn))
    #println("trial: ",t,"  avg mean: ",sn.average_attr_mean,"  avg coef var: ",sn.average_attr_coef_var )
    #println("mean history: ",sn.attr_mean_history)
    #println("coef_var history: ",sn.attr_coef_var_history)
  end
  run_sim = (sn.N > 1 || sn.use_population) ? simple_neutral_simulation : ces 
  sn_list_result = pmap( run_sim, sn_list_run )
  #for t = 1:num_trials
  for rsn in sn_list_result
    accumulate_results( rsn, csn )
  end
  #print_cummulative_neutral( csn )
  #print_results( csn )
  #=
  println("saved_populations: ")
  for i = 1:csn.num_trials
    println(csn.saved_populations[:,i])
  end
  =#
  if !csn.save_populations
    writeheader( STDOUT, csn )
    writerows( STDOUT, csn )
    open("$(simname).csv","w") do stream
      writeheader( stream, csn )
      writerows( stream, csn )
    end
  else
    writeheader_populations( STDOUT, csn )
    writerows_populations( STDOUT, csn )
    open("$(simname)_pops.csv","w") do stream
      writeheader_populations( stream, csn )
      writerows_populations( stream, csn )
    end
  end
end

if length(ARGS) == 0
  simname = "examples/sn_example1"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
println("simname: ",simname)
include("$(simname).jl")
println("simtype: ",simtype)
run_trials( simname )

