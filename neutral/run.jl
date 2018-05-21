# top-level functions to run neutral.jl without parallel map
# Example run:  julia run.jl examples/sn_example1
# Example run:  julia run.jl examples/sn_example1 <seed>    # run with specified random number seed
# Example run:  julia -p 4 run.jl examples/sn_example1 <seed>  # run with 4 processes (cores)
@everywhere include("NeutralEvolution.jl")

function run_trials( simname::AbstractString )
  global num_trials
  global use_population
  global save_populations
  global log_error    # use log normal error:  work in log space
  global wright_fisher_copy,  wright_fisher_copy_list
  global conformist_probability, conformist_probability_list
  global fit_slope
  global N, mutstddev, initial_value
  if !(simtype == 3 || simtype == 4)
    println("Warning! running neutral.jl on simtype == 2 file.")
  end
  # set defaults for parameters that may not be included in the config file.
  if !isdefined(:use_population)
    use_population = true
  end
  if !isdefined(:save_populations)
    save_populations = false
  end
  if !isdefined(:log_error)
    #println("log_error not defined")
    log_error = false    # default is to not use log error
  end
  if isdefined(:wright_fisher_copy_list)
    wright_fisher_copy = wright_fisher_copy_list[1]
  end
  if !isdefined(:wright_fisher_copy)
    #println("wright_fisher_copy not defined")
    wright_fisher_copy = true     # default is to use wright-fisher copy
  end
  if !isdefined(:wright_fisher_copy_list)
    wright_fisher_copy_list = [wright_fisher_copy]
  end
  if isdefined(:conformist_probability_list)
    conformist_probability = conformist_probability_list[1]
  end
  if !isdefined(:conformist_probability)
    #println("conformist_probability not defined")
    conformist_probability = 0.0     # default is no conformist copy
  end
  if !isdefined(:conformist_probability_list)
    conformist_probability_list = [conformist_probability]
  end
  if !isdefined(:fit_slope)
    fit_slope = 0.0
  end
  if !isdefined(:initial_value)
    initial_value = ideal
  end
  if !isdefined(:N) && isdefined(:N_list)
    N = N_list[1]
    println("Warning! Using N = N_list[1] ")
  end
  if !isdefined(:mutstddev) && isdefined(:mutation_stddev_list)
    mutstddev = mutation_stddev_list[1]
    println("Warning! Using mutstddev = mutation_stddev_list[1]")
  end
  sn = simple_neutral_init( simtype, N, mutstddev, ngens, initial_value, num_trials, record_interval, use_population,
      save_populations, log_error, wright_fisher_copy, conformist_probability, neutral, fit_slope )
  #print_simple_neutral_params( sn )
  csn = cummulative_neutral_init( sn )
  if !csn.save_populations
    writeheader( STDOUT, csn )
    open("$(simname).csv","w") do stream
      writeheader( stream, csn )
    end
  else
    writeheader_populations( STDOUT, csn )
    open("$(simname)_pops.csv","w") do stream
      writeheader_populations( stream, csn )
    end
  end
  for cf_prob in conformist_probability_list
    for wf_copy in wright_fisher_copy_list
      csn = cummulative_neutral_init( sn )
      sn.wright_fisher_copy = wf_copy
      sn.conformist_probability = cf_prob
      if nworkers() > 1
        sn_list_run = simple_neutral_type[]
        run_sim = (sn.N > 1 || sn.use_population) ? simple_neutral_simulation : ces
        for t = 1:num_trials
          Base.push!(sn_list_run,deepcopy(sn))
        end
        # TODO:  change map to pmap
        sn_list_result = pmap( run_sim, sn_list_run )
        #sn_list_result = map( run_sim, sn_list_run )
        for rsn in sn_list_result
          accumulate_results( rsn, csn )
        end
      else
        for t = 1:num_trials
          if sn.N > 1 || use_population
            sn = simple_neutral_simulation( sn )
          elseif sn.N == 1
            #sn = ces_log( sn )
            sn = ces( sn )
          else
            error(" sn.N must be positive ")
          end
          accumulate_results( sn, csn )
        end
      end
      if !csn.save_populations
        writerows( STDOUT, csn )
        open("$(simname).csv","a") do stream
          writerows( stream, csn )
        end
      else
        writerows_populations( STDOUT, csn )
        open("$(simname)_pops.csv","a") do stream
          writerows_populations( stream, csn )
        end
      end
    end
  end
end

function coef_var( lst )
  return std(lst)/mean(lst)
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

