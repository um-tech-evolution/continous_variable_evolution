# See run.jl for command-line example runs.
export run_trials, cont_var_result, print_cont_var_result, run_trial, writeheader, writerow, my_isdefined, check_parameters

@doc """ function run_trials(sim_record::ContVarEvolution.cont_var_result_type )
  Runs multiple trials of the simulation using the parameter file "\$(simname).jl" where  simname is the first command line argument.
  Trials can be run in parallel using the Julia parallel map (pmap) facility.
"""
function run_trials( simname::AbstractString, sim_record::ContVarEvolution.cont_var_result_type ) 
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  # collect these parameter records into the list   sim_record_list_run 
  sim_record_list_run = ContVarEvolution.cont_var_result_type[]
  trial=1
  #if my_isdefined(:mutation_stddev_list)
  if length(sim_record.mutation_stddev_list) > 0
    for N in sim_record.N_list
      if !(typeof(sim_record.burn_in) <: Int)
        sim_record.int_burn_in = Int(round(sim_record.burn_in*N+50.0))
        # println("sim_record.int_burn_in: ",sim_record.int_burn_in)
      end
      for num_attributes in sim_record.num_attributes_list
        for mutation_stddev in sim_record.mutation_stddev_list
          for trial = 1:sim_record.num_trials
            sim_record_run = deepcopy( sim_record )
            sim_record_run.N = N   
            sim_record_run.num_attributes = num_attributes   
            sim_record_run.mutation_stddev = mutation_stddev   
            Base.push!(sim_record_list_run, sim_record_run )
          end
        end
      end
    end
  #elseif my_isdefined(:N_mut_list)
  elseif length(sim_record.N_mut_list) > 0 
    for N in sim_record.N_list
      if !(typeof(sim_record.burn_in) <: Int)
        sim_record.int_burn_in = Int(round(sim_record.burn_in*N+50.0))
        # println("sim_record.int_burn_in: ",sim_record.int_burn_in)
      end
      for num_attributes in sim_record.num_attributes_list
        for N_mut in sim_record.N_mut_list
          for trial = 1:sim_record.num_trials
            sim_record_run = deepcopy( sim_record )
            sim_record_run.N = N   
            sim_record_run.num_attributes = num_attributes   
            sim_record_run.mutation_stddev = N_mut/N
            #println("N:",sim_record_run.N,"  num_attr:",sim_record_run.num_attributes,"  mutation_stddev: ",sim_record_run.mutation_stddev)
            Base.push!(sim_record_list_run, sim_record_run )
          end
        end
      end
    end
  else
    error("Either mutation_stddev_list or N_mut_list must be defined in the configuration file.")
  end
  println("===================================")
  # Run the simulation function "cont_var_simulation" on each parameter record in parallel
  # For each trial, "cont_var_simulation" adds the trial results to the parameter/result record
  # You may want to change "pmap" to "map" for debugging purposes if you are getting complicated error messages.
  sim_record_list_result = pmap(cont_var_simulation, sim_record_list_run )
  #sim_record_list_result = map(cont_var_simulation, sim_record_list_run )
  # The next line is for multigeneration output
  #println("    N, mutstd,   g, mean,median,coefvar,entropy")
  trial = 1
  writeheader( Base.stdout, sim_record )   # change STDOUT to stdout for julia v7 (but then will fail in julia v6)
  writeheader( stream, sim_record )  
  for sim_record_result in sim_record_list_result    # write 1 record per trial
    writerow(Base.stdout,trial,sim_record_result)  # change STDOUT to stdout for julia v7
    writerow(stream,trial,sim_record_result)
    trial += 1
  end
end    

@doc """ function cont_var_result()
   Constructs a cont_var_result record with the fields corresponding to parameters given in the argument list set, and the result fields set to 0 or 0.0.
"""
function cont_var_result( N_list::Vector{Int64}, num_attributes_list::Vector{Int64}, mutation_stddev_list::Vector{Float64}, N_mut_list::Vector{Float64},
        num_trials, N::Int64, num_subpops::Int64, num_attributes::Int64, ngens::Int64, burn_in::Number, int_burn_in::Int64,
        mutation_stddev::Float64, ideal::Float64, fit_slope::Float64, neutral::Bool, w::Float64  )
  if typeof(burn_in) == Int64
    int_burn_in = burn_in
  else
    #int_burn_in = Int(round(burn_in*N+50.0))
    int_burn_in = -1   # to be set later in run_trials
  end
  return cont_var_result_type( N_list, num_attributes_list, mutation_stddev_list, N_mut_list, num_trials, N, num_subpops, num_attributes, ngens, burn_in, int_burn_in, 
      mutation_stddev, ideal, fit_slope, neutral, w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0,0,0 )
end

@doc """ function writeheader()
  Writes the parameters and the header line for the output CSV file.
  The parameters that do not vary by trial are first written as comment lines where the comment character is "#".
  Then the header line is written.  The headers of course must correspond to the values written by the writerows function.
"""
function writeheader( stream::IO, sim_record::cont_var_result_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# N_list=$(sim_record.N_list)",
    "# num_attributes_list=$(sim_record.num_attributes_list)",
    "# mutation_stddev_list=$(sim_record.mutation_stddev_list)",
    "# N_mut_list=$(sim_record.N_mut_list)",
    "# num_trials=$(sim_record.num_trials)",
    #"# num_attributes=$(sim_record.num_attributes)",
    "# ngens=$(sim_record.ngens)",
    "# int_burn_in=$(sim_record.int_burn_in)",
    "# neutral=$(sim_record.neutral)",
    "# ideal=$(sim_record.ideal)",
    "# fit_slope=$(sim_record.fit_slope)"]
    #"# w=$(sim_record.w)"]  # Added as a column in the output table  8/9/19

  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "N",
    "N_mut",
    "mutation_stddev",
    "num_attributes",
    "int_burn_in",
    "w",
    "attribute_mean",
    "attribute_median",
    "attribute_coef_var",
    "attribute_entropy"
  ]
  fitness_heads = [   # not written if the simulation is neutral
    "fitness_mean",
    "fitness_median",
    "fitness_coef_var",
    "fit_diff_neg_fract",
    "fit_diff_neg_neutral",
    "fit_diff_pos_neutral",
    "fit_diff_pos_fract"
  ]
  if !sim_record.neutral
    heads = vcat(heads,fitness_heads)
  end
  write(stream,join(heads,","),"\n")
end

@doc """ function writerow()
  Writes one row of the CSV file, where a row corresponds to 1 trial of the simulation.
"""    
function writerow( stream::IO, trial::Int64, sim_record::cont_var_result_type )
  sum_fitdiff = Float64(sum( (sim_record.neg_count, sim_record.neg_neutral, sim_record.pos_neutral, sim_record.pos_count) ))
  values = Any[
          sim_record.N,
          sim_record.N*sim_record.mutation_stddev,
          sim_record.mutation_stddev,
          sim_record.num_attributes,
          sim_record.int_burn_in,
          sim_record.w,
          sim_record.attribute_mean,
          sim_record.attribute_median,
          sim_record.attribute_coef_var,
          sim_record.attribute_entropy
  ]
  fitness_values = [
          sim_record.fitness_mean,
          sim_record.fitness_median,
          sim_record.fitness_coef_var,
          sim_record.neg_count/sum_fitdiff,
          sim_record.neg_neutral/sum_fitdiff,
          sim_record.pos_neutral/sum_fitdiff,
          sim_record.pos_count/sum_fitdiff
  ]
  if !sim_record.neutral
    values = vcat(values,fitness_values)
  end
  write(stream,join(values,","),"\n")
end
