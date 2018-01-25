# See run.jl for command-line example runs.
export cont_var_result, print_cont_var_result, run_trial, writeheader, writerow
#include("types.jl")
  
function cont_var_result( num_trials, N::Int64, num_subpops::Int64, num_attributes::Int64, ngens::Int64, burn_in::Number,
     mutation_stddev::Float64, ideal::Float64, fit_slope::Float64, neutral::Bool=false )
  if typeof(burn_in) == Int64
    int_burn_in = burn_in
  else
    int_burn_in = Int(round(burn_in*N+50.0))
  end
  return cont_var_result_type( num_trials, N, num_subpops, num_attributes, ngens, int_burn_in,
      mutation_stddev, ideal, fit_slope, neutral, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0,0,0 )
end

# TODO This function is never called and so can be deleted
function print_cont_var_result( sim_record::cont_var_result_type )
  println("num_trials: ", sim_record.num_trials)
  println("N: ", sim_record.N)
  println("num_subpops: ", sim_record.num_subpops)
  println("num_attributes: ", sim_record.num_attributes)
  println("mutation_stddev: ", sim_record.mutation_stddev)
  println("ideal: ",sim_record.ideal)
  println("fit_slope: ",sim_record.fit_slope)
  println("ngens: ", sim_record.ngens)
  println("burn_in: ", sim_record.burn_in)
  println("neutral: ", sim_record.neutral )
  println("fitness_mean: ", sim_record.fitness_mean)
  println("fitness_median: ", sim_record.fitness_mean)
  println("fitness_coef_var: ", sim_record.fitness_coef_var)
  println("attiribute_mean: ", sim_record.attribute_mean)
  println("attiribute_median: ", sim_record.attribute_mean)
  println("attiribute_coef_var: ", sim_record.attribute_coef_var)
  println("fit diff neg count: ",sim_record.neg_count)
  println("fit diff neg neutral: ",sim_record.neg_neutral)
  println("fit diff pos neutral: ",sim_record.pos_neutral)
  println("fit diff pos count: ",sim_record.pos_count)
end

@doc """ function writeheader()
  Writes the parameters and the header line for the output CSV file.
  The parameters that do not vary by trial are first written as comment lines where the comment character is "#".
  Then the header line is written.  The headers of course must correspond to the values written by the writerows function.
"""
function writeheader( stream::IO, sim_record::cont_var_result_type )
  global mutation_stddev_list
  global N_mut_list
  #println("isdefined mutation_stddev_list: ",isdefined(:mutation_stddev_list))
  #println("isdefined N_mut_list: ",isdefined(:N_mut_list))
  if isdefined(:N_mut_list)
    N_mut_string = "# using N_mut_list"
  elseif isdefined(:mutation_stddev_list)
    N_mut_string = "# using mutation_stddev_list"
  end
  param_strings = [
    "# $(string(Dates.today()))",
    "# num_trials=$(sim_record.num_trials)",
    N_mut_string,
    #"# N=$(sim_record.N)",
    #"# num_attributes=$(sim_record.num_attributes)",
    "# ngens=$(sim_record.ngens)",
    "# neutral=$(sim_record.neutral)",
    "# ideal=$(sim_record.ideal)",
    "# fit_slope=$(sim_record.fit_slope)"]

  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "N",
    "N_mut",
    "mutation_stddev",
    "num_attributes",
    "int_burn_in",
    "attribute_mean",
    "attribute_median",
    "attribute_coef_var"
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
          sim_record.attribute_mean,
          sim_record.attribute_median,
          sim_record.attribute_coef_var
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

#= for testing purposes  TODO:  delete
if isdefined(:simtype)
  run_trials()
end
=#
