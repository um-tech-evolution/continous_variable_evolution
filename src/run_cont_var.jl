#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run_cv.jl configs/example1
=#
export cont_var_result, print_cont_var_result, run_trial, writeheader, writerow
#include("types.jl")
  
function cont_var_result( num_trials, N::Int64, num_subpops::Int64, num_attributes::Int64, ngens::Int64, burn_in::Number,
     mutation_stddev::Float64, ideal::Float64, fit_slope::Float64, wrap_attributes::Bool, additive_error::Bool, neutral::Bool=false )
  if typeof(burn_in) == Int64
    int_burn_in = burn_in
  else
    int_burn_in = Int(round(burn_in*N+50.0))
  end
  return cont_var_result_type( num_trials, N, num_subpops, num_attributes, ngens, int_burn_in,
      mutation_stddev, ideal, fit_slope, wrap_attributes, additive_error, neutral, 0.0, 0.0, 0.0, 0.0, 0,0,0,0 )
end

function print_cont_var_result( sr::cont_var_result_type )
  println("num_trials: ", sr.num_trials)
  println("N: ", sr.N)
  println("num_subpops: ", sr.num_subpops)
  println("num_attributes: ", sr.num_attributes)
  println("mutation_stddev: ", sr.mutation_stddev)
  println("ideal: ",sr.ideal)
  println("fit_slope: ",sr.fit_slope)
  println("ngens: ", sr.ngens)
  println("wrap attributes: ", sr.wrap_attributes)
  println("additive_error: ", sr.additive_error)
  println("burn_in: ", sr.burn_in)
  println("neutral: ", sr.neutral )
  println("fitness_mean: ", sr.fitness_mean)
  println("fitness_coef_var: ", sr.fitness_coef_var)
  println("attiribute_mean: ", sr.attribute_mean)
  println("attiribute_coef_var: ", sr.attribute_coef_var)
  println("fit diff neg count: ",sr.neg_count)
  println("fit diff neg neutral: ",sr.neg_neutral)
  println("fit diff pos neutral: ",sr.pos_neutral)
  println("fit diff pos count: ",sr.pos_count)
end

function writeheader( stream::IO, sr::cont_var_result_type )
  global mutation_stddev_list
  global N_mut_list
  println("isdefined mutation_stddev_list: ",isdefined(:mutation_stddev_list))
  println("isdefined N_mut_list: ",isdefined(:N_mut_list))
  if isdefined(:N_mut_list)
    N_mut_string = "# using N_mut_list"
  elseif isdefined(:mutation_stddev_list)
    N_mut_string = "# using mutation_stddev_list"
  end
  param_strings = [
    "# $(string(Dates.today()))",
    "# num_trials=$(sr.num_trials)",
    N_mut_string,
    #"# N=$(sr.N)",
    #"# num_attributes=$(sr.num_attributes)",
    "# ngens=$(sr.ngens)",
    "# wrap_attributes =$(sr.wrap_attributes)",
    "# additive_error=$(sr.additive_error)",
    "# int_burn_in=$(sr.int_burn_in)",
    "# neutral=$(sr.neutral)",
    "# ideal=$(sr.ideal)",
    "# fit_slope=$(sr.fit_slope)"]

  write(stream,join(param_strings,"\n"),"\n")
  heads = [
    "N",
    "N_mut",
    "mutation_stddev",
    #"num_emigrants",
    "num_attributes",
    "fitness_mean",
    "fitness_coef_var",
    "attribute_mean",
    "attribute_coef_var",
    "fit_diff_neg_fract",
    "fit_diff_neg_neutral",
    "fit_diff_pos_neutral",
    "fit_diff_pos_fract"]
  write(stream,join(heads,","),"\n")
end
    
function writerow( stream::IO, trial::Int64, sr::cont_var_result_type )
  sum_fitdiff = Float64(sum( (sr.neg_count, sr.neg_neutral, sr.pos_neutral, sr.pos_count) ))
  line = Any[
          sr.N,
          sr.N*sr.mutation_stddev,
          sr.mutation_stddev,
          sr.num_attributes,
          sr.fitness_mean,
          sr.fitness_coef_var,
          sr.attribute_mean,
          sr.attribute_coef_var,
          sr.neg_count/sum_fitdiff,
          sr.neg_neutral/sum_fitdiff,
          sr.pos_neutral/sum_fitdiff,
          sr.pos_count/sum_fitdiff]
  write(stream,join(line,","),"\n")
end


if isdefined(:simtype)
  run_trials()
end
