# Simple Wright-Fisher model for neutral continuous variable evolution
# Command line:  julia simple_neutral.jl
# Command line:  julia simple_neutral.jl <seed>
#    where seed is a random number seed.
# When run with a specific seed, the output should agree running the full simulation
#    for neutral multiplicative error with a single attribute (with the same seed).
using DataFrames
using CSV

type simple_neutral_type
  N::Int64   # popsize
  mutstddev::Float64
  ngens::Int64
  initial_value::Float64
  num_trials::Int64
  record_interval::Int64
  use_population::Bool     # must be true if N>1
  #use_log_update::Bool     # for now, must be false if N>1
  average_attr_mean::Float64
  average_attr_coef_var::Float64
  attr_mean_history::Vector{Float64}
  attr_median_history::Vector{Float64}
  attr_coef_var_history::Vector{Float64}
end

function simple_neutral_init( N::Int64, mutstddev::Float64, ngens::Int64, initial_value::Float64, num_trials::Int64, record_interval::Int64, use_population::Bool=true )
  num_results = Int(ceil(ngens/record_interval))
  simple_neutral_type( N, mutstddev, ngens, initial_value, num_trials, record_interval, use_population, 0.0, 0.0, 
    fill(0.0,num_results), 
    fill(0.0,num_results), 
    fill(0.0,num_results)) 
end


function print_simple_neutral_params( sn::simple_neutral_type )
  global use_population
  print("N: ",sn.N,"  mutstddev: ",sn.mutstddev," ngens: ",sn.ngens,"  initial_value: ",sn.initial_value,"  num_trials: ",sn.num_trials,
        " record_interval: ",sn.record_interval )
  if isdefined(:use_population)
    println("  use_population: ",use_population) 
  else 
    println()
  end
end

type cummulative_neutral_type
  N::Int64  # popsize
  mutstddev::Float64
  ngens::Int64
  initial_value::Float64
  num_trials::Int64
  record_interval::Int64
  use_population::Bool     # must be true if N>1
  count_trials::Int64
  gens_recorded::Vector{Int64}
  attr_mean_sum::Vector{Float64}
  attr_mean_sum_sq::Vector{Float64}
  attr_median_sum::Vector{Float64}
  attr_median_sum_sq::Vector{Float64}
  attr_coef_var_sum::Vector{Float64}
  attr_coef_var_sum_sq::Vector{Float64}
end

function cummulative_neutral_init( sn::simple_neutral_type )
  num_results = Int(ceil(ngens/record_interval))
  cummulative_neutral_type( sn.N, sn.mutstddev, sn.ngens, sn.initial_value, sn.num_trials, sn.record_interval, sn.use_population,
    0, fill(0,num_results), fill(0.0,num_results), fill(0.0,num_results), 
    fill(0.0,num_results), fill(0.0,num_results), fill(0.0,num_results), fill(0.0,num_results)) 
end


function print_cummulative_neutral( csn::cummulative_neutral_type )
  println("N: ",csn.N,"  mutstddev: ",csn.mutstddev," ngens: ",csn.ngens," record_interval: ",csn.record_interval,
      "  count_trials: ",csn.count_trials )
  println("mean sum   : ",csn.attr_mean_sum)
  println("mean sum sq: ",csn.attr_mean_sum_sq)
  println("cvar sum   : ",csn.attr_coef_var_sum)
  println("cvar sum sq: ",csn.attr_coef_var_sum_sq)
end

function copy_pop( sn::simple_neutral_type, pop::Vector{Float64} )
  i = 1
  for p in pop
    rn = randn()
    #println("rn: ",rn)
    mult = 1.0 + sn.mutstddev*rn
    while mult <= 1.e-6
      rn = randn()
      #println("rn: ",rn)
      mult = 1.0 + sn.mutstddev*rn
    end
    pop[i] = p*mult
    if pop[i] <= 0.0
      println("negative pop member at g= ",g)
    end
    i+= 1
  end
  return pop
end


function simple_neutral_simulation( sn::simple_neutral_type )
  int_burn_in = 0
  cumm_attr_mean = 0.0
  cumm_attr_coef_var = 0.0
  sum_gens = 0
  record_index = 1
  pop = fill( sn.initial_value, sn.N )
  for g = 1:(int_burn_in + sn.ngens )
    i = 1
    #println("before copy g: ",g,"  pop: ",pop)
    new_pop = copy_pop( sn, pop )
    #println("after copy g: ",g,"  new_pop: ",new_pop)
    pop = [ new_pop[ rand(1:sn.N) ] for j = 1:sn.N ]
    #println("after WF g: ",g,"  pop: ",pop)
    cumm_attr_mean += mean(pop)
    cumm_attr_coef_var += coef_var(pop)
    sum_gens += 1
    if (g-1) % record_interval == 0
      sn.attr_mean_history[record_index] = mean(pop)
      sn.attr_median_history[record_index] = median(pop)
      sn.attr_coef_var_history[record_index] = coef_var(pop)
      #println("g: ",g,"  attr_mean: ",sn.attr_mean_history[record_index],"  attr_coef_var: ",sn.attr_coef_var_history[record_index])
      record_index += 1
    end
    #pop = deepcopy(new_pop)
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= sn.ngens
  cumm_attr_coef_var /= sn.ngens
  #println("cumm_attr_mean: ",cumm_attr_mean)
  #println("cumm_attr_coef_var: ",cumm_attr_coef_var)
  sn.average_attr_mean = cumm_attr_mean
  sn.average_attr_coef_var = cumm_attr_coef_var
  return sn
end

# single individual copy error simulation according to Hamilton and Buchanan
# S is log of size
function ces_log( sn::simple_neutral_type )
  cumm_attr_mean = 0.0
  sum_gens = 0
  record_index = 1
  S = sn.initial_value
  for g = 1:sn.ngens
    S += log(1.0+sn.mutstddev*randn())
    sum_gens += 1
    #cumm_attr_mean += S
    cumm_attr_mean += exp(S)
    if (g-1) % record_interval == 0
      #sn.attr_mean_history[record_index] = S
      sn.attr_mean_history[record_index] = exp(S)
      sn.attr_coef_var_history[record_index] = 0.0
      record_index += 1
    end
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= sn.ngens
  sn.average_attr_mean = cumm_attr_mean
  sn.average_attr_coef_var = 0.0
  return sn
end
    
# Single individual copy error simulation according to Eerkins and Lipo
# Note that s is the size (not the log of size)
# See ces_log 
function ces( sn::simple_neutral_type )
  cumm_attr_mean = 0.0
  sum_gens = 0
  record_index = 1
  s = sn.initial_value
  for g = 1:sn.ngens
    mult = (1.0+sn.mutstddev*randn())
    while mult < 1.e-6
      mult = (1.0+sn.mutstddev*randn())
    end
    s *= mult
    sum_gens += 1
    #cumm_attr_mean += s
    cumm_attr_mean += s
    if (g-1) % record_interval == 0
      #sn.attr_mean_history[record_index] = s
      sn.attr_mean_history[record_index] = s
      sn.attr_coef_var_history[record_index] = 0.0
      record_index += 1
    end
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= sn.ngens
  sn.average_attr_mean = cumm_attr_mean
  sn.average_attr_coef_var = 0.0
  return sn
end
    
function accumulate_results( sn::simple_neutral_type, csn::cummulative_neutral_type )
  num_results = Int(ceil(sn.ngens/sn.record_interval))
  record_index = 1
  for g = 1:sn.ngens
    if (g-1) % record_interval == 0
      csn.gens_recorded[record_index] = g
      record_index += 1
    end
  end
  for i = 1:num_results
    csn.attr_mean_sum[i] += sn.attr_mean_history[i]
    csn.attr_mean_sum_sq[i] += sn.attr_mean_history[i]^2
    csn.attr_median_sum[i] += sn.attr_median_history[i]
    csn.attr_median_sum_sq[i] += sn.attr_median_history[i]^2
    csn.attr_coef_var_sum[i] += sn.attr_coef_var_history[i]
    csn.attr_coef_var_sum_sq[i] += sn.attr_coef_var_history[i]^2
  end
  csn.count_trials += 1
end

function print_results( csn::cummulative_neutral_type )
  num_results = Int(ceil(csn.ngens/csn.record_interval))
  mean_squared_list = map(x->x^2, csn.attr_mean_sum)
  mean_sq_list = csn.attr_mean_sum .* csn.attr_mean_sum
  @assert reduce(&, mean_squared_list .== mean_sq_list )
  if csn.N > 1 
    coef_var_squared_list = map(x->x^2, csn.attr_coef_var_sum)
    println("coef_var_squared_list: ",coef_var_squared_list)
    coef_var_sq_list = csn.attr_coef_var_sum .* csn.attr_coef_var_sum
    println("coef_var_squ_list: ",coef_var_sq_list)
    @assert reduce(&, coef_var_squared_list .== coef_var_sq_list )
    cv_stdevs = csn.attr_coef_var_sum_sq/(csn.count_trials-1) - coef_var_squared_list/csn.count_trials/(csn.count_trials-1)
  else
    cv_stdevs = fill(0.0,num_results)
  end
  println("Gens: ",csn.gens_recorded)
  println("means    : ",csn.attr_mean_sum/csn.count_trials)
  println("mn_stdvs : ",csn.attr_mean_sum_sq/(csn.count_trials-1) - mean_squared_list/csn.count_trials/(csn.count_trials-1))
  println("coef_vars: ",csn.attr_coef_var_sum/csn.count_trials)
  println("cv_stdevs: ",cv_stdevs)
  df = DataFrame()
  df[:Gens] = csn.gens_recorded
  df[:means] = csn.attr_mean_sum/csn.count_trials
  df[:mean_stdvs] = csn.attr_mean_sum_sq/(csn.count_trials-1) - mean_squared_list/csn.count_trials/(csn.count_trials-1)
  df[:medians] = csn.attr_median_sum/csn.count_trials
  df[:median_stdvs] = csn.attr_median_sum_sq/(csn.count_trials-1) - median_squared_list/csn.count_trials/(csn.count_trials-1)
  df[:coef_vars] = csn.attr_coef_var_sum/csn.count_trials
  df[:cv_stdevs] = cv_stdevs
  CSV.write("$(simname).csv",df)
end

function writeheader( stream::IO, csn::cummulative_neutral_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((simtype==3)?"neutral continuous var model":"invalid simtype")",
    "# N=$(csn.N)",
    "# mutation stddev=$(csn.mutstddev)",
    "# ngens=$(csn.ngens)",
    "# num_trials=$(csn.num_trials)",
    "# initial_value=$(csn.initial_value)",
    "# record_interval=$(csn.record_interval)",
    "# use_population=$(csn.use_population)"
    ]
  head_strings = [
    "Gen",
    "average mean",
    "stddev of mean",
    "average median",
    "stddev of median",
    "average coef var",
    "stddev of coef var"
    ]
  write(stream, join(param_strings, "\n"), "\n")
  write(stream, join(head_strings, ","), "\n")
end

function writerows( stream::IO, csn::cummulative_neutral_type )
  num_results = Int(ceil(csn.ngens/csn.record_interval))
  means = csn.attr_mean_sum/csn.count_trials
  mean_squared_list = map(x->x^2, csn.attr_mean_sum)
  mean_sq_list = csn.attr_mean_sum .* csn.attr_mean_sum
  @assert reduce(&, mean_squared_list .== mean_sq_list )
  mean_stdvs = csn.attr_mean_sum_sq/(csn.count_trials-1) - mean_squared_list/csn.count_trials/(csn.count_trials-1)
  medians = csn.attr_median_sum/csn.count_trials
  median_squared_list = map(x->x^2, csn.attr_median_sum)
  median_sq_list = csn.attr_median_sum .* csn.attr_median_sum
  @assert reduce(&, median_squared_list .== median_sq_list )
  median_stdvs = csn.attr_median_sum_sq/(csn.count_trials-1) - median_squared_list/csn.count_trials/(csn.count_trials-1)
  coef_vars = csn.attr_coef_var_sum/csn.count_trials
  if csn.N > 1
    coef_var_squared_list = map(x->x^2, csn.attr_coef_var_sum)
    #println("coef_var_squared_list: ",coef_var_squared_list)
    coef_var_sq_list = csn.attr_coef_var_sum .* csn.attr_coef_var_sum
    #println("coef_var_squ_list: ",coef_var_sq_list)
    @assert reduce(&, coef_var_squared_list .== coef_var_sq_list )
    cv_stdevs = csn.attr_coef_var_sum_sq/(csn.count_trials-1) - coef_var_squared_list/csn.count_trials/(csn.count_trials-1)
  else
    cv_stdevs = fill(0.0,num_results)
  end
  for i = 1:num_results
    values = Number[
      csn.gens_recorded[i],
      means[i],
      mean_stdvs[i],
      medians[i],
      median_stdvs[i],
      coef_vars[i],
      cv_stdevs[i]
      ] 
    write(stream, join(values,","), "\n")
  end
end
  

function run_trials( simname::AbstractString )
  global num_trials
  global use_population
  if !(simtype == 3 || simtype == 4)
    error("simtype must be 3 or 4 for simple neutral evolution!")
  end
  overall_avg_mean = 0.0
  overall_avg_coef_var = 0.0
  if isdefined(:use_population)
    sn = simple_neutral_init( N, mutstddev, ngens, initial_value, num_trials, record_interval, use_population )
  else
    sn = simple_neutral_init( N, mutstddev, ngens, initial_value, num_trials, record_interval )
  end
  csn = cummulative_neutral_init( sn )
  print_simple_neutral_params( sn )
  for t = 1:num_trials
    if sn.N > 1 || use_population
      sn = simple_neutral_simulation( sn )
    elseif sn.N == 1
      #sn = ces_log( sn )
      sn = ces( sn )
    else
      error(" sn.N must be positive ")
    end
    #println("trial: ",t,"  avg mean: ",sn.average_attr_mean,"  avg coef var: ",sn.average_attr_coef_var )
    #println("mean history: ",sn.attr_mean_history)
    #println("coef_var history: ",sn.attr_coef_var_history)
    accumulate_results( sn, csn )
  end
  overall_avg_mean /= num_trials
  overall_avg_coef_var /= num_trials
  #print_cummulative_neutral( csn )
  #print_results( csn )
  writeheader( STDOUT, csn )
  writerows( STDOUT, csn )
  open("$(simname).csv","w") do stream
    writeheader( stream, csn )
    writerows( stream, csn )
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

