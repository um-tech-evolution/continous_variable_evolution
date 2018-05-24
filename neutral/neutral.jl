# Simple Wright-Fisher model for neutral continuous variable evolution
# Example runs with parallel map over trais.
# Example run:  julia -L NeutralEvolution.jl nrun.jl examples/sn_example1
# Example run:  julia -L NeutralEvolution.jl nrun.jl examples/sn_example1 <seed>
# Example run:  julia -p 4 -L NeutralEvolution.jl nrun.jl examples/sn_example1
#    where seed is a random number seed.
# When run with a specific seed, the output should agree running the full simulation
#    for neutral multiplicative error with a single attribute (with the same seed).
export simple_neutral_type, simple_neutral_init, cummulative_neutral_type, cummulative_neutral_init,
      print_simple_neutral_params, simple_neutral_simulation, ces, accumulate_results, writeheader,
      writerows, writeheader_populations, writerows_populations, coef_var, inverse_fitness

try
  using Dates
  using Printf
catch
end

mutable struct simple_neutral_type
  simtype::Int64
  N::Int64   # popsize
  mutstddev::Float64
  ngens::Int64
  initial_value::Float64
  num_trials::Int64
  record_interval::Int64
  use_population::Bool     # must be true if N>1
  save_populations::Bool    # Save the final population for each trial as a column of the CSV file
  log_error::Bool           # use log normally distributed error to work in log space
  wright_fisher_copy::Bool  # use Wright-Fisher copy to enable population level drift (default true)
  conformist_probability::Float64   # Probability of copying the mean population value rather than Wright-Fisher copy
  neutral::Bool             # if true, no selection 
  fit_slope::Float64        # Slope of the inverse fitness function.  0.0 means neutral
  fit_funct::Function       # Fitness function
  average_attr_mean::Float64
  average_attr_median::Float64
  average_attr_coef_var::Float64
  attr_mean_history::Vector{Float64}
  attr_median_history::Vector{Float64}
  attr_coef_var_history::Vector{Float64}
  saved_population::Vector{Float64}
end

function inverse_fitness( a::Float64, initial_value::Float64, fit_slope::Float64=1.0 )
  return 1.0/(abs(fit_slope*a-initial_value)+1.0)
end

function simple_neutral_init( simtype::Int64,N::Int64, mutstddev::Float64, ngens::Int64, initial_value::Float64, num_trials::Int64, 
    record_interval::Int64, use_population::Bool=true, save_population::Bool=false, 
    log_error::Bool=false, wright_fisher_copy::Bool=true, conformist_probability::Float64=0.0, neutral::Bool=true,
    fit_slope::Float64=1.0, fit_funct::Function=inverse_fitness)
  num_results = Int(ceil(ngens/record_interval))
  saved_pop_empty = try  Array{Float64,1}(undef,0) catch; Array{Float64,1}(0) end
  sn = simple_neutral_type( simtype, N, mutstddev, ngens, initial_value, num_trials, record_interval, use_population, save_population,
    log_error, wright_fisher_copy, conformist_probability, neutral, fit_slope, inverse_fitness,
    0.0, 0.0, 0.0,
    fill(0.0,num_results),
    fill(0.0,num_results), 
    fill(0.0,num_results),
    saved_pop_empty)
  #println("log_error: ",sn.log_error)
  sn
end

function print_simple_neutral_params( sn::simple_neutral_type )
  global use_population
  println("N: ",sn.N,"  mutstddev: ",sn.mutstddev," ngens: ",sn.ngens,"  initial_value: ",sn.initial_value,"  num_trials: ",sn.num_trials,
        " record_interval: ",sn.record_interval )
  if isdefined(:use_population)
    print("  use_population: ",sn.use_population)
  end
  if isdefined(:save_populations)
    print("  save_populations: ",sn.save_populations)
  end
  if isdefined(:conformist_probability)
    println("  conformist_probability: ",sn.conformist_probability) 
  end
  if isdefined(:fit_slope)
    println("  fit_slope: ",sn.fit_slope) 
  end
  println()
  #println("log_error: ",sn.log_error,"  wright_fisher_copy: ",sn.wright_fisher_copy)
end

mutable struct cummulative_neutral_type
  simtype::Int64
  N::Int64  # popsize
  mutstddev::Float64
  ngens::Int64
  initial_value::Float64
  num_trials::Int64
  record_interval::Int64
  use_population::Bool     # must be true if N>1
  save_populations::Bool    # Save the final population for each trial as a column of the CSV file
  log_error::Bool
  wright_fisher_copy::Bool
  conformist_probability::Float64
  neutral::Bool
  fit_slope::Float64
  count_trials::Int64
  gens_recorded::Vector{Int64}
  attr_mean_sum::Vector{Float64}
  attr_mean_sum_sq::Vector{Float64}
  attr_median_sum::Vector{Float64}
  attr_median_sum_sq::Vector{Float64}
  attr_coef_var_sum::Vector{Float64}
  attr_coef_var_sum_sq::Vector{Float64}
  saved_populations::Array{Float64,2}
end

function cummulative_neutral_init( sn::simple_neutral_type )
  num_results = Int(ceil(sn.ngens/sn.record_interval))
  saved_pops_empty = try  Array{Float64,2}(undef,0,0) catch; Array{Float64,2}(0,0) end
  csn = cummulative_neutral_type( sn.simtype, sn.N, sn.mutstddev, sn.ngens, sn.initial_value, sn.num_trials, sn.record_interval,
      sn.use_population, sn.save_populations, sn.log_error, sn.wright_fisher_copy, sn.conformist_probability, sn.neutral, sn.fit_slope,
      0, fill(0,num_results), fill(0.0,num_results), fill(0.0,num_results), fill(0.0,num_results),
      fill(0.0,num_results), fill(0.0,num_results), fill(0.0,num_results), saved_pops_empty)
  if sn.save_populations
    csn.saved_populations = fill(0.0, (sn.N, sn.num_trials) )
  end
  csn
end

function print_cummulative_neutral( csn::cummulative_neutral_type )
  println("N: ",csn.N,"  mutstddev: ",csn.mutstddev," ngens: ",csn.ngens," record_interval: ",csn.record_interval,
      "  count_trials: ",csn.count_trials )
  println("mean sum   : ",csn.attr_mean_sum)
  println("mean sum sq: ",csn.attr_mean_sum_sq)
  println("cvar sum   : ",csn.attr_coef_var_sum)
  println("cvar sum sq: ",csn.attr_coef_var_sum_sq)
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
    new_pop = mutate_pop( sn, pop )
    #println("after mutate g: ",g,"  new_pop: ",new_pop)
    if sn.neutral || sn.fit_slope == 0.0   # neutral, no selection
      if sn.wright_fisher_copy && (sn.neutral || sn.fit_slope == 0.0)
        #pop = [ new_pop[ sn.N>1 ? rand(1:sn.N) : 1 ] for j = 1:sn.N ]
        if sn.N > 1
          r = rand(1:sn.N,sn.N)
        else
          r = fill(1,sn.N)
        end
        #println("r: ",r)
        pop = new_pop[r]
      else
        #pop = deepcopy(new_pop)   # not sure if deepcopy is necessary
        pop = new_pop
      end
    elseif !sn.neutral && sn.fit_slope > 0.0
      pop = propsel_funct( new_pop, sn.fit_funct, sn.initial_value, sn.fit_slope )
    end
    #println("after WF g: ",g,"  pop: ",pop)
    if sn.conformist_probability > 0.0  # unclear if conformist copy should happen after wf copy or propsel
      conformist_copy!( pop, sn )
    end
    #println("after CC g: ",g,"  pop: ",pop)
    cumm_attr_mean += mean(pop)
    cumm_attr_coef_var += coef_var(pop)
    sum_gens += 1
    if (g-1) % sn.record_interval == 0
      sn.attr_mean_history[record_index] = mean(pop)
      sn.attr_median_history[record_index] = median(pop)
      sn.attr_coef_var_history[record_index] = coef_var(pop)
      #println("g: ",g,"  attr_mean: ",sn.attr_mean_history[record_index],"  attr_coef_var: ",sn.attr_coef_var_history[record_index])
      record_index += 1
    end
    if g == (int_burn_in + sn.ngens)
      #println(" pop: ",pop)
      sn.saved_population = pop
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

function wright_fisher_copy_function( pop::Vector{Float64}, sn::simple_neutral_type )
  if sn.N > 1
    r = rand(1:sn.N,sn.N)
  else
    r = fill(1,sn.N)
  end
  pop[r]
#=
  if sn.conformist_probability > 0.0
    pop_mean = mean(pop)
    for i = 1:sn.N
      if rand() < sn.conformist_probability
        new_pop[i] = pop_mean
      end
    end
  end
  return new_pop 
=#
end  

function conformist_copy!( pop::Vector{Float64}, sn::simple_neutral_type )
  for i = 1:sn.N
    if rand() < sn.conformist_probability
      pop[i] = mean(pop)
    end
  end
end

function mutate_pop( sn::simple_neutral_type, pop::Vector{Float64} )
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
    #println("mult: ",mult)
    if sn.log_error
      pop[i] = p+log(mult)
    else
      pop[i] = p*mult
      @assert pop[i] > 0.0
    end
    i+= 1
  end
  return pop
end

# Single individual copy error simulation according to Eerkins and Lipo
# Note that s is the size (not the log of size)
function ces( sn::simple_neutral_type )
  cumm_attr_mean = 0.0
  cumm_attr_median = 0.0
  sum_gens = 0
  record_index = 1
  s = sn.initial_value
  for g = 1:sn.ngens
    rn = randn()
    #println("rn: ",rn)
    mult = 1.0+sn.mutstddev*rn
    #println("mult: ",mult)
    while mult < 1.e-6
      rn = randn()
      #println("rn: ",rn)
      mult = 1.0+sn.mutstddev*rn
    end
    if sn.log_error
      s += log(mult)
    else
      s *= mult
      @assert s > 0.0
    end
    #println("s: ",s)
    sum_gens += 1
    cumm_attr_mean += s
    cumm_attr_median += s
    if (g-1) % record_interval == 0
      sn.attr_mean_history[record_index] = s
      sn.attr_median_history[record_index] = s
      sn.attr_coef_var_history[record_index] = 0.0
      record_index += 1
    end
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= sn.ngens
  sn.average_attr_mean = cumm_attr_mean
  cumm_attr_median /= sn.ngens
  sn.average_attr_median = cumm_attr_median
  sn.average_attr_coef_var = 0.0
  return sn
end

# single individual copy error simulation according to Hamilton and Buchanan
# S is log of size
function ces_log( sn::simple_neutral_type )
  cumm_attr_mean = 0.0
  cumm_attr_median = 0.0
  sum_gens = 0
  record_index = 1
  S = sn.initial_value
  for g = 1:sn.ngens
    S += log(1.0+sn.mutstddev*randn())
    sum_gens += 1
    #cumm_attr_mean += S
    cumm_attr_mean += exp(S)
    cumm_attr_median += exp(S)
    if (g-1) % record_interval == 0
      #sn.attr_mean_history[record_index] = S
      sn.attr_mean_history[record_index] = exp(S)
      sn.attr_median_history[record_index] = exp(S)
      sn.attr_coef_var_history[record_index] = 0.0
      record_index += 1
    end
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= sn.ngens
  sn.average_attr_mean = cumm_attr_mean
  cumm_attr_median /= sn.ngens
  sn.average_attr_median = cumm_attr_median
  sn.average_attr_coef_var = 0.0
  return sn
end
    
function accumulate_results( sn::simple_neutral_type, csn::cummulative_neutral_type )
  num_results = Int(ceil(sn.ngens/sn.record_interval))
  csn.N = sn.N
  csn.wright_fisher_copy=sn.wright_fisher_copy
  csn.conformist_probability=sn.conformist_probability
  record_index = 1
  for g = 1:sn.ngens
    if (g-1) % sn.record_interval == 0
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
  if csn.save_populations
    csn.saved_populations[:,csn.count_trials] = sn.saved_population
  end
end

function writeheader( stream::IO, csn::cummulative_neutral_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((csn.simtype==3) ? "neutral continuous var model simtype==3" : "invalid simtype")",
    #"# N=$(csn.N)",
    "# mutation stddev=$(csn.mutstddev)",
    "# ngens=$(csn.ngens)",
    "# num_trials=$(csn.num_trials)",
    "# initial_value=$(csn.initial_value)",
    "# record_interval=$(csn.record_interval)",
    "# use_population=$(csn.use_population)",
    "# save_populations=$(csn.save_populations)",
    "# log_error=$(csn.log_error)",
    "# fit_slope=$(csn.fit_slope)",
    "# neutral=$(csn.neutral)"
    #"# wright_fisher_copy=$(csn.wright_fisher_copy)"
    #"# conformist_probability=$(csn.conformist_probability)"
    ]
  head_strings = [
    "N",
    "wf_copy",
    "conf_prob",
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

function writeheader_populations( stream::IO, csn::cummulative_neutral_type )
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((csn.simtype==3) ? "neutral continuous var model simtype==3" : "invalid simtype")",
    "# N=$(csn.N)",
    "# mutation stddev=$(csn.mutstddev)",
    "# ngens=$(csn.ngens)",
    "# num_trials=$(csn.num_trials)",
    "# initial_value=$(csn.initial_value)",
    "# record_interval=$(csn.record_interval)",
    "# use_population=$(csn.use_population)",
    "# save_populations=$(csn.save_populations)",
    "# log_error=$(csn.log_error)",
    #"# wright_fisher_copy=$(csn.wright_fisher_copy)",
    "# conformist_probability=$(csn.conformist_probability)"
    ]
  head_strings = vcat(["wf_copy"],[ "P$i" for i = 1:csn.num_trials ])
  write(stream, join(param_strings, "\n"), "\n")
  write(stream, join(head_strings, ","), "\n")
end

function writerows_populations( stream::IO, csn::cummulative_neutral_type )
  for j = 1:csn.N
    @printf(stream,"%6s,",csn.wright_fisher_copy ? "true" : "false")
    for i = 1:(csn.num_trials-1)
      @printf(stream,"%.4f,",csn.saved_populations[j,i])
    end
    @printf(stream,"%.4f\n",csn.saved_populations[j,csn.num_trials])
  end
end

function writerows( stream::IO, csn::cummulative_neutral_type )
  num_results = Int(ceil(csn.ngens/csn.record_interval))
  means = csn.attr_mean_sum/csn.count_trials
  mean_squared_list = map(x->x^2, csn.attr_mean_sum)
  mean_sq_list = csn.attr_mean_sum .* csn.attr_mean_sum
  @assert reduce(&, mean_squared_list .== mean_sq_list )
  mean_stdvs = sqrt.(csn.attr_mean_sum_sq/(csn.count_trials-1) - mean_squared_list/csn.count_trials/(csn.count_trials-1))
  medians = csn.attr_median_sum/csn.count_trials
  median_squared_list = map(x->x^2, csn.attr_median_sum)
  median_sq_list = csn.attr_median_sum .* csn.attr_median_sum
  @assert reduce(&, median_squared_list .== median_sq_list )
  median_stdvs = sqrt.(csn.attr_median_sum_sq/(csn.count_trials-1) - median_squared_list/csn.count_trials/(csn.count_trials-1))
  coef_vars = csn.attr_coef_var_sum/csn.count_trials
  if csn.N > 1
    coef_var_squared_list = map(x->x^2, csn.attr_coef_var_sum)
    #println("coef_var_squared_list: ",coef_var_squared_list)
    coef_var_sq_list = csn.attr_coef_var_sum .* csn.attr_coef_var_sum
    #println("coef_var_squ_list: ",coef_var_sq_list)
    @assert reduce(&, coef_var_squared_list .== coef_var_sq_list )
    cv_stdevs = sqrt.(csn.attr_coef_var_sum_sq/(csn.count_trials-1) - coef_var_squared_list/csn.count_trials/(csn.count_trials-1))
  else
    cv_stdevs = fill(0.0,num_results)
  end
  for i = 1:num_results
    values = Number[
      csn.N,
      csn.wright_fisher_copy,
      csn.conformist_probability,
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

function coef_var( lst )
  return std(lst)/mean(lst)
end
