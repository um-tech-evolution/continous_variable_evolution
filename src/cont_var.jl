# ContVar structure simulation with horizontal transfer
export cont_var_simulation, fitness
using DataStructures
#=
Recommended command line to run:
>  julia -L ContVarEvolution.jl run_cv.jl configs/example1
=#

@doc """ function cont_var_simulation()
  Wright-Fisher model simulation (as opposed to Moran model)
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m
    ngens number of generations after burn in
    num_attributes   number of quantitative attributes of a variant
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
"""
function cont_var_simulation( sr::ContVarEvolution.cont_var_result_type )
  fit_diff_counter = DataStructures.counter(Int64)
  variant_table = Dict{Int64,variant_type}()
  #int_burn_in = Int(round(sr.burn_in*sr.N+50.0))  # reduce for testing
  id = Int[1]
  n = Int(floor(sr.N/sr.num_subpops))    # size of subpopulations
  println("N: ",sr.N,"  mutation_stddev: ",sr.mutation_stddev,"  num_attributes: ",sr.num_attributes,"  int_burn_in: ",sr.int_burn_in)
  cumm_fitness_means = zeros(Float64,sr.num_subpops)
  cumm_fitness_coef_vars = zeros(Float64,sr.num_subpops)
  cumm_attr_means = [ zeros(Float64,sr.num_attributes) for i in 1:sr.num_subpops]
  cumm_attr_coef_vars = [ zeros(Float64,sr.num_attributes) for i in 1:sr.num_subpops]
  count_gens = 0
  subpops = PopList()
  for j = 1:sr.num_subpops
    Base.push!( subpops, Population() )
    for i = 1:n
      Base.push!( subpops[j], new_innovation( id, sr.ideal, sr.num_attributes, variant_table, sr.neutral ) )
    end
  end
  previous_variant_id = 1
  current_variant_id = id[1]
  previous_subpops = deepcopy(subpops)
  for g = 1:sr.ngens+sr.int_burn_in
    #println("before g: ",g,"  pop: ",subpops[1],"  pop attr: ",[ variant_table[subpops[1][i]].attributes[1] for i = 1:n ])
    after_burn_in = g > sr.int_burn_in
    #println("g: ",g,"  after_burn_in: ",after_burn_in)
    previous_previous_variant_id = previous_variant_id
    previous_variant_id = current_variant_id
    current_variant_id = id[1]
    for j = 1:sr.num_subpops
      for i = 1:n
        cp = copy_parent( previous_subpops[j][i], id, sr.ideal, variant_table, sr, after_burn_in, fit_diff_counter )
        subpops[j][i] = cp
      end
      #println("g: ",g,"  pop: ",subpops[j],"  pop attr: ",[ variant_table[subpops[j][i]].attributes[1] for i = 1:n ])
      #subpops[j] = propsel( subpops[j], n, variant_table )
      if sr.neutral
        new_subpop = deepcopy(subpops[j])
        subpops[j] = [ new_subpop[ rand(1:n) ] for j = 1:n ]
        #println("g: ",g," j: ",j," new_subpop: ",subpops[j])
      else
        subpops[j] = propsel( subpops[j], n, variant_table )
      end
      #println("g: ",g,"  pop: ",subpops[j],"  pop attr: ",[ variant_table[subpops[j][i]].attributes[1] for i = 1:n ])
    end
    previous_subpops = deepcopy(subpops)
    if after_burn_in
      cumm_fitness_means += [ mean( [variant_table[v].fitness for v in s]) for s in subpops]
      cumm_fitness_coef_vars += [ coef_var( [variant_table[v].fitness for v in s]) for s in subpops]
      # cumm_attr_means[s][i] is the mean of attribute i for subpop s, where the mean is over elements of s
      cumm_attr_means += [ [ mean( [ variant_table[v].attributes[i] for v in s]) for i =1:sr.num_attributes ] for s in subpops]
      # cumm_attr_coef_vars[s][i] is the coefficient of variation of attribute i for subpop s, where the mean is over elements of s
      cumm_attr_coef_vars += [ [ coef_var( [ variant_table[v].attributes[i] for v in s]) for i =1:sr.num_attributes ] for s in subpops]
      count_gens += 1
    end
    clean_up_variant_table(previous_previous_variant_id,previous_variant_id,variant_table)
  end  # for g
  @assert count_gens == sr.ngens
  cumm_fitness_means /= sr.ngens
  cumm_fitness_coef_vars /= sr.ngens
  cumm_attr_means /= sr.ngens
  cumm_attr_coef_vars /= sr.ngens
  # The next 2 means are over subpops
  sr.fitness_mean = mean(cumm_fitness_means)
  sr.fitness_coef_var = mean(cumm_fitness_coef_vars)
  # The next 2 means are over attributes and subpops
  sr.attribute_mean = mean(mean(cumm_attr_means))
  sr.attribute_coef_var = mean(mean(cumm_attr_coef_vars))
  (sr.neg_count, sr.neg_neutral, sr.pos_neutral, sr.pos_count ) = summarize_bins( fit_diff_counter )
  return sr
end

function fitness( attributes::Vector{Float64}, ideal::Vector{Float64}, neutral::Bool )
  const fit_slope = 1.0
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  if neutral
    return 1.0
  end
  dis = 0.0
  for k = 1:length(attributes)
    dis += abs( attributes[k] - ideal[k] )
  end
  #result = 1.0-dis/length(attributes)
  result = 1.0/(fit_slope*dis+1.0)
  if result < 0.0
    #println("negative fitness")
    #println("fitness: attributes: ",attributes,"  ideal: ",ideal," fit: ",result)
    result = 0.0
  end
  @assert result >= 0.0
  return result
end

function new_innovation( id::Vector{Int64}, ideal::Float64, num_attributes::Int64, variant_table::Dict{Int64,variant_type}, neutral::Bool )
  i = id[1]
  variant_table[i] = variant_type( 0.0, fill( ideal, num_attributes ) )
  variant_table[i].fitness = fitness( variant_table[i].attributes, fill( ideal, num_attributes), neutral )  
  id[1] += 1
  #println("new innovation attributes: ",variant_table[i].attributes)
  i
end


@doc """  copy_parent()
"""
function copy_parent( v::Int64, id::Vector{Int64}, 
    ideal,
    variant_table::Dict{Int64,ContVarEvolution.variant_type}, 
    sr::ContVarEvolution.cont_var_result_type, after_burn_in::Bool, fit_diff_counter::DataStructures.Accumulator{Int64,Int64} )
  i = id[1]
  vt = variant_table[v]
  new_attributes = mutate_attributes( vt.attributes, sr.mutation_stddev, sr.wrap_attributes, sr.additive_error )
  new_fit = fitness( new_attributes, fill( ideal, sr.num_attributes), sr.neutral )
  if after_burn_in
    increment_bins( fit_diff_counter, new_fit-vt.fitness, 1.0/sr.N )
  end
  variant_table[i] = deepcopy(vt)
  variant_table[i].fitness = new_fit
  variant_table[i].attributes = new_attributes
  id[1] += 1
  return i
end  

function mutate_attributes( attributes::Vector{Float64}, mutation_stddev::Float64, wrap_attributes::Bool, additive_error::Bool )
  #println("mutate attributes  attributes: ",attributes)
  new_attributes = deepcopy(attributes)
  if wrap_attributes
    for i = 1:length(new_attributes)
      #println("B new_attributes[",i,"]: ",new_attributes[i])
      new_attributes[i] += +mutation_stddev*randn()
      if new_attributes[i] < 0
          new_attributes[i] += 1.0
          #println("wrapped up: ",new_attributes[i])
      end
      if new_attributes[i] > 1.0
          new_attributes[i] -= 1.0
          #println("wrapped down: ",new_attributes[i])
      end
      new_attributes[i] = min(1.0,max(0.0,new_attributes[i]))
      #println("A new_attributes[",i,"]: ",new_attributes[i])
    end
    #println("new_attributes: ",new_attributes)
    return new_attributes
  else
    if additive_error
      for i = 1:length(new_attributes)
        new_attributes[i] += +mutation_stddev*randn()
      end
    else
      for i = 1:length(new_attributes)
        if new_attributes[i] <= 0.0
          println("neg attribute: ",new_attributes[i])
          new_attributes[i] = 1.0e-6
        end
        multiplier = (1.0+mutation_stddev*randn())
        while multiplier <= 1.0e-6
          println("neg multiplier")
          multiplier = (1.0+mutation_stddev*randn())
        end
        new_attributes[i] *= multiplier
        if new_attributes[i] < 0.0
          println("negative attribute with i=",i,": attribute: ",new_attributes[i])
        end
      end
    end
    return new_attributes
  end
end

function print_subpop( subpop::Vector{Int64}, variant_table::Dict{Int64,variant_type} )
  "[ $(join([ @sprintf(" %5d:%5.4f",vt,variant_table[vt].fitness) for vt in subpop ]))]"
end

function print_pop( stream::IO, subpops::PopList, variant_table::Dict{Int64,variant_type} )
  for sp in subpops
    print(stream,print_subpop(sp,variant_table))
  end
  println(stream)
end

using DataFrames
# compute and save statistics about subpopulations and populations

function means( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_means( subpops::PopList, variant_table::Dict{Int64,variant_type}, num_attributes::Int64 )
  #num_attributes = length(variant_table[1].attributes)
  #println("attr_vars: num_attributes: ",num_attributes)
  ave_means = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_means = [ mean([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_means: ",att_means)
    ave_means[i] = mean(att_means)
    i += 1
  end
  println("ave_means: ",ave_means)
  return ave_means
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type}, num_attributes::Int64 )
  #num_attributes = length(variant_table[1].attributes)
  #println("attr_vars: num_attributes: ",num_attributes)
  ave_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_vars = [ var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    #println(s," att_vars: ",att_vars)
    ave_vars[i] = mean(att_vars)
    i += 1
  end
  #println("ave_vars: ",ave_vars)
  return ave_vars
end

function fit_loc_index(N,num_subpops,num_fit_locs,j,i)
  #=
  n = Int(ceil(N/num_subpops))
  mult = Int(ceil(num_fit_locs/num_subpops))
  div = Int(ceil(n*num_subpops/num_fit_locs))
  return mult*(j-1) + Int(floor((i-1)/div))+1
  =#
  return 1
end

function coef_var( lst )
  return std(lst)/mean(lst)
end

#= Duplicates a function in bin_data.jl
# Create a dictionary that bins a real-valued vector
function create_bins( vect::Vector{Float64}, cutoff::Float64 )
  bins = DataStructures.counter(Int64)
  for v in vect
    index = Int(floor(v/cutoff))
    push!(bins, index)
  end
  bins
end
=#

function clean_up_variant_table( previous_variant_id::Int64, previous_previous_variant_id::Int64,
    variant_table::Dict{Int64,variant_type} )
  #println("clean up:  ppv: ",previous_previous_variant_id,"  pv: ",previous_variant_id)
  for v = previous_variant_id:previous_previous_variant_id-1
    delete!(variant_table,v)
  end
end
 

