# ContVar structure simulation with horizontal transfer
export cont_var_simulation, fitness
try    # "using Printf" fails in julia v6, but avoids a deprecation error in julia v7
  using Printf
catch
end
#using DataStructures
import DataStructures.counter
import DataStructures.Accumulator
#=
Recommended command line to run:
>  julia run.jl examples/example1
=#

@doc """ function cont_var_simulation()
  Agent-based Henrich 2004 model
  Parameters:
    N     MetaPopulation size
    m     number of subpopulations   # for now, subpopulation size = N/m;  m=1 for Henrich model
    ngens number of generations after burn in
    num_attributes   number of quantitative attributes of a variant;  num_attributes=1 for Henrich model
    variant_table Keeps track fitnesses and variant parent and innovation ancestor
"""
function henrich_simulation( simrecord::ContVarEvolution.cont_var_result_type )
  fit_diff_counter = counter(Int64)
  variant_table = Dict{Int64,variant_type}()
  #int_burn_in = Int(round(simrecord.burn_in*simrecord.N+50.0))  # reduce for testing
  id = Int[1]
  n = Int(floor(simrecord.N/simrecord.num_subpops))    # size of subpopulations
  #println("N: ",simrecord.N,"  mutation_stddev: ",simrecord.mutation_stddev,"  num_attributes: ",simrecord.num_attributes,"  int_burn_in: ",simrecord.int_burn_in)
  cumm_fitness_means = zeros(Float64,simrecord.num_subpops)
  cumm_fitness_medians = zeros(Float64,simrecord.num_subpops)
  cumm_fitness_coef_vars = zeros(Float64,simrecord.num_subpops)
  cumm_attr_means = [ zeros(Float64,simrecord.num_attributes) for i in 1:simrecord.num_subpops]
  cumm_attr_medians = [ zeros(Float64,simrecord.num_attributes) for i in 1:simrecord.num_subpops]
  cumm_attr_coef_vars = [ zeros(Float64,simrecord.num_attributes) for i in 1:simrecord.num_subpops]
  count_gens = 0
  subpops = PopList()
  for j = 1:simrecord.num_subpops
    Base.push!( subpops, Population() )
    for i = 1:n
      Base.push!( subpops[j], new_innovation( id, simrecord.ideal, simrecord.num_attributes, variant_table, simrecord.neutral ) )
    end
  end
  previous_variant_id = 1
  current_variant_id = id[1]
  previous_subpops = deepcopy(subpops)
  for g = 1:simrecord.ngens+simrecord.int_burn_in
    #println("before g: ",g,"  pop: ",subpops[1],"  pop attr: ",[ variant_table[subpops[1][i]].attributes[1] for i = 1:n ])
    after_burn_in = g > simrecord.int_burn_in
    previous_previous_variant_id = previous_variant_id
    previous_variant_id = current_variant_id
    current_variant_id = id[1]
    for j = 1:simrecord.num_subpops
      for i = 1:n
        subpops[j][i] = mutate_attributes( previous_subpops[j][i], id, variant_table, simrecord, after_burn_in, fit_diff_counter )
      end
      #println("B fits: ",[variant_table[subpops[j][i]].attributes[1] for i = 1:n ])
      if simrecord.neutral
        r = rand(1:n,n)
        subpops[j] = subpops[j][ r ]
      else
        subpops[j] = power_sel( subpops[j], n, simrecord.fit_power, variant_table )
      end
      if simrecord.renormalize    # rescale fitness so that the maximum is 1.0
        renormalize( subpops[j], n, variant_table )
      end
      #println("A fits: ",[variant_table[subpops[j][i]].attributes[1] for i = 1:n ])
    end
    previous_subpops = deepcopy(subpops)
    #println("g: ",g,"  fitness: ", [ mean( [variant_table[v].fitness for v in s]) for s in subpops][1] )
    if after_burn_in
      cumm_fitness_means += [ mean( [variant_table[v].fitness for v in s]) for s in subpops]
      cumm_fitness_medians += [ median( [variant_table[v].fitness for v in s]) for s in subpops]
      cumm_fitness_coef_vars += [ coef_var( [variant_table[v].fitness for v in s]) for s in subpops]
      # cumm_attr_means[s][i] is the mean of attribute i for subpop s, where the mean is over elements of s
      cumm_attr_means += [ [ mean( [ variant_table[v].attributes[i] for v in s]) for i =1:simrecord.num_attributes ] for s in subpops]
      # cumm_attr_means[s][i] is the median of attribute i for subpop s, where the median is over elements of s
      cumm_attr_medians += [ [ median( [ variant_table[v].attributes[i] for v in s]) for i =1:simrecord.num_attributes ] for s in subpops]
      # cumm_attr_coef_vars[s][i] is the coefficient of variation of attribute i for subpop s, where the mean is over elements of s
      cumm_attr_coef_vars += [ [ coef_var( [ variant_table[v].attributes[i] for v in s]) for i =1:simrecord.num_attributes ] for s in subpops]
      #println("fitness_mean: ", [ mean( [variant_table[v].fitness for v in s]) for s in subpops])
      #println("attr_coef_var: ",  [ [ coef_var( [ variant_table[v].attributes[i] for v in s]) for i =1:simrecord.num_attributes ] for s in subpops])
      #println("attr mean: ", [ [ mean( [ variant_table[v].attributes[i] for v in s]) for i =1:simrecord.num_attributes ] for s in subpops])
      count_gens += 1
    end
    clean_up_variant_table(previous_previous_variant_id,previous_variant_id,variant_table)
  end  # for g
  @assert count_gens == simrecord.ngens
  cumm_fitness_means /= simrecord.ngens
  cumm_fitness_medians /= simrecord.ngens
  cumm_fitness_coef_vars /= simrecord.ngens
  cumm_attr_means /= simrecord.ngens
  cumm_attr_medians /= simrecord.ngens
  cumm_attr_coef_vars /= simrecord.ngens
  # The next 3 means are over subpops
  simrecord.fitness_mean = mean(cumm_fitness_means)
  simrecord.fitness_median = mean(cumm_fitness_medians)
  simrecord.fitness_coef_var = mean(cumm_fitness_coef_vars)
  # The next 3 means are over attributes and subpops
  simrecord.attribute_mean = mean(mean(cumm_attr_means))
  simrecord.attribute_median = mean(mean(cumm_attr_medians))
  simrecord.attribute_coef_var = mean(mean(cumm_attr_coef_vars))
  (simrecord.neg_count, simrecord.neg_neutral, simrecord.pos_neutral, simrecord.pos_count ) = summarize_bins( fit_diff_counter )
  return simrecord
end

function hfitness( attributes::Vector{Float64}, ideal::Vector{Float64}, neutral::Bool )
  if length(attributes) != length(ideal)
    error("length(attributes) must equal length(ideal) in fitness")
  end
  if neutral
    return 1.0
  end
  return maximum(attributes)
end

@doc """ function new_innovation()
  Create a new variant_table record whose attributes are all set to ideal (specificied as an input parameter)
"""
function new_innovation( id::Vector{Int64}, ideal::Float64, num_attributes::Int64, variant_table::Dict{Int64,variant_type}, neutral::Bool )
  i = id[1]
  variant_table[i] = variant_type( 0.0, fill( ideal, num_attributes ) )
  variant_table[i].fitness = hfitness( variant_table[i].attributes, fill( ideal, num_attributes), neutral )  
  id[1] += 1
  i
end


@doc """ function mutate_attributes()
  Mutate the attributes corresponding the variant_table[v].
  The mutation is actually done by the function mutate().
"""
function mutate_attributes( v::Int64, id::Vector{Int64}, variant_table::Dict{Int64,ContVarEvolution.variant_type}, 
    simrecord::ContVarEvolution.cont_var_result_type, after_burn_in::Bool, fit_diff_counter::Accumulator{Int64,Int64} )
  i = id[1]
  vt = variant_table[v]
  new_attributes = mutate( vt.attributes, simrecord.mutation_stddev, simrecord.mutation_bias )
  new_fit = hfitness( new_attributes, fill( simrecord.ideal, simrecord.num_attributes), simrecord.neutral )
  if after_burn_in
    increment_bins( fit_diff_counter, new_fit-vt.fitness, 1.0/simrecord.N )
  end
  variant_table[i] = deepcopy(vt)
  variant_table[i].fitness = new_fit
  variant_table[i].attributes = new_attributes
  id[1] += 1
  return i
end  

@doc """ function mutate()
  Multiplicatively mutate the attributes array.
  mutation_bias  should be positive and less than 1.0
"""
function mutate( attributes::Vector{Float64}, mutation_stddev::Float64, 
     mutation_bias::Float64  )
  new_attributes = deepcopy(attributes)
      for i = 1:length(new_attributes)
        if new_attributes[i] <= 0.0
          println("neg attribute: ",new_attributes[i])
          new_attributes[i] = 0.1
        end
        multiplier = (mutation_bias+mutation_stddev*randn())
        #println("mutation_bias: ",mutation_bias,"  multiplier: ",multiplier)
        while multiplier <= 0.1
          #println("neg multiplier")
          multiplier = (mutation_bias+mutation_stddev*randn())
        end
        new_attributes[i] *= multiplier
        if new_attributes[i] < 0.0
          println("negative attribute with i=",i,": attribute: ",new_attributes[i])
        end
      end
    return new_attributes
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

# compute and save statistics about subpopulations and populations

function means( subpops::PopList, variant_table::Dict{Int64,variant_type} )
  fit(v) = variant_table[v].fitness
  means = [ mean(map(fit,s)) for s in subpops ]
  vars = [ var(map(fit,s)) for s in subpops ]
  return means, vars
end

function attr_means( subpops::PopList, variant_table::Dict{Int64,variant_type}, num_attributes::Int64 )
  ave_means = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_means = [ mean([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    ave_means[i] = mean(att_means)
    i += 1
  end
  #println("ave_means: ",ave_means)
  return ave_means
end

function attr_vars( subpops::PopList, variant_table::Dict{Int64,variant_type}, num_attributes::Int64 )
  ave_vars = zeros(Float64,length(subpops))
  i = 1
  for s in subpops
    att_vars = [ var([variant_table[v].attributes[j] for v in s]) for j =1:num_attributes]
    ave_vars[i] = mean(att_vars)
    i += 1
  end
  return ave_vars
end

function coef_var( lst )
  return StatsBase.std(lst)/mean(lst)
end

function clean_up_variant_table( previous_variant_id::Int64, previous_previous_variant_id::Int64,
    variant_table::Dict{Int64,variant_type} )
  #println("clean up:  ppv: ",previous_previous_variant_id,"  pv: ",previous_variant_id)
  for v = previous_variant_id:previous_previous_variant_id-1
    delete!(variant_table,v)
  end
end
 

