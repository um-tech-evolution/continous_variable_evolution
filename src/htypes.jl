#  Defines the types "Population" and "PopList".
#  Defines the composite types "cont_var_result_type" and "variant_type"
export variant_type, fitness_location_type
const Population = Array{Int64,1}
const PopList = Array{Population,1}
using StatsBase
try  # for julia v7
  using Distributed
  using Random
  using Dates
catch
end


mutable struct variant_type
  fitness::Float64    # The fitness of this variant
  attributes::Vector{Float64}   # attributes of the variant
end

mutable struct cont_var_result_type
  N_list::Vector{Int64}
  #num_attributes_list::Vector{Int64}
  mutation_stddev_list::Vector{Float64}  # list of mutation_stddev values
  N_mut_list::Vector{Float64}  # list of mutation_stddev values
  num_trials::Int64  # Number of times to repeat simulation for each setting of the parameters
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  num_attributes::Int64  # number of attributes of a variant
  ngens::Int64  # number of generations after burn-in
  int_burn_in::Int64
  mutation_stddev::Float64  # standard deviation of mutation distribution of mutation perturbations
  mutation_bias::Float64  # Multiplicate bias of mutation. Must be positve.  Normally < 1.0
  ideal::Float64             #ideal value 
  fit_power::Float64         # use proportional selection applied to fitness^fit_power.
  renormalize::Bool          # multiplicatively renormalize fitesses so that max fitness is 1.0 on every generation
  neutral::Bool              # If true, fitness = 1
  fitness_mean::Float64      # average of fitnesses over subpopulations and generations
  fitness_median::Float64      # average of fitnesses over subpopulations and generations
  fitness_coef_var::Float64  # average coefficent of variation of fitnesss over subpopulations and generations
  attribute_mean::Float64    # average of attributes over attributes, subpopulations, and generations
  attribute_median::Float64    # average of attributes over attributes, subpopulations, and generations
  attribute_coef_var::Float64 # average coefficent of variation of attributes over attributes, subpopulations and generations
  neg_count::Int64        # Number of fitness differences < -1.0/N
  neg_neutral::Int64      # Number of fitness differences >= -1.0/N  and < 0.0
  pos_neutral::Int64     # Number of fitness differences >= 0.0 and < 1.0/N
  pos_count::Int64        # Number of fitness differences >= 1.0/N
end

