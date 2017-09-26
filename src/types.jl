export variant_type, fitness_location_type
#using Distributions
const Population = Array{Int64,1}
const PopList = Array{Population,1}

type variant_type
  fitness::Float64    # The fitness of this variant
  attributes::Vector{Float64}   # attributes of the variant
end

type cont_var_result_type
  num_trials::Int64  # Number of times to repeat simulation for each setting of the parameters
  N::Int64   # meta-population size
  num_subpops::Int64   # number of subpopulations
  num_attributes::Int64  # number of attributes of a variant
  #mu::Float64     # innovation rate
  ngens::Int64  # number of generations after burn-in
  int_burn_in::Int64
  mutation_stddev::Float64  # standard deviation of mutation distribution of mutation perturbations
  ideal::Float64      #maximum ideal value 
  wrap_attributes::Bool      # Keep attributes within the unit interval by wrapping around (so that attribute space is toroidal)
  additive_error::Bool       # If true, use additive rather than multiplicative copy error
  neutral::Bool              # If true, fitness = 1
  fitness_mean::Float64      # average of fitnesses over subpopulations and generations
  fitness_coef_var::Float64  # average coefficent of variation of fitnesss over subpopulations and generations
  attribute_mean::Float64    # average of attributes over attributes, subpopulations, and generations
  attribute_coef_var::Float64 # average coefficent of variation of attributes over attributes, subpopulations and generations
  neg_count::Int64        # Number of fitness differences < -1.0/N
  neg_neutral::Int64      # Number of fitness differences >= -1.0/N  and < 0.0
  pos_neutral::Int64     # Number of fitness differences >= 0.0 and < 1.0/N
  pos_count::Int64        # Number of fitness differences >= 1.0/N
end

