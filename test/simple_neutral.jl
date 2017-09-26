# Simple Wright-Fisher model for neutral continuous variable evolution
# Command line:  julia simple_neutral.jl
# Command line:  julia simple_neutral.jl <seed>
#    where seed is a random number seed.
# When run with a specific seed, the output should agree running the full simulation
#    for neutral multiplicative error with a single attribute (with the same seed).

N=10
mutstddev = 0.04
burn_in = 3
ngens = 1
ideal = 1.0
num_trials = 2

type simple_neutral_type
  N::Int64
  mutstddev::Float64
  ngens::Int64
  ideal::Float64
  int_burn_in::Int64
  average_attr_mean::Float64
  average_attr_coef_var::Float64
end

function print_simple_neutral( sn::simple_neutral_type )
  println("N: ",sn.N,"  mutstddev: ",sn.mutstddev," ngens: ",sn.ngens," int_burn_in: ",sn.int_burn_in)
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
  cumm_attr_mean = 0.0
  cumm_attr_coef_var = 0.0
  sum_gens = 0
  pop = fill( sn.ideal, sn.N )
  for g = 1:(sn.int_burn_in + sn.ngens )
    i = 1
    #println("before copy g: ",g,"  pop: ",pop)
    new_pop = copy_pop( sn, pop )
    #println("after copy g: ",g,"  new_pop: ",new_pop)
    pop = [ new_pop[ rand(1:sn.N) ] for j = 1:sn.N ]
    #println("after WF g: ",g,"  pop: ",pop)
    if g > sn.int_burn_in
      cumm_attr_mean += mean(pop)
      cumm_attr_coef_var += coef_var(pop)
      sum_gens += 1
      #println("cumm_attr_mean: ",cumm_attr_mean,"cumm_attr_coef_var: ",cumm_attr_coef_var)
    end
    #pop = deepcopy(new_pop)
  end # for g
  @assert sn.ngens == sum_gens
  cumm_attr_mean /= ngens
  cumm_attr_coef_var /= ngens
  #println("cumm_attr_mean: ",cumm_attr_mean)
  #println("cumm_attr_coef_var: ",cumm_attr_coef_var)
  sn.average_attr_mean = cumm_attr_mean
  sn.average_attr_coef_var = cumm_attr_coef_var
  return sn
end

function run_trials( num_trials::Int64, N::Int64, mutstddev::Float64, ngens::Int64, burn_in::Number )
  overall_avg_mean = 0.0
  overall_avg_coef_var = 0.0
  sn = simple_neutral_type( N, mutstddev, ngens, ideal, 0, 0.0, 0.0 )
  if typeof(burn_in) == Int64
    sn.int_burn_in = burn_in
  else
    sn.int_burn_in = Int(round(burn_in*sn.N+50.0)) 
  end
  print_simple_neutral( sn )
  for t = 1:num_trials
    sn = simple_neutral_simulation( sn )
    println("trial: ",t,"  avg mean: ",sn.average_attr_mean,"  avg coef var: ",sn.average_attr_coef_var )
    overall_avg_mean += sn.average_attr_mean
    overall_avg_coef_var += sn.average_attr_coef_var
  end
  overall_avg_mean /= num_trials
  overall_avg_coef_var /= num_trials
  println("overall avg mean: ",overall_avg_mean,"  overall avg coef_var: ",overall_avg_coef_var)
end

function coef_var( lst )
  return std(lst)/mean(lst)
end

if length(ARGS) >= 1   # second command-line argument is random number seed
  seed = parse(Int,ARGS[1])
  println("seed: ",seed)
  srand(seed)
end

run_trials( num_trials, N, mutstddev, ngens, burn_in )
