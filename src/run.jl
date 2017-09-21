# Example run:  julia -L ContVarEvolution.jl run.jl examples/example1
using ContVarEvolution

function run_trials( simname::AbstractString ) 
  #circular_variation = extreme_variation = false
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  sr = ContVarEvolution.cont_var_result(num_trials,N_list[1],num_subpops,num_attributes_list[1], ngens, burn_in,
      N_mut_list[1]/N_list[1]/100, ideal, wrap_attributes, additive_error, neutral )
  sr_list_run = ContVarEvolution.cont_var_result_type[]
  trial=1
  #=
  for N in N_list
    for num_attributes in num_attributes_list
      for mutation_stddev in mutation_stddev_list
            sr = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
               mutation_stddev, ideal, wrap_attributes, additive_error, neutral )
            Base.push!(sr_list_run, sr )
      end
    end
  end
  =#
  for N in N_list
    for N_mut in N_mut_list
      for num_attributes in num_attributes_list
        for trial = 1:num_trials
          mutation_stddev = N_mut/N
          println("N: ",N,"  N_mut ",N_mut,"  mutation stddev: ",mutation_stddev)
          sr = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
               mutation_stddev, ideal, wrap_attributes, additive_error, neutral )
          Base.push!(sr_list_run, sr )
        end
      end
    end
  end
  println("===================================")
  sr_list_result = pmap(cont_var_simulation, sr_list_run )
  trial = 1
  writeheader( STDOUT, sr )
  writeheader( stream, sr )
  for sr_result in sr_list_result
    writerow(stream,trial,sr_result)
    writerow(STDOUT,trial,sr_result)
    trial += 1
  end
end    

if length(ARGS) == 0
  simname = "configs/example2"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
include("$(simname).jl")
println("simname: ",simname)
#println("simtype: ",simtype)
run_trials( simname )
