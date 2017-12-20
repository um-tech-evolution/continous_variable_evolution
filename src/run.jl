# Example run:  julia -L ContVarEvolution.jl run.jl examples/example1
# Example run:  julia -p :4 -L ContVarEvolution.jl run.jl examples/example1
using ContVarEvolution

function run_trials( simname::AbstractString ) 
  global mutation_stddev_list
  global N_mut_list
  #circular_variation = extreme_variation = false
  stream = open("$(simname).csv","w")
  println("stream: ",stream)
  println("isdefined mutation_stddev_list: ",isdefined(:mutation_stddev_list))
  println("isdefined N_mut_list: ",isdefined(:N_mut_list))
  if isdefined(:mutation_stddev_list)
    sr = ContVarEvolution.cont_var_result(num_trials,N_list[1],num_subpops,num_attributes_list[1], ngens, burn_in,
       mutation_stddev_list[1], ideal, fit_slope, wrap_attributes, additive_error, neutral )
  elseif isdefined(:N_mut_list)
    sr = ContVarEvolution.cont_var_result(num_trials,N_list[1],num_subpops,num_attributes_list[1], ngens, burn_in,
       N_mut_list[1]/N_list[1]/100, ideal, fit_slope, wrap_attributes, additive_error, neutral )
  end
  sr_list_run = ContVarEvolution.cont_var_result_type[]
  trial=1
  if isdefined(:mutation_stddev_list)
    for N in N_list
      for num_attributes in num_attributes_list
        for mutation_stddev in mutation_stddev_list
          for trial = 1:num_trials
              sr = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
                 mutation_stddev, ideal, fit_slope, wrap_attributes, additive_error, neutral )
              Base.push!(sr_list_run, sr )
          end
        end
      end
    end
  elseif isdefined(:N_mut_list)
    for N in N_list
      for N_mut in N_mut_list
        for num_attributes in num_attributes_list
          for trial = 1:num_trials
            mutation_stddev = N_mut/N
            sr = ContVarEvolution.cont_var_result(num_trials,N,num_subpops,num_attributes, ngens, burn_in,
                 mutation_stddev, ideal, fit_slope, wrap_attributes, additive_error, neutral )
            int_burn_in = Int(round(burn_in*sr.N+50.0)) 
            #println("N: ",N,"  N_mut ",N_mut,"  mutation stddev: ",mutation_stddev,"  int_burn_in: ",int_burn_in)
            Base.push!(sr_list_run, sr )
          end
        end
      end
    end
  end
  println("===================================")
  sr_list_result = map(cont_var_simulation, sr_list_run )
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
  simname = "examples/example2"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
include("$(simname).jl")
#println("simname: ",simname)
println("simtype: ",simtype)
println("num_trials: ",num_trials)
run_trials( simname )
