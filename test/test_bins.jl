using DataFrames, Dates, CSV
import DataStructures.counter
import DataStructures.Accumulator

include("../../continuous_variable_evolution/src/types.jl")
include("../../continuous_variable_evolution/src/cont_var.jl")
include("../../continuous_variable_evolution/src/run_cont_var.jl")
include("../../continuous_variable_evolution/src/bin_data.jl")
include("../../continuous_variable_evolution/src/propsel.jl")

function test_bins()
  fdc = counter(Int64)
  cutoff = 0.125
  vals= 3.0*rand(20) .- 1.5
  for v in vals
    increment_bins( fdc, v, cutoff )
  end
  println("bins: ",fdc)
  bins_to_vector( fdc, cutoff )
end

function test_run_simulation(mutation_stddev::Float64, N::Int64, ngens::Int64, bins_cutoff::Float64;
    csvfile::String="" )
  #N = 100; 
  N_list=[N]
  num_attributes=1; num_attributes_list = [num_attributes]
  #mutation_stddev=0.10;  
  mutation_stddev_list = [mutation_stddev]
  N_mut_list = []
  num_trials=1
  num_subpops=1
  num_attributes=1
  #ngens=10
  ideal=1.0
  fit_slope=1.0
  int_burn_in=100
  neutral=:false
  #bins_cutoff = 0.02
  simrecord = cont_var_result_type( N_list, num_attributes_list, mutation_stddev_list, N_mut_list, num_trials, 
      N, num_subpops, num_attributes, ngens, int_burn_in,
      mutation_stddev, ideal, fit_slope, neutral, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0,0,0,0 )
  (simrec,bins,lbvec)=cont_var_simulation( simrecord, bins_to_vect=:true, bins_cutoff=bins_cutoff )
  println("bins:  ",bins)
  println("lvvec: ",lbvec)
  df = DataFrame()
  df.bins=bins
  df.lbvec=lbvec
  if length(csvfile) > 0
    open( csvfile, "w" ) do f
      hostname = chomp(open("/etc/hostname") do f read(f,String) end)
      println(f,"# date and time: ",Dates.now())
      println(f,"# host: ",hostname," with ",nprocs()-1,"  processes: " )
      println(f,"# N: ",N)
      println(f,"# mutation_stddev: ",mutation_stddev)
      println(f,"# ngens: ",ngens)
      println(f,"# bins_cutoff: ",bins_cutoff)
      CSV.write( f, df, append=true, writeheader=true )
    end
  end
  (bins, lbvec)
end

