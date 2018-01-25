# Create a dictionary that bins a real-valued vector
# Used to compute fit-diff statistics
#  Example:
#  julia> w = [-2.,-1,0.,1.]
#  
#  4-element Array{Float64,1}:
#   -2.0
#   -1.0
#    0.0
#    1.0
#  
#  julia> wbins = create_bins(w,1.0)
#  DataStructures.Accumulator{Int64,Int64} with 4 entries:
#    0  => 1
#    -2 => 1
#    -1 => 1
#    1  => 1
#  
#  julia> summarize_bins(wbins)
#  (1,1,1,1)

import DataStructures.Accumulator
import DataStructures.counter
export create_bins, increment_bins, summarize_bins

@doc """ function create_bins()
 Starting with a list of floats, count the number that are in the intervals [(i-1)/coutff,i/cutoff] for all integers i.
"""
function create_bins( vect::Vector{Float64}, cutoff::Float64 )
  bins = counter(Int64)
  for v in vect
    index = Int(floor(v/cutoff))
    push!(bins, index)
  end
  bins
end

@doc """ function increment_bins()
  Add value  x  to bins. 
  The cutoff value should be the value used to create the bins.
"""
function increment_bins( bins::Accumulator{Int64,Int64}, x::Float64, cutoff::Float64 )
  index = Int(floor(x/cutoff))
  push!(bins, index)
end  

@doc """ function create_bins_by_filter()
 Create bins by using the Julia filter functon.
 Used to check the correctness of create_bins()
"""
function create_bins_by_filter( vect::Vector{Float64}, cutoff::Float64 )
  bins = Dict{Int64,Int64}()
  Min = minimum(vect)
  Max = maximum(vect)
  for i = collect(Int(floor(Min/cutoff)):(Int(floor(Max/cutoff))))
    len = length(filter(x->(((i*cutoff))<=x && x <(i+1)*cutoff),vect))
    println("i: ",i,"  lb: ",i*cutoff,"  ub: ",(i+1)*cutoff,"  len: ",len)
    #println("  len: ", length(filter(x->(((i*cutoff))<=x && x <(i+1)/cutoff),vect)) )
    if len > 0
      bins[i] = len
    end
  end
  bins
end

@doc """ check_bins()
  Creates bins dictionary, and then uses create_bins_by_filter to check correctness
"""
function check_bins( vect::Vector{Float64}, cutoff::Float64 )
  bins = create_bins( vect, cutoff )
  fbins =  create_bins_by_filter( vect, cutoff )
  for k in keys(bins)
    @assert bins[k] == fbins[k]
  end
  for k in keys(fbins)
    @assert bins[k] == fbins[k]
  end
end

@doc """ function summarize_bins()
  Returns a 4-tuple of vect values in the following intervals: 
  (-infinity,-cutoff)
  [-cutoff,0.0)
  [0.0,cutoff)
  [cutoff,infinity)
"""
function summarize_bins( bins::Accumulator{Int64,Int64} )
  min = minimum(keys(bins))
  max = maximum(keys(bins))
  neg_neutral = bins[-1]
  pos_neutral = bins[0]
  neg_count = min <= -2 ? sum( bins[k] for k = min:-2 ) : 0
  pos_count = 1 <= max ? sum( bins[k] for k = 1:max ) : 0
  return (neg_count,neg_neutral,pos_neutral,pos_count)
end
