# Create a dictionary that bins a real-valued vector
import DataStructures.Accumulator
import DataStructures.counter
export create_bins, increment_bins, summarize_bins

# Starting with a list of floats, count the number that are in the intervals [(i-1)/coutff,i/cutoff] for all integers i.
function create_bins( vect::Vector{Float64}, cutoff::Float64 )
  bins = counter(Int64)
  for v in vect
    index = Int(floor(v/cutoff))
    push!(bins, index)
  end
  bins
end

function increment_bins( bins::Accumulator{Int64,Int64}, x::Float64, cutoff::Float64 )
  index = Int(floor(x/cutoff))
  push!(bins, index)
end  

# Create bins by using the Julia filter functon.
# Used to check the correctness of create_bins()
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

function summarize_bins( bins::Accumulator{Int64,Int64} )
  min = minimum(keys(bins))
  max = maximum(keys(bins))
  neg_neutral = bins[-1]
  pos_neutral = bins[0]
  neg_count = min <= -2 ? sum( bins[k] for k = min:-2 ) : 0
  pos_count = 1 <= max ? sum( bins[k] for k = 1:max ) : 0
  return (neg_count,neg_neutral,pos_neutral,pos_count)
end
