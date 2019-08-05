# My definition of entropy for a population C = [c_1,...,c_n] where each c_i is a real value.  
# The technique is to use a step function for each c_i:
#   f_{c_i}(x) = abs(x-c_i)/2 < w ? 1.0/w : 0.0
# More precisely, we assume that all c_i are in an interval [a,b]
#   f_{c_i}(x) = 1/(u-v) ? v <= x < u : 0 
#        where v = max(a,c_i-w/2) and u = min(b,c_i+w/2).
# This function is a probability distribution.
# For the populaition,  f_C(x) = 1/n \sum_{i=1}^n f_{c_i}(x)
export cont_entropy


w = 0.1   # for testing

function step(x::Float64, c::Float64, w::Float64, a::Float64=0.0, b::Float64=1.0) 
  v = max(a,c-w/2)
  u = min(b,c+w/2)
  return x >= v && x < u  ? 1.0/(u-v) : 0.0
end

function steps(x::Float64, C::Vector{Float64}, w::Float64, a::Float64=0.0, b::Float64=1.0)
  if x < a || b < x
    return 0.0
  end
  sum( step(x,c,w) for c in C )/length(C)
end

function truncate(x::Float64, a::Float64, b::Float64)
  result1 = x >= a ? x : a
  return result1 <= b ? result1 : b
end

# Returns a 2D array whose rows are the lower and upper bounds of the intervals correspond
#   to the elements of C.  These boundes are truncated to be between a and b
function steps_bounds( C::Vector{Float64}, w::Float64, a::Float64, b::Float64)
  n = length(C)
  result = zeros(Float64,n,3)
  i = 1
  for c in C
    lb = truncate(c-0.5*w,a,b)
    ub = truncate(c+0.5*w,a,b)
    val = 1.0/(ub-lb)/n
    result[i,:] = [lb,ub,val]
    #result[i,:] = [truncate(c-0.5*w,a,b), truncate(c+0.5*w,a,b)]
    #println("i: ",i,"  ",result[i,:])
    i += 1
  end
  result
end 

function bounds_list( step_bds::Array{Float64,2}, a::Float64, b::Float64 )
  n = size(step_bds)[1]
  interval_bds = 
    sort(vcat([(step_bds[i,1],:lower,step_bds[i,3]) for i=1:n],
              [(step_bds[i,2],:upper,step_bds[i,3]) for i=1:n]), by=x->x[1])
  bounds_list = vcat((a,:global,0.0),interval_bds,(b,:global,0.0))
  bounds_list
end
  
function values_list( bounds_lst::Array{Tuple{Float64,Symbol,Float64},1} )
  N = size(bounds_lst)[1]
  n = div(N-2,2)
  val = 0.0
  d = fill((0.0,0.0),N-1)
  sum = 0.0
  sum_entropy = 0.0
  for i = 1:2*n+1
    #if bounds_lst[i] == bounds_lst[i+1] break end
    if bounds_lst[i][2] == :lower
      val += bounds_lst[i][3]
    elseif bounds_lst[i][2] == :upper
      val -= bounds_lst[i][3]
    end
    d[i] = (bounds_lst[i+1][1],val)
    suminc = (bounds_lst[i+1][1]-bounds_lst[i][1])*d[i][2]
    sum += suminc    
    sum_entropy -= suminc > 0.0 ? suminc * log2(suminc) : 0.0
    #print("i: ",i,"  bl[i+1]: ",bounds_lst[i+1],"  d[i]: ",d[i],"  suminc: ",suminc)
    #println("  suminc*log2(suminc): ",suminc*log2(suminc))
  end
  #println("sum: ",sum,"  sum_entropy: ",sum_entropy)
  #d
  sum_entropy
end

function cont_entropy( C::Vector{Float64}, w::Float64, a::Float64=0.0, b::Float64=1.0 )
  c_ent = values_list(bounds_list(steps_bounds(C,w,a,b),a,b))
  #println("cont_entropy: C: ",C,"  w: ",w,"  (a,b):",(a,b),"  c_ent: ",c_ent)
  return c_ent
end

function step_bounds_check( sb )
  n = size(sb)[1]
  sum = 0.0
  for i = 1:n
    sum += (sb[i,2]-sb[i,1])*sb[i,3]
  end
  sum
end
