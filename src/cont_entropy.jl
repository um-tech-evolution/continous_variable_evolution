# A function to compute the entropy of a population of continous variables which decreases entropy when population members are close.
# The w parameter defines what "close" means.  
# For a population point x, if there are other points in the population which are closer than w to x, then entropy is reduced.
# For a population where all points are distance at least w from each other, this definition reduces to discrete entropy.
# This applies if the population has exact duplicates.
export cont_entropy

function cont_entropy( x::Vector{Float64}, w::Float64=0.25 )
  n = length(x)
  ent_sum = 0.0
  for i = 1:n
    d = 0.0
    for j = 1:n
      if i != j && abs(x[j] - x[i]) < w
        d += (w - abs(x[j] - x[i]))/w
        #println("(i,j): ",(i,j)," dinc: ",(w - abs(x[j] - x[i]))/w,"  d: ",d)
      end
    end
    #println("i: ",i,"  d: ",d,"  ent_inc: ",log2(n/(d+1)) / n)
    ent_sum += log2(n/(d+1))/n 
  end
  ent_sum
end

