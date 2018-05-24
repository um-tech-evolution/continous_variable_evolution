# Proportional selection using the stohcastic acceptance algorithm of Lipowski and Lipowska

@doc """ propsel_funct()
Proportional selection applied to a population of attribute values.
Fitness is given by: function  fitness( x::Real64 )
"""
function propsel_funct( attributes::Vector{Float64}, fitness::Function, initial_value::Float64, fit_slope::Float64 )
  #println("propsel funct: initial value: ",initial_value,"  fit_slope: ",fit_slope)
  #println("  fit: ",fitness(attributes[1],initial_value,fit_slope))
  #println("propsel funct: pop: ",attributes)
  n = length(attributes)
  fitvec = map( a->fitness(a,initial_value,fit_slope), attributes )
  #println("fitvec: ",fitvec)
  fmax = maximum(fitvec)
  if fmax == 0.0
    # all elements have fitness zero
    println("Warning: all elements have fitness zero")
    return deepcopy(attributes)  # do not modify attributes
  end
  #new_attributes = zeros(Float64,n)
  fitvec /= fmax
  #println("fitvec: ",fitvec)
  selected = zeros(Int64, n)
  k = 0
  while k < n
    i = rand(1:n)
    w = fitvec[i]
    if rand() < w
      selected[k + 1] = i
      k += 1
    end
  end
  attributes[selected] 
end

@doc """ test_propsel_funct()
  Does a simple test of the propsel_funct function.
  Will fail with an assertion error if the test fails.
"""
function test_propsel_funct()
  n = 10
  # create an attribute array [1.0, 2.0, . . . ,10.0]
  attr = [ convert(Float64,i) for i=1:n]
  # x->(x==5.0 ? :2.0 : 0.0)  is an anonymous function that gives only 5.0 fitness 2.0
  #   and all other numbers fitness 0.0
  new_attr = propsel_funct( attr, x->(x==5.0 ? :2.0 : 0.0) )
  @assert( new_attr == fill(5.0,n) )  # new_attr should be [5.0, 5.0, . . . ,5.0]
end
  
