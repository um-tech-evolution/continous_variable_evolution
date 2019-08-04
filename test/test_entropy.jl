# test ../src/cont_entropy.jl

C = [0.1,0.15,0.3]
w = 0.1
a = 0.0
b = 0.5

sb = steps_bounds(C,w,a,b)
@assert isapprox(1.0,step_bounds_check(sb))
bl = bounds_list(sb,a,b)
vl = values_list(bl)
cont_ent = values_list(bounds_list(steps_bounds(C,w,a,b),a,b))
@assert isapprox(cont_ent, 1.9182958340544896)
c_ent = cont_entropy(C,w,a,b)
@assert isapprox(c_ent, 1.9182958340544896)

C2 = [0.02,0.04,0.06,0.44,0.47,0.49]
sb2 = steps_bounds(C2,w,a,b)
@assert isapprox(1.0,step_bounds_check(sb2))
bl2 = bounds_list(sb2,a,b)
vl2 = values_list(bl2)
@assert isapprox(vl2, 2.3921335657167524)
c_ent = cont_entropy(C2,w,a,b)
@assert isapprox(c_ent, 2.3921335657167524)
