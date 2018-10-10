# Tests automatic differentiation on s2.jl:
include("../src/s2.jl")

using Test
using ForwardDiff
using DiffResults


function s2_grad_func(r::T,b::T) where {T <: Real}
  # Computes the derivative of s_n(r,b) with respect to r, b.
  # Create a vector for use with ForwardDiff
  x=[r,b]
  # Now, define a wrapper of s2 for use with ForwardDiff:
  function diff_s2(x::Array{T,1}) where {T <: Real}
  # x should be a two-element vector with values [r,b]
  rc,bc = x
  s2c= s2(rc,bc)
  return s2c
  end

  # Set up a type to store s_n and it's Jacobian with respect to x:
  out = DiffResults.GradientResult(x)
  # Compute the Jacobian (and value):
  out = ForwardDiff.gradient!(out,diff_s2,x)
  # Place the value in the s_2 vector:
  s_2 = DiffResults.value(out)
  # And, place the Jacobian in an array:
  s2_gradient= DiffResults.gradient(out)
return s_2,s2_gradient
end

#r = 0.1; b= 0.95
#r = 0.1; b= 1.0-r
#r = 0.1; b= r
#r = 100.0; b=100.5

function s2_grad_num(r::T,b::T) where {T <: Real}
# Now, carry out finite-difference:
dq = big(1e-18)

# Allocate an array for s_2:
# Make BigFloat versions of r & b:
r_big = big(r); b_big = big(b)
# Compute s_n to BigFloat precision:
s_2_big= s2(r_big,b_big)
# Now, compute finite differences:
s2_grad_big= zeros(BigFloat,2)
s2_plus= s2(r_big+dq,b_big)
s2_minus= s2(r_big-dq,b_big)
s2_grad_big[1] = (s2_plus-s2_minus)*.5/dq
s2_plus=s2(r_big,b_big+dq)
s2_minus=s2(r_big,b_big-dq)
s2_grad_big[2] = (s2_plus-s2_minus)*.5/dq
return convert(Float64,s_2_big),convert(Array{Float64,1},s2_grad_big)
end

# In the following lines I test the derivatives for the few special
# cases needed for computing s_2.  Note that in the paper I report
# \Lambda, while this gets multiplied by -\pi to produce s_2.

function runtest_s2(r,b,test_name)
@testset "$test_name" begin
s_2,s2_gradient= s2_grad_func(r,b)
s2_big,s2_grad_numeric = s2_grad_num(r,b)
s2_grad_ana = zeros(2)
diff1 = s2_grad_numeric-s2_gradient
s_2,Eofk,Em1mKdm=s2!(r,b,s2_grad_ana)
diff2 = s2_grad_ana-s2_gradient
println("b : ",b," r: ",r," flux diff: ",s_2-s2_big," grad diff(num-auto): ",diff1," diff(ana-auto): ",diff2 )
#println("gradient: ",s2_gradient," expected: ",-[2,-2./3.])
@test isapprox(s2_grad_numeric[1],s2_gradient[1],atol=1e-8)
@test isapprox(s2_grad_numeric[2],s2_gradient[2],atol=1e-8)
@test isapprox(s2_grad_ana[1],s2_gradient[1],atol=1e-8)
@test isapprox(s2_grad_ana[2],s2_gradient[2],atol=1e-8)
end
return
end

# Try r=b=1/2 special case:
runtest_s2(0.5,0.5,"r=0.5;b=0.5")

# Try b=0 special case:
runtest_s2(0.3,0.0,"r=0.3;b=0.0")

# Try r+b=1 special case:
runtest_s2(0.2,0.8,"r=0.2;b=0.8")

# Try r=b <1/2 special case:
runtest_s2(0.3,0.3,"r=0.3;b=0.3")

# Try r=b > 1/2 special case:
runtest_s2(3.0,3.0,"r=3.0;b=3.0")

# Now, try a random case with b+r < 1:
runtest_s2(0.435,0.232,"r+b < 1")

# Now, try a random case with b+r > 1:
runtest_s2(11.434,11.113,"r+b > 1")

@testset "s2 grid" begin
# Now, try out hard cases - ingress/egress of small planets:
r0 = [0.01,0.01,1.0,10.0,100.0]
nb = 200
#l_max = 20
l_max = 10
n_max = l_max^2+2*l_max

using PyPlot
fig,axes = subplots(1,length(r0))
get_cmap("plasma")
epsilon = 1e-12; delta = 1e-3
i=1
for i=1:length(r0)
  r=r0[i]
  if r < 1.0
    b = [range(1e-15,stop=epsilon,length=nb); range(epsilon,stop=delta,length=nb); range(delta,stop=r-delta,length=nb);
     r .- 10 .^ range(log10(delta),stop=log10(epsilon),length=nb); range(r-epsilon,stop=r+epsilon,length=nb); r .+ 10 .^ range(log10(epsilon),stop=log10(delta),length=nb);
     range(r+delta,stop=1-r-delta,length=nb); (1 - r) .- 10 .^ range(log10(delta),stop=log10(epsilon),length=nb); range(1-r-epsilon,stop=1-r+epsilon,length=nb);
     (1 - r) .+ 10 .^ range(log10(epsilon),stop=log10(delta),length=nb); range(1-r+delta,stop=1+r-delta,length=nb); (1 + r) .- 10 .^ range(log10(delta),stop=log10(epsilon),length=nb);range(1+r-epsilon,stop=1+r-1e-15,length=nb)]
  else
    b = [range(r-1+1e-15,stop=r-1+epsilon,length=nb); (r - 1) .+ 10 .^ range(log10(epsilon),stop=log10(delta),length=nb); range(r-1+delta,stop=r-delta,length=nb);
     r .- 10 .^ range(log10(delta),stop=log10(epsilon),length=nb); range(r-epsilon,stop=r+epsilon,length=nb); r .+ 10 .^ range(log10(epsilon),stop=log10(delta),length=nb);
     range(r+delta,stop=r+1-delta,length=nb); (r + 1) .- 10 .^ range(log10(delta),stop=log10(epsilon),length=nb); range(r+1-epsilon,stop=r+1-1e-15,length=nb)]
  end
  igrid=range(1,stop=length(b),length=length(b)) .- 1
  s2_jac_grid = zeros(length(b),2)
  s2_grid = zeros(length(b))
  s2_grid_big = zeros(length(b))
  s2_jac_grid_num = zeros(length(b),2)
  s2_grad_ana = zeros(2)
  for j=1:length(b)
#    println("r: ",r," b: ",b[j])
#    s_2,s2_gradient= s2_grad(r,b[j])
    s_2,Eofk,Em1mKdm= s2!(r,b[j],s2_grad_ana)
    s2_grid[j]=s_2
    s2_jac_grid[j,:]=s2_grad_ana
    s2_big,s2_grad_numeric = s2_grad_num(r,b[j])
    s2_grid_big[j]=s2_big
    s2_jac_grid_num[j,:]= s2_grad_numeric
#    println("r: ",r," b: ",b[j]," s2: ",s_2," s2_big: ",s2_big)
    @test isapprox(s_2,s2_big,atol=1e-8)
    @test isapprox(s2_grad_ana[1],s2_grad_numeric[1],atol=1e-8)
    @test isapprox(s2_grad_ana[2],s2_grad_numeric[2],atol=1e-8)
#  println("r: ",r," b: ",b[j]," ds2/dr: ",s2_jac_grid[j,1]," ",s2_jac_grid[j,1]-s2_jac_grid_num[j,1])
#  println("r: ",r," b: ",b[j]," ds2/db: ",s2_jac_grid[j,2]," ",s2_jac_grid[j,2]-s2_jac_grid_num[j,2])
  end
# Now, make plots:
  ax = axes[i]
  ax[:semilogy](b,abs.(s2_grid[:,1]-s2_grid_big[:,1]) .+ 1e-18,lw=1,label="s 2 ")
  ax[:semilogy](b,abs.(s2_jac_grid[:,1]-s2_jac_grid_num[:,1]) .+ 1e-18,lw=1,label="ds2/dr")
  ax[:semilogy](b,abs.(s2_jac_grid[:,2]-s2_jac_grid_num[:,2]) .+ 1e-18,lw=1,label="ds2/db")
#  ax[:semilogy](b,abs.(asinh.(s2_jac_grid[:,1])-asinh.(s2_jac_grid_num[:,1])),lw=1,label="ds2/dr")
#  ax[:semilogy](b,abs.(asinh.(s2_jac_grid[:,2])-asinh.(s2_jac_grid_num[:,2])),lw=1,label="ds2/db")
  ax[:legend](loc="upper right",fontsize=6)
  ax[:set_xlabel]("b values")
  ax[:set_ylabel]("Derivative Error")
#  ax[:axis]([0,length(b),1e-16,1])
end

savefig("test_s2_gradient.pdf", bbox_inches="tight")

end
