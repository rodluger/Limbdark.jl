# Tests automatic differentiation on s2.jl:
# include("../src/transit_structure.jl")
# include("../src/s2.jl")
#= if VERSION >= v"0.7"
  using Test
else
  using Base.Test
end=#

import Limbdark: s2, s2!

# Randomizer seed
Random.seed!(42)

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
b=r=2.0
while b+r > 1
  global r = 2rand(); global b= 2rand()
end
runtest_s2(r,b,"r+b < 1")

# Now, try a random case with b+r > 1:
b=r=0.
while b+r < 1
  global r = 2rand(); global b= 2rand()
end
runtest_s2(r,b,"r+b > 1")

@testset "s2 grid" begin
# Now, try out hard cases - ingress/egress of small planets:
r0 = [0.01,0.01,1.0,10.0,100.0]
nb = 200
#l_max = 20
l_max = 10
n_max = l_max^2+2*l_max

#=
using PyPlot
fig,axes = subplots(1,length(r0))
get_cmap("plasma")
epsilon = 1e-12; delta = 1e-3
for i=1:length(r0)
  global r=r0[i]
  if r < 1.0
    bgrid = [linearspace(1e-15,epsilon,nb); linearspace(epsilon,delta,nb); linearspace(delta,r-delta,nb);
     -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
     linearspace(r+delta,1-r-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (1-r); linearspace(1-r-epsilon,1-r+epsilon,nb);
     logarithmspace(log10(epsilon),log10(delta),nb) .+ (1-r); linearspace(1-r+delta,1+r-delta,nb); - logarithmspace(log10(delta),log10(epsilon),nb) .+ (1+r);
     linearspace(1+r-epsilon,1+r-1e-15,nb)]
  else
    bgrid = [linearspace(r-1+1e-15,r-1+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ (r-1); linearspace(r-1+delta,r-delta,nb);
     -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
     linearspace(r+delta,r+1-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (r+1); linearspace(r+1-epsilon,r+1-1e-15,nb)]
  end
  nbgrid = length(bgrid)
  s2_jac_grid = zeros(nbgrid,2)
  s2_grid = zeros(nbgrid)
  s2_grid_big = zeros(nbgrid)
  s2_jac_grid_num = zeros(nbgrid,2)
  s2_grad_ana = zeros(2)
  for j=1:nbgrid
    s_2,Eofk,Em1mKdm= s2!(r,bgrid[j],s2_grad_ana)
    s2_grid[j]=s_2
    s2_jac_grid[j,:]=s2_grad_ana
    s2_big,s2_grad_numeric = s2_grad_num(r,bgrid[j])
    s2_grid_big[j]=s2_big
    s2_jac_grid_num[j,:]= s2_grad_numeric
    @test isapprox(s_2,s2_big,atol=1e-8)
    @test isapprox(s2_grad_ana[1],s2_grad_numeric[1],atol=1e-8)
    @test isapprox(s2_grad_ana[2],s2_grad_numeric[2],atol=1e-8)
  end
# Now, make plots:
  ax = axes[i]
  ax[:semilogy](bgrid,abs.(s2_grid[:,1]-s2_grid_big[:,1]) .+1e-18,lw=1,label="s 2 ")
  ax[:semilogy](bgrid,abs.(s2_jac_grid[:,1]-s2_jac_grid_num[:,1]) .+1e-18,lw=1,label="ds2/dr")
  ax[:semilogy](bgrid,abs.(s2_jac_grid[:,2]-s2_jac_grid_num[:,2]) .+1e-18,lw=1,label="ds2/db")
  ax[:legend](loc="upper right",fontsize=6)
  ax[:set_xlabel]("b values")
  ax[:set_ylabel]("Derivative Error")
end

savefig("test_s2_gradient.pdf", bbox_inches="tight")
=#
end
