# Tests analytic differentiation on transit_poly.jl:
#include("../src/transit_poly.jl")
# include("../src/transit_poly_struct.jl")
#include("../src/dJv_seriesdk.jl")
#=using PyPlot
if VERSION >= v"0.7"
  using Test
  if !(eval(:(@isdefined skip_plots)))
    skip_plots = false
  end
else
  using Base.Test
  if !(isdefined(:skip_plots))
    skip_plots = false
  end
end
=#
if !(eval(:(@isdefined skip_plots)))
  skip_plots = true
end

# Now, carry out finite-difference derivative computation:
function transit_poly_grad_num(r::T,b::T,u_n::Array{T,1},dq0::T) where {T <: Real}
  dq = big(dq0)
# Make BigFloat versions of r, b & u_n:
  r_big = big(r)
  b_big = big(b); u_big = big.(u_n)
# Compute flux to BigFloat precision:
  tp=transit_poly(r_big,b_big,u_big)
# Now, compute finite differences:
#  tp_grad_big= zeros(BigFloat,2+length(u_n))
  tp_grad_big= zeros(BigFloat,2)
  tp_plus = transit_poly(r_big+dq,b_big,u_big)
  tp_minus = transit_poly(r_big-dq,b_big,u_big)
  tp_grad_big[1] = (tp_plus-tp_minus)*.5/dq
  if r == 1.0 && length(u_n) > 2
    tp_grad_big[1] = (tp_plus-tp)/dq
  end
  tp_plus = transit_poly(r_big,b_big+dq,u_big)
  if b_big > dq
    tp_minus = transit_poly(r_big,b_big-dq,u_big)
    width = 2dq
  else
    tp_minus = tp
    width = dq
  end
  tp_grad_big[2] = (tp_plus-tp_minus)/width
#  for i=1:length(u_n)
#    u_tmp = copy(u_big); u_tmp[i] += dq
#    tp_plus = transit_poly(r_big,b_big,u_tmp)
#    u_tmp[i] -= 2dq
#    tp_minus = transit_poly(r_big,b_big,u_tmp)
#    tp_grad_big[2+i] = (tp_plus-tp_minus)*.5/dq
#  end
return convert(Float64,tp),convert(Array{Float64,1},tp_grad_big)
end

include("test_transit_poly_gradient_case.jl")

@testset "transit_poly_gradient" begin
# No limb-darkening:
#global u_n = [0.0]
#global test_name = "uniform"
#include("test_transit_poly_gradient_case.jl")
#fgrid0 = fgrid
fgrid0 = test_transit_poly_gradient_case([0.0],"uniform")
#read(STDIN,Char)
# Linear limb-darkening:
#global u_n = [1.0]
#global test_name = "linear"
#include("test_transit_poly_gradient_case.jl")
#fgrid1 = fgrid
fgrid1 = test_transit_poly_gradient_case([1.0],"linear")
#read(STDIN,Char)
# Quadratic limb-darkening:
#global u_n = [2.0,-1.0]
#global test_name = "quadratic"
#include("test_transit_poly_gradient_case.jl")
#fgrid2 = fgrid
fgrid2 = test_transit_poly_gradient_case([2.0,-1.0],"quadratic")
#read(STDIN,Char)
# Cubic limb-darkening:
#global u_n = [3.0,-3.0,1.0]
#global test_name = "cubic"
#include("test_transit_poly_gradient_case.jl")
#fgrid3 = fgrid
fgrid3 = test_transit_poly_gradient_case([3.0,-3.0,1.0],"cubic")
fgrid4 = test_transit_poly_gradient_case(ones(4)*0.1,"4th order")
fgrid5 = test_transit_poly_gradient_case(ones(5)*0.1,"5th order")
# Now, for arbitrary limb-darkening:
#global u_n = ones(10)*0.1
#global test_name = "10th order"
#include("test_transit_poly_gradient_case.jl")
#fgridn = fgrid
fgridn = test_transit_poly_gradient_case(ones(10)*0.1,"10th order")

if ~skip_plots
  r0 = [0.01,0.1,0.5,0.999999,1.0,1.000001,2.0,10.,100.0]
  r0_name =["0.01","0.1","0.5","0.99","1","1.01","2","10","100"]
  nb = 50

  epsilon = 1e-9; delta = 1e-3
  fig,axes = subplots(3,3)
  for i=1:length(r0)
    global r=r0[i]
    ax = axes[i]
    if r < 1.0
      bgrid = abs.([linearspace(0.0,epsilon,nb); linearspace(epsilon,delta,nb); linearspace(delta,r-delta,nb);
       -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
       linearspace(r+delta,1-r-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (1-r); linearspace(1-r-epsilon,1-r+epsilon,nb);
       logarithmspace(log10(epsilon),log10(delta),nb) .+ (1-r); linearspace(1-r+delta,1+r-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (1+r);linearspace(1+r-epsilon,1+r-1e-13,nb)])
       nticks = 14
       xticknames=[L"$10^{-13}$",L"$10^{-12}$",L"$10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
       L"$1-r-10^{-3}$",L"$1-r-10^{-12}$",L"$1-r+10^{-12}$",L"$1-r+10^{-3}$",L"$1+r-10^{-3}$",L"$1+r-10^{-12}$",L"$1+r-10^{-13}$"]
    else
      bgrid = [r-1+1e-13; logarithmspace(log10(epsilon),log10(delta),nb) .+ (r-1); linearspace(r-1+delta,r-delta,nb);
       -logarithmspace(log10(delta),log10(epsilon),nb) .+ r; linearspace(r-epsilon,r+epsilon,nb); logarithmspace(log10(epsilon),log10(delta),nb) .+ r;
       linearspace(r+delta,r+1-delta,nb); -logarithmspace(log10(delta),log10(epsilon),nb) .+ (r+1); r+1-1e-13]
       nticks = 8
       xticknames=[L"$r-1+10^{-13}$",L"$r-1+10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
       L"$r+1-10^{-3}$",L"$r+1-10^{-13}$"]
    end
    nbgrid = length(bgrid)
    ax[:plot](bgrid,fgrid0[i,1:nbgrid],label="Uniform")
    ax[:plot](bgrid,fgrid1[i,1:nbgrid],label="Linear")
    ax[:plot](bgrid,fgrid2[i,1:nbgrid],label="Quadratic")
    ax[:plot](bgrid,fgrid3[i,1:nbgrid],label="Cubic")
    ax[:plot](bgrid,fgrid4[i,1:nbgrid],label="N=4")
    ax[:plot](bgrid,fgrid5[i,1:nbgrid],label="N=5")
    ax[:plot](bgrid,fgridn[i,1:nbgrid],label="N=10")
    if mod(i,3) == 0
      ax[:set_xlabel]("b")
    end
    if i <= 3
      ax[:set_ylabel]("F")
    end
    ax[:legend](loc="upper left",fontsize=6)
    ax[:set_title](string("r = ",r0[i]),fontsize=6)
  end
end
end
