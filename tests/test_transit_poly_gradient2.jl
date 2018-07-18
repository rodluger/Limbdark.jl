# Tests analytic differentiation on transit_poly.jl:
#include("../src/transit_poly.jl")
include("../src/transit_poly_struct.jl")
#include("../src/dJv_seriesdk.jl")
using PyPlot


# Now, carry out finite-difference derivative computation:
function transit_poly_grad_num(r::T,b::T,u_n::Array{T,1},dq0::T) where {T <: Real}
  dq = big(dq0)
# Make BigFloat versions of r, b & u_n:
  r_big = big(r); b_big = big(b); u_big = big.(u_n)
# Compute flux to BigFloat precision:
  tp=transit_poly(r_big,b_big,u_big)
# Now, compute finite differences:
  tp_grad_big= zeros(BigFloat,2+length(u_n))
  tp_plus = transit_poly(r_big+dq,b_big,u_big)
  tp_minus = transit_poly(r_big-dq,b_big,u_big)
  tp_grad_big[1] = (tp_plus-tp_minus)*.5/dq
  tp_plus = transit_poly(r_big,b_big+dq,u_big)
  if b_big > dq
    tp_minus = transit_poly(r_big,b_big-dq,u_big)
    width = 2dq
  else
    tp_minus = tp
    width = dq
  end
  tp_grad_big[2] = (tp_plus-tp_minus)/width
  for i=1:length(u_n)
    u_tmp = copy(u_big); u_tmp[i] += dq
    tp_plus = transit_poly(r_big,b_big,u_tmp)
    u_tmp[i] -= 2dq
    tp_minus = transit_poly(r_big,b_big,u_tmp)
    tp_grad_big[2+i] = (tp_plus-tp_minus)*.5/dq
  end
return convert(Float64,tp),convert(Array{Float64,1},tp_grad_big)
end

# No limb-darkening:
u_n = [0.0]
test_name = "uniform"
include("test_transit_poly_gradient_case.jl")
fgrid0 = fgrid
#read(STDIN,Char)
# Linear limb-darkening:
u_n = [1.0]
test_name = "linear"
include("test_transit_poly_gradient_case.jl")
fgrid1 = fgrid
#read(STDIN,Char)
# Quadratic limb-darkening:
u_n = [2.0,-1.0]
test_name = "quadratic"
include("test_transit_poly_gradient_case.jl")
fgrid2 = fgrid
#read(STDIN,Char)
# Cubic limb-darkening:
u_n = [3.0,-3.0,1.0]
test_name = "cubic"
include("test_transit_poly_gradient_case.jl")
fgrid3 = fgrid
# Now, for arbitrary limb-darkening:
u_n = ones(10)*0.1
test_name = "10th order"
include("test_transit_poly_gradient_case.jl")
fgridn = fgrid

if ~skip_plots
  r0 = [0.01,.1,0.5,0.99,1.0,1.01,2.0,10.,100.0]; n_u = length(u_n)
  r0_name =["0.01","0.1","0.5","0.99","1","1.01","2","10","100"]
  nb = 50

  epsilon = 1e-12; delta = 1e-3
  fig,axes = subplots(3,3)
  for i=1:length(r0)
    r=r0[i]
    ax = axes[i]
    if r < 1.0
#    b = [linspace(1e-15,epsilon,nb); linspace(epsilon,delta,nb); linspace(delta,r-delta,nb);
#    b = [linspace(1e-13,epsilon,nb); linspace(epsilon,delta,nb); linspace(delta,r-delta,nb);
      b = [linspace(0.0,epsilon,nb); linspace(epsilon,delta,nb); linspace(delta,r-delta,nb);
       r-logspace(log10(delta),log10(epsilon),nb); linspace(r-epsilon,r+epsilon,nb); r+logspace(log10(epsilon),log10(delta),nb);
       linspace(r+delta,1-r-delta,nb); 1-r-logspace(log10(delta),log10(epsilon),nb); linspace(1-r-epsilon,1-r+epsilon,nb);
#     1-r+logspace(log10(epsilon),log10(delta),nb); linspace(1-r+delta,1+r-delta,nb); 1+r-logspace(log10(delta),log10(epsilon),nb);linspace(1+r-epsilon,1+r,nb)]
       1-r+logspace(log10(epsilon),log10(delta),nb); linspace(1-r+delta,1+r-delta,nb); 1+r-logspace(log10(delta),log10(epsilon),nb);linspace(1+r-epsilon,1+r-1e-13,nb)]
       b = abs.(b)
       nticks = 14
       xticknames=[L"$10^{-13}$",L"$10^{-12}$",L"$10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
       L"$1-r-10^{-3}$",L"$1-r-10^{-12}$",L"$1-r+10^{-12}$",L"$1-r+10^{-3}$",L"$1+r-10^{-3}$",L"$1+r-10^{-12}$",L"$1+r-10^{-13}$"]
    else
      b = [r-1+1e-13;r-1+logspace(log10(epsilon),log10(delta),nb); linspace(r-1+delta,r-delta,nb);
       r-logspace(log10(delta),log10(epsilon),nb); linspace(r-epsilon,r+epsilon,nb); r+logspace(log10(epsilon),log10(delta),nb);
       linspace(r+delta,r+1-delta,nb); r+1-logspace(log10(delta),log10(epsilon),nb);r+1-1e-13]
       nticks = 8
       xticknames=[L"$r-1+10^{-13}$",L"$r-1+10^{-3}$",L"$r-10^{-3}$",L"$r-10^{-12}$",L"$r+10^{-12}$",L"$r+10^{-3}$",
       L"$r+1-10^{-3}$",L"$r+1-10^{-13}$"]
    end
    ax[:plot](b,fgrid0[i,1:length(b)],label="Uniform")
    ax[:plot](b,fgrid1[i,1:length(b)],label="Linear")
    ax[:plot](b,fgrid2[i,1:length(b)],label="Quadratic")
    ax[:plot](b,fgrid3[i,1:length(b)],label="Cubic")
    ax[:plot](b,fgridn[i,1:length(b)],label="N=10")
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
#return
#end
