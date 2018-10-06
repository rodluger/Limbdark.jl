# Benchmarks the version of transit_poly which uses a Type
# for holding all of the variables related to the transit
# computation.

# Defines the function that computes the transit:
include("../../../src/transit_poly_struct.jl")

nu = [1,2,3,5,8,13,21,34,55,89]; nnu = length(nu)
nb = [100,316,1000,3160,10000,31600,100000,316000,1000000]; nnb = length(nb)
timing_ratio = zeros(9,nnu)

# Define a function that loops over impact parameter
# and defines arrays to hold the results:

function profile_transit_poly(trans,nb)
fgrad = zeros(2+trans.n)
timing = zeros(length(nb))
tmean = zeros(9)
for j=1:length(nb)
  nb0 = nb[j]
  b = zeros(nb0)
  flux = zeros(nb0)
  for k=1:9
    tic()
    for i=1:nb0
      b[i] = sqrt(((i-float(nb0)/2)*2/float(nb0)*(1.0+2.0*trans.r))^2)
      trans.b = b[i]
      flux[i]= transit_poly!(trans)
    end
    tmean[k] = toq()
  end
  timing[j] = median(tmean)
#  println("nb = ",nb0," min(b): ",minimum(b)," max(b): ",maximum(b)," min(flux): ",minimum(flux)," max(flux): ",maximum(flux))
end
return timing
end

# Instantiate a transit with r=0.1, b=0.5, u_n=[0.2,0.3]:
r = 0.1; b = 0.5
# Load in PyPlot (Julia matplotlib wrapper):
using PyPlot
clf()
for iu = 1:nnu
  u_n = ones(iu)/iu
  trans = transit_init(r,b,u_n,true)
# Transform from u_n to c_n coefficients:
#compute_c_n_grad!(trans)  # Now included in transit_init function

# Call the function:
  @time timing=profile_transit_poly(trans,nb)
  if iu == 1
    timing_ratio[:,iu]=timing
  else
    timing_ratio[:,iu]=timing./timing_ratio[:,1]
  end 
  loglog(nb,timing)
  loglog(nb,timing,"o",label=string("n: ",nu[iu]))
  legend(loc="upper left")
end
xlabel("Number of points")
ylabel("Timing [sec]")
timing_ratio[:,1]=1.0
#clf()
#plot(b,flux-1,label="flux-1")
#plot(b,trans.r*flux_grad[:,1],label="r*dfdr")
#plot(b,trans.r*flux_grad[:,2],label="r*dfdb")
#plot(b,flux_grad[:,3],label="dfdu1")
#plot(b,flux_grad[:,4],label="dfdu2")
#legend(loc="lower right")
#println("b: ",minimum(b)," ",maximum(b))
#println("f: ",minimum(flux)," ",maximum(flux))

savefig("benchmark_poly_transit.pdf", bbox_inches="tight")

clf()
tmed = vec(median(timing_ratio,1))
semilogx(nu,tmed,"o")
semilogx(nu,tmed,label="Measured",linewidth=2)
alp = log(tmed[nnu])/log(nu[nnu])
semilogx(nu,nu.^alp,linestyle="--",label=L"$n^{0.237}$")
#plot(nu,(nu).^(1./3.),linestyle="--",label=L"$n^{1/3}$",linewidth=2)
#plot(nu,(nu).^(3./8.),linestyle="--",label=L"$n^{3/8}$",linewidth=2)
#plot(nu,(nu).^(3./8.),linestyle="--",label=L"$n^{3/8}$",linewidth=2)
xlabel("Number of limb-darkening coefficients")
ylabel("Relative timing")
legend(loc="upper left")
savefig("benchmark_limbdark_timing.pdf", bbox_inches="tight")
