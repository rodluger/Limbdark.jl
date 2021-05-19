import Limbdark: Mn_series_coeff!

# Randomizer seed
Random.seed!(42)

#@testset "transit_poly" begin

# Now, integrate by hand:
function transit_poly_int(r,b,u,ns)
n = length(u)
if b < (1.0+r)
  s_1 = maximum([0.0,b-r]); s_2=minimum([1.0,b+r])
  ds=(s_2-s_1)/ns
  s = linearspace(s_1+.5*ds,s_2-.5*ds,ns)
  fobs = zero(r)
  for j=1:ns
    if s[j] < r-b
      dphi = 2pi
    else
      dphi = 2*acos((s[j]^2+b^2-r^2)/(2*b*s[j]))
    end
    z = sqrt(1-s[j]^2)
    imu = 1.0
    for i=1:n
      imu -= u[i]*(1-z)^i
    end
    fobs += s[j]*ds*dphi*imu
  end
  norm = 1.
  for i=1:n
    norm -= 2*u[i]/(2+3*i+i^2)
  end
  fobs /= pi*norm
  fobs = 1-fobs
  return fobs
else
  return 1.0
end
end

r0=[0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,1.0,1.0,1.0,10.,10.,10.,100.,100.,100.]
b0=[0.,0.01,0.99,1.0,1.01, 0.,0.1,0.9,1.0,1.1, 0.,1.0,2.0, 9.,10.,11., 99.,100.,101.]

n = 2+ceil(Int64,rand()*20); u = rand(n); u *= rand()/sum(u)

ns = 5000
trans = transit_init(r0[1],b0[1],u,false)
for i=1:length(r0)
  trans.r = r0[i]; trans.b = b0[i]
  flux = transit_poly!(trans)
  @time flux = transit_poly!(trans)
  f_num = transit_poly_int(r0[i],b0[i],u,ns)
  @time f_num = transit_poly_int(r0[i],b0[i],u,ns)
  println("r: ",r0[i]," b: ",b0[i]," f_an: ",flux," f_num: ",f_num," diff: ",flux-f_num)
  @test isapprox(f_num,flux,atol=1e-5)
end

# Now, compute light curves and plot:

npts = 1000
iu = 30
u_n = ones(iu) / iu
b0 = zeros(npts)
for i=1:npts
    b0[i] = 3.0 * (i / float(npts) - 0.5)
end

lc_ana = zeros(npts); lc_num=zeros(npts)
rtest = 0.1
trans = transit_init(rtest,b0[1],u_n,false)

for i=1:length(b0)
#  if mod(i,10) == 0
#    println("i: ",i," b: ",b)
#  end
#  lc_ana[i] = transit_poly(r,abs(b),u_n)
  trans.b = abs(b0[i])
  lc_ana[i] = transit_poly!(trans)
  lc_num[i] = transit_poly_int(rtest,abs(b0[i]),u_n,ns)
end

using PyPlot

plot(b0,lc_ana, linewidth=2, label="Analytic")
plot(b0,lc_num, linewidth=1, label="Numeric")

println("Maximum difference of lightcurve for N=",iu,": ",maximum(abs,lc_ana-lc_num))
@test isapprox(lc_ana,lc_num,atol=1e-5)

#end
