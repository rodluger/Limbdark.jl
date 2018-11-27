# This is test of the code for computing a transit model integrated over a time step.

include("../src/integrate_transit.jl")

# Test it out:

include("../test/loglinspace.jl")

t1 = -1.2; t2 = 1.2; nt = 1000; dt = (t2-t1)/nt*5
t = linearspace(t1,t2,nt)
favg0 = zeros(nt)
favg1 = zeros(nt); depth1 = zeros(Int64,nt); diff1 = zeros(nt)
favg2 = zeros(nt); depth2 = zeros(Int64,nt); diff2 = zeros(nt)
favg3 = zeros(nt); depth3 = zeros(Int64,nt); diff3 = zeros(nt)
favg4 = zeros(nt); depth4 = zeros(Int64,nt); diff4 = zeros(nt)
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.2,0.2,0.2,0.2,0.2]
trans = transit_init(r,b0,u_n,false)
param = [0.0,1.0,r,b0]
for i=1:nt
  trans.b = sqrt(param[4]^2+(param[2]*(t[i]-param[1]))^2)
  favg0[i] = transit_poly_d!(trans)
  favg1[i],depth1[i],diff1[i] = integrate_timestep_track(param,trans,t[i],dt,1e-5*r^2,4)
  favg2[i],depth2[i],diff2[i]  = integrate_timestep_track(param,trans,t[i],dt,1e-5*r^2,256)
  favg3[i],depth3[i],diff3[i]  = integrate_timestep_track(param,trans,t[i],dt,1e-10*r^2,4)
  favg4[i],depth4[i],diff4[i]  = integrate_timestep_track(param,trans,t[i],dt,1e-10*r^2,256)
end
favg1 /= dt
favg2 /= dt
favg3 /= dt
favg4 /= dt
# The depth arrays can be plotted to show how the maximum depth used varies across
# the light curve

# Now, carry out speed test:

@time for i=1:nt
        trans.b = sqrt(param[4]^2+(param[2]*(t[i]-param[1]))^2)
        favg0[i] = transit_poly_d!(trans)
      end
dtinv = inv(dt)
@time for i=1:nt
        favg1[i] = integrate_timestep(param,trans,t[i],dt,1e-5*r^2,8)*dtinv
      end

# Estimate the difference between time-integrated and point estimates of the flux
# by integrating the second order term in the Taylor-series expansion.

d2fdb2 = zeros(nt)
dfdb = zeros(nt)
trans_big = transit_init(big(r),big(b0),big.(u_n),true)
bigparam = [big(0.0),big(1.0),big(r),big(b0)]
dq = 1e-18
for i=1:nt
  trans_big.b = sqrt(bigparam[4]^2+(bigparam[2]*(t[i]-bigparam[1]))^2)
  trans_big.b -= dq
  transit_poly_d!(trans_big)
  dfdbm = trans_big.dfdrb[2]
  trans_big.b += 2dq
  transit_poly_d!(trans_big)
  dfdbp = trans_big.dfdrb[2]
  d2fdb2[i] = convert(Float64,(dfdbp-dfdbm)/(2*dq))
  dfdb[i] = convert(Float64,trans_big.dfdrb[2])
end

# Now plot the estimate of the difference:
plot(favg4-favg0)
b = sqrt.(param[4]^2+(param[2]*(t-param[1])).^2)
plot(dt^2*(d2fdb2.*(param[2]*(t-param[1])./b).^2+dfdb.*param[2]^2./b)/22.)
