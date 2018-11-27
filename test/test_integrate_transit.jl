# This is code for computing a transit model integrated over
# a time step.

include("../src/integrate_transit.jl")

# Test it out:

include("../test/loglinspace.jl")

t1 = -1.2; t2 = 1.2; nt = 1000; dt = (t2-t1)/nt*5
t = linearspace(t1,t2,nt)
favg0 = zeros(nt)
favg1 = zeros(nt)
depth1 = zeros(Int64,nt)
diff1 = zeros(nt)
favg2 = zeros(nt)
depth2 = zeros(Int64,nt)
diff2 = zeros(nt)
favg3 = zeros(nt)
depth3 = zeros(Int64,nt)
diff3 = zeros(nt)
favg4 = zeros(nt)
depth4 = zeros(Int64,nt)
diff4 = zeros(nt)
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

# Now, carry out speed test:

@time for i=1:nt
        favg1[i] = integrate_timestep(param,trans,t[i],dt,1e-5*r^2,8)
      end


