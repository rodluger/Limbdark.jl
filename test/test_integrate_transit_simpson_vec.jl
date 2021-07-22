
# This is test of the code for computing a transit model integrated over a time step.
# Computes derivatives over the timestep.
# Test it out:

import Limbdark: Transit_Struct

#@testset "test_integrate" begin

# Create a time array:
t1 = -1.5; t2 = 1.5; nt = 10001; dt = (t2-t1)/10
#t1 = -1.5; t2 = 1.5; nt = 40001; dt = (t2-t1)/10
t = zeros(nt)
t .= linearspace(t1,t2,nt)
# The following compares different tolerances and maxdepths:
r = 0.1; b0 = 0.5; u_n = [0.3,0.3]; nu = length(u_n)
favg0 = zeros(nt,5+nu)
neval0= ones(Int64,nt)
favg1 = zeros(nt,5+nu)
favg2 = zeros(nt,5+nu)
neval1= zeros(Int64,nt)
neval2= zeros(Int64,nt)
depthmax1= zeros(Int64,nt)
depthmax2= zeros(Int64,nt)
# Initialize transit structure, trans:
trans = transit_init(r,b0,u_n,true)
# Set parameters for the transit light curve:
param = [0.0,1.0,b0]   # [t_0,v,b_0]

# Now, carry out tests.  First, compute the light curve over the finely spaced grid:
function compute_lightcurve!(trans::Transit_Struct{T},param::Array{T,1},t::Array{T,1},favg0::Array{T,2},nt::Int64,neval::Array{Int64,1}) where {T <: Real}
@inbounds for i=1:nt
  if neval[i] == 0
    neval[i] = 1
  end
  for j=1:neval[i]
    trans.b = sqrt(param[3]^2+(param[2]*(t[i]-param[1]))^2)
    favg0[i,1] =  transit_poly!(trans)
    binv = inv(trans.b)
    favg0[i,2] =  trans.dfdrb[1]
    favg0[i,3] =  trans.dfdrb[2]*binv*param[2]^2*(param[1]-t[i])
    favg0[i,4] =  trans.dfdrb[2]*param[2]*binv*(t[i]-param[1])^2
    favg0[i,5] =  trans.dfdrb[2]*param[3]*binv
    favg0[i,6:5+trans.n] = trans.dfdu
  end
end
return
end

compute_lightcurve!(trans,param,t,favg0,nt,neval0)
@time compute_lightcurve!(trans,param,t,favg0,nt,neval0)
# Subtract 1 from the flux:
favg0[:,1] .-= 1

# Now, integrate this lightcurve by simply averaging over this grid using trapezoidal rule:
fsmooth = copy(favg0)
# Each time step is 3e-4, which happens to be 1/1000 of the integration time:
nint = 500
# Each time step is 0.75e-4, which happens to be 1/4000 of the integration time:
#nint = 2000
for i=1:nt
  fsum = zeros(5+nu)
  # Starting edge has weight of 1/2:
  if i-nint < 1
#    fsum[1] += 0.5
  else
    fsum += favg0[i-nint,:]*0.5
  end
  for j=i-nint+1:i+nint-1
    if j < 1 || j > nt
      # Outside of range set the flux to one:
#      fsum[1] += 1.0
      # Outside of range set the flux-1 to zero:
      # Outside of range set the derivatives to zero:
    else
      fsum += favg0[j,:]
    end
  end
  # Ending edge has weight of 0.5:
  if i+nint > nt
#    fsum[1] += 0.5
  else
    fsum += favg0[i+nint,:]*0.5
  end
  # Divide by the number of substeps:
  fsum /= nint*2
  fsmooth[i,:] = fsum
end

integrate_lightcurve!(trans,param,t,dt,favg1,nt,1e-11,40,neval1,depthmax1)
integrate_lightcurve!(trans,param,t,dt,favg2,nt,1e-13,60,neval2,depthmax2)

# Plot the adaptive-integrated result:
for j=1:5+nu
clf()
plot(t,favg1[:,j])
# Plot the unbinned results:
plot(t,favg0[:,j])
# Plot the grid-integrated result:
plot(t,fsmooth[:,j],linestyle=":")
#if VERSION >= v"0.7"
#  read(stdin,Char)
#else
#  read(STDIN,Char)
#end
end

# Now, plot differences between smoothed versions:
clf()
labels = ["flux","r","t0","v","b0","u1","u2","u3"]
for j=1:5+nu
#for j=1:1
#  plot(t,(favg1[j,:]-fsmooth[j,:])/maximum(abs.(favg1[j,:])),label=labels[j])
#  plot(t,(favg1[j,:]-favg2[j,:])/maximum(abs.(favg2[j,:])),label=labels[j],linestyle=":")
  plot(t,(favg1[:,j]-favg2[:,j]),label=labels[j],linestyle=":")
#  if j != 1
    @test (maximum(abs.(favg1[:,j]-fsmooth[:,j])) < maximum(abs.(favg1[:,j]))/nint^2*10.0)
    @test (maximum(abs.(favg1[:,j]-favg2[:,j])) < maximum(abs.(favg2[:,j]))*1e-7)
    println(j," ",labels[j]," ",maximum(abs.(favg1[:,j]-favg2[:,j]))," ",maximum(abs.(favg2[:,j]))," ",maximum(abs.(favg1[:,j]-favg2[:,j]))/maximum(abs.(favg2[:,j])))
#  else
#    # For the flux, subtract 1 when taking fractional difference:
#    @test (maximum(abs.(favg1[j,:]-fsmooth[j,:])) < maximum(abs.(favg1[j,:]-1.0))/nint^2*10.0)
#    @test (maximum(abs.(favg1[j,:]-favg2[j,:])) < maximum(abs.(favg2[j,:]))*1e-7)
#    println(j," ",labels[j]," ",maximum(abs.(favg1[j,:]-favg2[j,:]))," ",maximum(abs.(favg2[j,:]-1.0))," ",maximum(abs.(favg1[j,:]-favg2[j,:]))/maximum(abs.(favg2[j,:]-1.0)))
#  end
#  println(j," ",labels[j]," ",maximum(abs.(favg1[j,:]-favg2[j,:]))/maximum(abs.(favg2[j,:])))
end
legend(loc = "lower left")
#end
