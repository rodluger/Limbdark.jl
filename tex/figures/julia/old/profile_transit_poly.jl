# Profiles transit_poly.jl

r = 0.1  # radius ratio is 0.1
u = [0.2,0.3]; n = length(u)
N_c = n+2

const c_n = zeros(typeof(r),n+3)
const dfdrbc = zeros(typeof(r),n+3)
const a_n = zeros(typeof(r),n+1)
const dadu = zeros(typeof(r),n+1,n)
const dcdu = zeros(typeof(r),n+3,n)
const sn = zeros(typeof(r),N_c+1)
const dsndr = zeros(typeof(r),N_c+1)
const dsndb = zeros(typeof(r),N_c+1)
if iseven(N_c)
  v_max = round(Int64,N_c/2)+2
else
  v_max = round(Int64,(N_c-1)/2)+2
end
const Iv = zeros(typeof(r),v_max+1)
const s2_grad = zeros(typeof(r),2)
const Jv = zeros(typeof(r),v_max+1)
const dIvdk = zeros(typeof(r),v_max+1)
const dJvdk = zeros(typeof(r),v_max+1)

#include("../../../src/transit_poly_prealloc.jl")
include("../../../src/transit_poly_struct.jl")
#include("../../../src/dJv_seriesdk.jl")
using PyPlot

# Now, run the same with transit_poly:

function profile_transit_poly()
nb = 1000001; b=zeros(nb)
flux = zeros(nb); fgrad = zeros(4)
flux_big = zeros(nb); fgrad_big = zeros(BigFloat,4)
flux_grad = zeros(nb,4)
flux_grad_big = zeros(nb,4)
tic()
for i=0:1000000
  b[i+1] = sqrt(((i-500000.0)/500000.0*(1.0+2.0*trans.r))^2)
  flux[i+1]= transit_poly!(r,b[i+1],u,fgrad)
  flux_grad[i+1,:] = fgrad
end
toc()
plot(b,flux)
println("b: ",minimum(b)," ",maximum(b))
println("f: ",minimum(flux)," ",maximum(flux))
return
end
