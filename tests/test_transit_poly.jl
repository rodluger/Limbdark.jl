include("../src/transit_poly.jl")
#nu = 2+ceil(Int64,rand()*20); r=rand(); b=rand()*(1+r); u = rand(nu); u *= rand()/sum(u)
#r=rand(); b=rand()*(1+r); u = [0.,0.,0.,0.,1.0]; n=length(u)

@testset "transit_poly" begin

# Now, integrate by hand:
function transit_poly_int(r,b,u)
s_1 = maximum([0.0,b-r]); s_2=minimum([1.0,b+r])
ns = 10000; ds=(s_2-s_1)/ns
s = linspace(s_1+.5*ds,s_2-.5*ds,ns)
fobs = zero(r)
for j=1:ns
  if s[j] < r-b
    dphi = 2pi
  else
    dphi = 2*acos((s[j]^2+b^2-r^2)/(2*b*s[j]))
  end
  z = sqrt(1-s[j]^2)
  imu = 1.0
  for n=1:n
    imu -= u[n]*(1-z)^n
  end
  fobs += s[j]*ds*dphi*imu
end
norm = 1.
for n=1:n
  norm -= 2*u[n]/(2+3*n+n^2)
end
fobs /= pi*norm
fobs = 1-fobs
return fobs
end

r0=[0.01,0.01,0.01,0.01,0.01,0.1,0.1,0.1,0.1,0.1,1.0,1.0,1.0,10.,10.,10.,100.,100.,100.]
b0=[0.,0.01,0.99,1.0,1.01, 0.,0.1,0.9,1.0,1.1, 0.,1.0,2.0, 9.,10.,11., 99.,100.,101.]

r = r0[1]
n = 2+ceil(Int64,rand()*20); u = rand(n); u *= rand()/sum(u)
#N_c = n+2
#const c_n = zeros(typeof(r),n+3)
#const dfdrbc = zeros(typeof(r),n+3)
#const a_n = zeros(typeof(r),n+1)
#const dadu = zeros(typeof(r),n+1,n)
#const dcdu = zeros(typeof(r),n+3,n)
#const sn = zeros(typeof(r),N_c+1)
#const dsndr = zeros(typeof(r),N_c+1)
#const dsndb = zeros(typeof(r),N_c+1)
#if iseven(N_c)
#  v_max = round(Int64,N_c/2)+2
#else
#  v_max = round(Int64,(N_c-1)/2)+2
#end
#const Iv = zeros(typeof(r),v_max+1)
#const Jv = zeros(typeof(r),v_max+1)
#const dIvdk = zeros(typeof(r),v_max+1)
#const dJvdk = zeros(typeof(r),v_max+1)

for i=1:length(r0)
  r=r0[i]; b=b0[i]
  flux = transit_poly(r,b,u)
  @time flux = transit_poly(r,b,u)
  f_num = transit_poly_int(r,b,u)
  @time f_num = transit_poly_int(r,b,u)
  @test isapprox(f_num,flux,atol=1e-5)
  println("r: ",r," b: ",b," f_an: ",flux," f_num: ",f_num," diff: ",flux-f_num)
end

end
