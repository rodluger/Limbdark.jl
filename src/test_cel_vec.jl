# Tests the vectorized version of cel:
#include("cel_bulirsch.jl")
include("/Users/ericagol/Computer/Julia/Transit/cel_bulirsch.jl")

p = -rand(); a=rand(3); b=rand(3)
kc = rand(); k2 = 1.-kc^2
ell1 = zeros(3); ell2=zeros(3)
p0 = sqrt((kc*kc-p)/(1-p))
a0 = (a[1]-b[1])/(1-p)
b0 = -(b[1]-a[1]*p)/(1-p)^2*(1-kc^2)/p0+a0*p0

ell1[1] = cel_bulirsch(k2,kc,p,a[1],b[1])
ell1_alt = cel_bulirsch(k2,kc,p0^2,a0,b0*p0)
ell1[2] = cel_bulirsch(k2,kc,1.0,a[2],b[2])
ell1[3] = cel_bulirsch(k2,kc,1.0,a[3],b[3])

nphi = 10000; dphi = .5*pi/nphi
phi0 = asin(sqrt(1./(1.-p)))
phi = linspace(.5*dphi,.5*pi-dphi,nphi)
ell3 = zeros(3); cphi2 = cos.(phi).^2; sphi2=sin.(phi).^2
den = dphi./sqrt.(cphi2+kc*kc*sphi2)
ell3[1] = sum((a[1]*cphi2+b[1]*sphi2)./(cphi2+p*sphi2).*den)
ell3[2] = sum((a[2]*cphi2+b[2]*sphi2).*den)
ell3[3] = sum((a[3]*cphi2+b[3]*sphi2).*den)
#ell2 = cel_bulirsch(k2,kc,p,a,b)
println("cel_bul: ",ell1)
println("cel_alt: ",ell1_alt)
println("numeric: ",ell3)
#println("sca-vec: ",ell1-ell2)
println("sca-num: ",ell1-ell3)
