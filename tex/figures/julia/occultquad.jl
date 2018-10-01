function occultquad(z0::Float64,u1::Float64,u2::Float64,p::Float64)
#
# Julia implementation of Mandel & Agol (2002).
# Computes the limb-darkened transit (or eclipse) of a star.
# Assumes quadratic limb-darkening:  I(\mu) \propto 1-u1*(1-\mu)-u2*(1-\mu)^2
#
# If you make use of this routine in your research, please cite
# Mandel & Agol (2002).
#
# Input:
#  z0  : separation between center of planet & star (in units of stellar radius)
#  u1  : linear limb-darkening coefficient
#  u2  : quadratic limb-darkening coefficient
#  p   : radius ratio of planet to star
#
# Output:
#  mu  : flux of the star
#  mu0 : flux of the star ignoring limb-darkening
#
  z = abs(z0)
# Go through the cases:
# Case 1- No planet or no obscuration:
  if (p <= 0.0) || (z >= 1.0+p) 
    return 1.0
  end
# Case 11- complete obscuration:
  if (p >= 1.0) &&  ((0.0 <= z) && (z <= (p-1.0)))
    return 0.0
  end

# Define some constants used in the lightcurve computation:
  lambdad = 0.0
  etad = 0.0
  lambdae = 0.0
  omega = 1-u1/3-u2/6
  a = (z-p)^2
  b = (z+p)^2
  q =  p^2-z^2

# Case 2: ingress/egress
  if (0.5+abs(p-0.5)) < z < (1.0+p)
    k2 = (1-a)/(4z*p)
    ek,kk = ellke(k2)
    lambdad = (((1-b)*(2b+a-3)-3q*(b-2))*kk + 4p*z*(z^2+7p^2-4)*ek-
             3q/a*ellpic_bulirsch(1.0/a-1.0,k2))/9/pi/sqrt(p*z)
    coskap0=0.5*(p^2+z^2-1)/p/z
    abs(coskap0) <= 1 ? kap0=acos(coskap0) : coskap0 < 0.0 ? kap0 = pi : kap0 = 0
    coskap1=0.5*(1-p^2+z^2)/z
    abs(coskap1) <= 1 ? kap1=acos(coskap1) : coskap1 < 0.0 ? kap1 = pi : kap1 = 0
    etad = (kap1 + p^2*(p^2+2z^2)*kap0 - 0.25*(1+5p^2+z^2)*sqrt((1-a)*(b-1)))/(2pi)
    sqt=0.25*(4z^2-(1+z^2-p^2)^2)
    lambdae = (p^2*kap0+kap1-sqrt(sqt*(sqt >= 0)))/pi
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println(z,' ',"Case 2")
  end
# Cases 3/9: planet within disk of star
  if ((p < 0.5) && (p < z < 1-p)) || ((0 < p < 1) && (0 < z < 0.5-abs(p-0.5)))
    k2 = 4z*p/(1-a)
    ek,kk = ellke(k2)
#    println("k ",k," ek ",ek," kk ",kk)
    lambdad = 2//9/pi/sqrt(1-a)*((1-5z^2+p^2+q^2)*kk + (1-a)*(z^2+7p^2-4)*ek-
             3q/a*ellpic_bulirsch(b/a-1,k2))
#    println("lambdad ",lambdad)
    etad = 0.5*p^2*(p^2+2z^2)
#    println("etad ",etad)
    lambdae = p^2
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println("lambdae ",lambdae)
#    println(z,' ',"Cases 3/9")
  end
# Case 4: 2nd/3rd contact
  if (p < 0.5) && (z == 1-p)
    lambdad = 2//3/pi*acos(1-2p)-4//9/pi*(3+2p-8p^2)*sqrt(p*(1-p))-2/3*(p > 0.5)
    etad = 0.5*p^2*(p^2+2z^2)
    lambdae = p^2
    return 1.0 - ((1.-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println("Case 4")
  end
# Case 5: edge of planet at center of star
  if (p < 0.5) && (z == p)
    k2 = 4p^2
    ek,kk = ellke(k2)
    lambdad = 1//3 + 2//9/pi*(4*(2p^2-1)*ek+(1-4p^2)*kk)
    etad = 0.5*p^2*(p^2+2z^2)
    lambdae = p^2
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println("Case 5")
  end
# Case 6: edge of planet at center of star, p=1/2
  if (p == z) && (p == 0.5)
    lambdad = 1//3 - 4//9/pi
    etad = 3//32
    lambdae = p^2
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println("Case 6")
  end
# Cases 7/8 
  if p > 0.5
#  Case 7: edge of planet at origin of star (and crossing limb of star)
    if z == p
      k2 = 0.25/p^2
      ek,kk = ellke(k2)
      lambdad = 1//3 + 16//9/pi*p*(2p^2-1)*ek-(1-4p^2)*(3-8p^2)/9/pi/p*kk
#      println("Case 7")
    end
# Case 8: ingress/egress
    if abs(1-p) <= z < p
      k2 = (1-a)/(4z*p)
      ek,kk = ellke(k2)
      lambdad = ((1-b)*(2b+a-3)-3q*(b-2))*kk + 4p*z*(z^2+7p^2-4)*ek
      lambdad -= 3q/a*ellpic_bulirsch(1.0/a-1,k2)
      lambdad *= 1//9/pi/sqrt(p*z)
#      println("Case 8")
    end
    coskap0=0.5*(p^2+z^2-1)/p/z
    abs(coskap0) <= 1 ? kap0=acos(coskap0) : coskap0 < 0.0 ? kap0 = pi : kap0 = 0
    coskap1=0.5*(1-p^2+z^2)/z
    abs(coskap1) <= 1 ? kap1=acos(coskap1) : coskap1 < 0.0 ? kap1 = pi : kap1 = 0
#    println("z ",z," p ",p," (1-a)*(b-1) ",(1-a)*(b-1))
    etad = (kap1 + p^2*(p^2+2z^2)*kap0 - 0.25*(1+5p^2+z^2)*sqrt((1-a)*(b-1)))/(2pi)
    lambdae = (p^2*kap0+kap1-sqrt(0.25*(4z^2-(1+z^2-p^2)^2)))/pi
#    println("z ",z," p ",p," lambdae ",lambdae," kap0 ",kap0," kap1 ",kap1)
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
  end
# Case 10: planet centered on star
  if (0 < p < 1) && z == 0
    lambdad = -2//3*(1-p^2)^(3//2)
    etad = 0.5*p^2*(p^2+2z^2)
    lambdae = p^2
    return 1.0 - ((1.0-u1-2u2)*lambdae + (u1+2u2)*(lambdad+2//3*(p > z)) + u2*etad)/omega
#    println("Case 10")
  end
# Should have already returned
  return Inf
end

function ellke(k2)
  return ellec_bulirsch(k2),ellk_bulirsch(k2)
end

function ellec_bulirsch(k2::Float64)
# Computes the complete elliptic integral of the second
# kind, using the approach of Bulirsch (1965)
m=1.0
a=1.0
b=1.0-k2
kc=sqrt(b)
c=a
a+=b

b=2.0*(c*kc+b)
c=a
m0=m
m+=kc
a+=b/m
while abs(m0-kc) > (1e-14*m0)
  kc=2.0*sqrt(kc*m0)
  b=2.0*(c*kc+b)
  c=a
  m0=m
  m+=kc
  a+=b/m
#  println(a,' ',b,' ',c,' ',m0,' ',kc)
end
pi*0.25*a/m
end

function ellk_bulirsch(k2::Float64)
# Computes the complete elliptic integral of the first
# kind, using the approach of Bulirsch (1965)
kc=sqrt(1.0-k2)
h=1.0
m=1.0
while abs(h-kc) > (1e-14*h)
  h=m
  m+=kc
  kc=sqrt(h*kc)
  m=0.5*m
#  println(h,' ',m,' ',kc)
end
h=m
m+=kc
pi/m
end

function ellpic_bulirsch( n::Float64, k2::Float64 )
# Computes the complete elliptical integral of the third kind using
# the algorithm of Bulirsch (1965):
# Note that n is defined with opposite sign in the ellpic_bulirsch
# than it is in Mandel & Agol (2002)

#  println(nn)
kc = sqrt(1.0-k2)
@assert(kc*(n+1.0) != 0.0)
ee = kc
m0 = 1.0
if n > -1.0
  c = 1.0
  p = sqrt(n+1.0)
  d = 1.0/p
else
  g = -n
  f = -k2-n
  p = sqrt(f/g)
  d = -k2/(g*p)
  c = 0.0
end

f = c
c = d/p + f
g = ee/p
d = (f*g+d)*2.0
p += g
g = m0
m0 += kc

while abs(g-kc) > (1e-14*g)
  kc = 2.0*sqrt(ee)
  ee = kc*m0
  f = c
  c = d/p + f
  g = ee/p
  d = (f*g+d)*2.0
  p += g
  g = m0
  m0 += kc
end
ellpic= 0.5*pi*(c*m0+d)/(m0*(m0+p))
end
