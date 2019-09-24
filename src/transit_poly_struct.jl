
include("define_constants.jl")
# Include linear algebara:
if VERSION >= v"0.7"
  using LinearAlgebra
end

# Include definition of Transit structure type:
include("transit_structure.jl")
# Include code which computes linear limb-darkening term:
include("s2.jl")
# Include code which computes M_n and N_n:
include("Mn_compute.jl")
include("Nn_compute.jl")

"""
    compute_uniform(t)

Computes the coefficient for the uniform disk case, `s_0`, given
a `TransitStruct` instance `t`.
"""
function compute_uniform!(t::Transit_Struct{T}) where {T <: Real}
  r=t.r; b=t.b; r2=r*r; b2=b*b
  # Uniform disk case:
  # Compute sn[1] and its derivatives:
  if b <= 1-r  # k^2 > 1
    t.sn[1] = pi*(1-r2)
    if t.grad
      t.dsndr[1] = -2*pi*r
      t.dsndb[1] = 0.
    end
    t.kap0 = convert(T,pi); t.kck = zero(T)
  else
    # Twice area of kite-shaped region connecting centers of circles & intersection points:
  #  t.kite_area2 = sqrt(sqarea_triangle(one(T),b,r)) 
    t.kite_area2 = sqrt(t.sqarea)
    # Angle of section for occultor:
    t.kap0  = atan(t.kite_area2,(r-1)*(r+1)+b2)
    # Angle of section for source:
    t.pimkap1 = atan(t.kite_area2,(r-1)*(r+1)-b2)
    # Flux of visible uniform disk:
    t.sn[1] = t.pimkap1 - r2*t.kap0 + .5*t.kite_area2
    if t.grad
      t.dsndr[1]= -2*r*t.kap0
      t.dsndb[1]= t.kite_area2/b
    end
    t.kck = t.kite_area2*t.fourbrinv
  end
  return
end

"""
    sqarea_triangle(a,b,c)

Function which computes sixteen times the square of the area
of a triangle with sides a, b and c using Kahan method.
How to compute (sixteen times the) area squared of triangle with 
high precision (Goldberg 1991).
"""
function sqarea_triangle(a::T,b::T,c::T) where {T <: Real}

  # First, do a quick sort of three numbers:
  if c > b
    tmp = c; c = b; b = tmp
  end
  if c > a
    tmp = a; a = c; c = tmp
  end
  if b > a
    tmp = b; b = a; a = tmp
  end
  area = (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
  return area

end

"""
    transit_poly_g(t)

Given a `TransitStruct` instance `t`, computes a limb-darkened transit light 
curve without gradient.


"""
function transit_poly_g(t::Transit_Struct{T}) where {T <: Real}
  # Number of limb-darkening components to include (beyond 0 and 1):
  # We are parameterizing these with the function:
  # g_n = g_n [(n+2) z^n - n z^{n-2}] for n >= 2
  # while g_{0} = g_0 z^0 (uniform source) and g_{1} = g_1 z^1 (linear limb-darkening)
  # which gives a Green's integral of:
  # P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
  # for which we have a solution in terms of M_n.
  # The variable "t" is a structure which contains transit parameters
  # and intermediate quantities computed from these:
  r=t.r; b=t.b; n = t.n; r2=r*r; b2 = b*b
  # Set up a vector for storing results of P(G_n)-Q(G_n); note that
  # this is a different vector than the Starry case:

  # Check for different cases:
  if b >= 1+r
    # unobscured - return one:
    return one(T)
  end
  if r >= 1+b
    # full obscuration - return zero:
    return zero(T)
  end
  if b == 0.0
    # Annular eclipse - integrate around the full boundary of both bodies:
    onemr2 = 1-r2
    flux = zero(T); t.sqrt1mr2 = sqrt(onemr2)
    flux = (t.g_n[1]*onemr2+t.twothird*t.g_n[2]*t.sqrt1mr2^3)
    fac= 2r2*onemr2
    @inbounds for i=2:t.n
      flux += -t.g_n[i+1]*fac
      fac *= t.sqrt1mr2
    end
    return flux*pi*t.den
  else
  # Next, compute k^2 = m:
    t.onembmr2=(r+1-b)*(1-r+b); t.fourbr = 4b*r; t.fourbrinv = inv(t.fourbr)
    t.sqbr = sqrt(b*r); t.sqbrinv = inv(t.sqbr)
    t.onembmr2inv=inv(t.onembmr2); t.sqonembmr2 = sqrt(t.onembmr2)
    t.onembpr2 = (1-r-b)*(1+b+r)
    t.sqarea = sqarea_triangle(one(T),r,b)
  #  t.k2 = t.onembmr2*t.fourbrinv
    t.k2 = t.onembpr2*t.fourbrinv+1
    t.onemr2mb2 = (1.0-r)*(1.0+r)-b2
    if t.k2 > 0
      t.k = sqrt(t.k2)
    else
      println("negative k2: ",t.k2," r: ",r," b: ",b)
      t.k2 = 0.0; t.k = 0.0
    end
    if t.k2 > 1
      if t.k2 > 2.0
        t.kc2 = 1.0-inv(t.k2)
        t.kc = sqrt(t.kc2)
      else
        t.kc2 = t.onembpr2/((1-b+r)*(1-r+b))
        t.kc = sqrt(t.kc2)
      end
    else
      if t.k2 > 0.5
        t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
        t.kc = sqrt(t.kc2)
      else
        t.kc2 = 1.0-t.k2
        t.kc = sqrt(t.kc2)
      end
    end
  end

  # Compute uniform case, sn[1]:
  compute_uniform!(t)
  # If uniform model, then return:
  if t.n == 0
    return t.g_n[1]*t.sn[1]*t.den
  end

  # Compute linear case, sn[2]:
  #s2!(t)
  t.sn[2],t.Eofk,t.Em1mKdm = s2_ell(r,b)
  #sn2,Eofk,Em1mKdm = s2_ell(r,b)
  # Now, compare results:
  #if abs(t.sn[2]-sn2) > 1e-8*abs(sn2)
  #  println("r: ",r," b: ",b," sn[2]: ",t.sn[2]," sn2: ",sn2)
  #end
  # Special case of linear limb-darkening:
  if t.n == 1
    flux = t.g_n[1]*t.sn[1]+t.g_n[2]*t.sn[2]
    flux *= t.den  # for g_2 and above, the flux is zero.
    return flux
  end

  if t.n >= 2
  # Special case of quadratic limb-darkening:
  # Transformed expressions from Mandel & Agol:
    eta2 = r2*(r2+2*b2)
    if t.k2 >= 1
      four_pi_eta = 2pi*(eta2-1.0)
    else
      four_pi_eta = 2*(-t.pimkap1+eta2*t.kap0-0.25*t.kite_area2*(1.0+5r2+b2))
    end
    t.sn[3] = 2*t.sn[1]+four_pi_eta
    if t.n == 2
      flux = t.g_n[1]*t.sn[1]+t.g_n[2]*t.sn[2]+t.g_n[3]*t.sn[3]
      flux *= t.den  # for g_2 and above, the flux is zero.
      return flux
    end
  end

  # Compute the M_n functions:
  if t.k2 > 0
  #  if (t.k2 < 0.5 || t.k2 > 20.0) && t.n > 3
    if (t.k2 < 0.5) && t.n > 3
  # This computes M_n for the largest four m, and then works down to smaller values:
      Mn_lower!(t)
    else
  # This computes Mn for m=0 to 3, and then works upward to larger m:
      Mn_raise!(t)
    end
  end

  # Add up first three terms in flux numerator:
  flux = t.g_n[1]*t.sn[1]+t.g_n[2]*t.sn[2]+t.g_n[3]*t.sn[3]
  # Next, loop over the Green's function components:
  @inbounds for n=3:t.n
  #  pofgn_M = (1+(r-b)*(r+b))*t.Mn[n+1]-t.Mn[n+3]
  #  pofgn_M = t.onemr2mb2*t.Mn[n+1]-t.Mn[n+3]
  #  pofgn_M = 2*r^2*t.Mn[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mn[n+1]+(1-(b-r)^2)*((b+r)^2-1)*t.Mn[n-1])
  #  pofgn_M = 2*r^2*t.Mn[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mn[n+1]+sqarea_triangle(one(T),r,b)*t.Mn[n-1])
  #  pofgn_M = 2*r^2*t.Mn[n+1]-n/(n+2)*((1-r^2-b^2)*t.Mn[n+1]+t.sqarea*t.Mn[n-1])
  #  pofgn_M = 2*r2*t.Mn[n+1]-n/(n+2)*(t.onemr2mb2*t.Mn[n+1]+t.sqarea*t.Mn[n-1])
    pofgn_M = 2*r2*t.Mn[n+1]-n*t.ninv[n+2]*(t.onemr2mb2*t.Mn[n+1]+t.sqarea*t.Mn[n-1])
  # Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
  # boundary for n > 0.
  # Compute sn[n]:
    t.sn[n+1] = -pofgn_M
      flux += t.g_n[n+1]*t.sn[n+1]
  end
  # That's it!
  # flux = sum(t.g_n.*t.sn)/(pi*(t.g_n[1]+t.twothird*t.g_n[2]))  # for g_2 and above, the flux is zero.
  # Divide by denominator:
  flux *= t.den  # for g_2 and above, the flux is zero.
  return flux
end


"""
    transit_poly!(r,b,u_n,dfdrbu)

Given a radius-ratio, impact parameter, vector of limb-darkening coefficients
of size N, and pre-allocated vector for derivatives of size N+2, returns
the flux (normalized to one for unocculted star) for a limb-darkened transit 
light curve and its gradient with respect to the input parameters is returned.

# Arguments
- `r::Real`: The radius of the occultor in units of the radius of the occulted body.
- `b::Real`: The (initial) impact parameter of the occultation.
- `u_n::Array{Real,1}`: The array of limb darkening coefficients.  Size N.
- `dfdrbu::Array{Real,1}`: Gradient of the flux with respect to all input parameters. Size is N+2
"""
function transit_poly!(r::T,b::T,u_n::Array{T,1},dfdrbu::Array{T,1}) where {T <: Real}
  # Initialize `TransitStruct`:
  t = transit_init(r,b,u_n,true)
  # Pass g_n (without last two dummy values):
  flux = transit_poly_g!(t)
  # Now, transform derivaties from g to u:
  fill!(dfdrbu,zero(T))
  # u_n derivatives:
  dfdrbu[1] = t.dfdrb[1]
  dfdrbu[2] = t.dfdrb[2]
  @inbounds for i=1:t.n, j=0:t.n
    dfdrbu[i+2] += t.dfdg[j+1]*t.dgdu[j+1,i]
  end
  return flux
end

"""
    transit_poly(r,b,u_n)

Given a radius-ratio, impact parameter, vector of limb-darkening coefficients
of size N, returns the flux (normalized to one for unocculted star) for a 
limb-darkened transit.

# Arguments
- `r::Real`: The radius of the occultor in units of the radius of the occulted body.
- `b::Real`: The (initial) impact parameter of the occultation.
- `u_n::Array{Real,1}`: The array of limb darkening coefficients.  Size N.
"""
function transit_poly(r::T,b::T,u_n::Array{T,1}) where {T <: Real}
  t = transit_init(r,b,u_n,false)
  # Pass g_n (without last two dummy values):
  return transit_poly_g(t) 
end

"""
    transit_poly!(t)

Given a `TransitStruct` instance `t`, with Green's coefficients already initialized,
computes a limb-darkened transit light curve with optional gradient with
respect to r, b and u coefficients.

"""
function transit_poly!(t::Transit_Struct{T}) where {T <: Real}
  # This function assumes that g_n has already been computed from u_n
  # (this can be used to save compute time when limb-darkening is fixed for
  # a range of radii/impact parameters).
  # Pass transit structure, and compute flux:
  if t.grad
    flux = transit_poly_g!(t)
    # Now, transform derivaties from g to u:
    fill!(t.dfdu,zero(T))
    BLAS.gemv!('T',1.0,t.dgdu,t.dfdg,0.0,t.dfdu)
  #  @inbounds for i=1:t.n, j=0:t.n
  #    t.dfdu[i] += t.dfdg[j+1]*t.dgdu[j+1,i]
  #  end
    return flux
  else
    return transit_poly_g(t)
  end
  return
end

"""
    transit_poly_g!(t)

Given a `TransitStruct` instance `t`, with Green's coefficients already initialized,
computes a limb-darkened transit light curve including gradient with
respect to r, b and g coefficients.

"""
function transit_poly_g!(t::Transit_Struct{T}) where {T <: Real}
  r = t.r; b=t.b; n = t.n; r2 =r*r; b2=b*b; bcut = 1e-3
  @assert((length(t.g_n)) == length(t.dfdg))
  @assert(r > 0)
  # Number of limb-darkening components to include (beyond 0 and 1):
  # We are parameterizing these with the function:
  # g_n = g_n [(n+2) z^n - n z^{n-2}] for n >= 2
  # while g_{0} = g_0 z^0 (uniform source) and g_{1} = g_1 z^1 (linear limb-darkening)
  # which gives a Green's integral of:
  # P(G_n) = \int_{\pi-\phi}^{2\pi+\phi} (1-r^2-b^2-2br s_\varphi)^{n/2} (r+b s_\varphi) d\varphi
  # for which we have a solution in terms of M_n.
  # Compute the derivative of the flux with respect to the different coefficients.

  # Set up a vector for storing results of P(G_n)-Q(G_n); note that
  # this is a different vector than the Starry case:
  # Check for different cases:
  if b >= 1+r || r ==  0.0
    # unobscured - return one, and zero derivatives:
    fill!(t.dfdrb,zero(T))
    fill!(t.dfdg,zero(T))
    return one(T)
  end
  if r >= 1+b
    # full obscuration - return zero, and zero derivatives:
    fill!(t.dfdrb,zero(T))
    fill!(t.dfdg,zero(T))
    return zero(T)
  end
  if b == 0.0
    # Annular eclipse - integrate around the full boundary of both bodies:
    flux = zero(T); onemr2 = 1-r2; t.sqrt1mr2 = sqrt(onemr2)
    fill!(t.dfdg,zero(T))
    flux = (t.g_n[1]*onemr2+t.twothird*t.g_n[2]*t.sqrt1mr2^3)*pi*t.den
    fac  = 2r2*onemr2*pi*t.den
    facd = -2r*pi*t.den
    t.dfdrb[1] = t.g_n[1]*facd + t.g_n[2]*facd*t.sqrt1mr2
    @inbounds for i=2:t.n
      flux -= t.g_n[i+1]*fac
      t.dfdrb[1] += t.g_n[i+1]*facd*(2*onemr2-i*r2)
      t.dfdg[i+1] -= fac
      fac *= t.sqrt1mr2
      facd *= t.sqrt1mr2
    end
    #  dfdrb[2]=0 since the derivative with respect to b is zero.
    t.dfdg[1] = (onemr2-flux)*pi*t.den
    t.dfdg[2] = t.twothird*(t.sqrt1mr2^3-flux)*pi*t.den
    return flux
  else
  # Next, compute k^2 = m:
    t.onembmr2=(r-b+1)*(1-r+b); t.fourbr = 4b*r; t.fourbrinv = inv(t.fourbr)
    t.sqbr = sqrt(b*r); t.sqbrinv = inv(t.sqbr)
    t.onembmr2inv=inv(t.onembmr2); t.sqonembmr2 = sqrt(t.onembmr2)
    t.onembpr2 = (1-r-b)*(1+b+r)
    t.sqarea = sqarea_triangle(one(T),r,b)
    t.k2 = t.onembmr2*t.fourbrinv; 
    t.onemr2mb2 = (1.0-r)*(1.0+r)-b2
    if t.k2 > 0
      t.k = sqrt(t.k2)
    else
      println("negative k2: ",t.k2," r: ",r," b: ",b)
      t.k2 = 0.0; t.k = 0.0
    end
    if t.k2 >= 1
      if t.k2 > 2.0
        t.kc2 = 1.0-inv(t.k2)
        t.kc = sqrt(abs(t.kc2))
      else
        t.kc2 = t.onembpr2/((1+r-b)*(1-r+b))
        t.kc = sqrt(abs(t.kc2))
      end
    else
      if t.k2 > 0.5
        t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
        t.kc = sqrt(abs(t.kc2))
      else
        t.kc2 = (r-1+b)*(b+r+1)*t.fourbrinv
        t.kc = sqrt(abs(t.kc2))
      end
    end
  end

  # Compute uniform case:
  compute_uniform!(t)

  # Compute sn[2] and its derivatives:
  #if t.n >= 1
  #  s2!(t)
  #  # debug:
  #  s2_grad = zeros(T,2)
  #  sn2,Eofk,Em1mKdm = s2!(r,b,s2_grad)
  #  if (abs(t.sn[2]-sn2) > 1e-8) || (abs(t.s2_grad[1]-s2_grad[1]) > 1e-8) || (abs(t.s2_grad[2]-s2_grad[2]) > 1e-8)
  #    println("r: ",r,"b: ",b," s2 new: ",t.sn[2]," s2 old: ",sn2)
  #    println("ds2/dr new: ",t.s2_grad[1]," ds2/dr old: ",s2_grad[1])
  #    println("ds2/db new: ",t.s2_grad[2]," ds2/db old: ",s2_grad[2])
  #    read(STDIN,Char)
  #  end
  #end
  if t.n >= 1
    t.sn[2],t.Eofk,t.Em1mKdm = s2!(r,b,t.s2_grad)
    t.dsndr[2] = t.s2_grad[1]
    t.dsndb[2] = t.s2_grad[2]
  end


  # Special case of quadratic limb-darkening:
  if t.n >= 2
  #  if t.n == 2
  # Transformed expressions from Mandel & Agol for n=2:
    r2pb2 = r2+b2
    eta2 = r2*(r2pb2+b2)
    deta2dr =  2*r*r2pb2
    deta2db = 2*b*r2
    if t.k2 > 1
      four_pi_eta = 2pi*(eta2-1.0)
      detadr = 4pi*deta2dr
      detadb = 4pi*deta2db
    else
      four_pi_eta = 2*(-t.pimkap1+eta2*t.kap0-0.25*t.kite_area2*(1.0+5r2+b2))
      detadr = 8r*(r2pb2*t.kap0-t.kite_area2)
      detadb = 2/b*(4*b2*r2*t.kap0-(1+r2pb2)*t.kite_area2)
    end
    t.sn[3] = 2*t.sn[1]+four_pi_eta
    t.dsndr[3] = 2*t.dsndr[1]+detadr
    t.dsndb[3] = 2*t.dsndb[1]+detadb
  #  else
    if t.n > 2
      # Compute the M_n functions:
      if t.k2 > 0
  #      if (t.k2 < 0.5 || t.k2 > 20.0) && t.n > 3
        if (t.k2 < 0.5) && t.n > 3
      # This computes Mn for largest four m, and then works down to smaller values:
          Mn_lower!(t)
          if b < bcut
            Nn_lower!(t)
          end
        else
      # This computes Mn and then works upward to larger m:
          Mn_raise!(t)
          if b < bcut
            Nn_raise!(t)
          end
        end
      end
      # Next, loop over the Green's function components:
      binv = inv(b)
  #    @inbounds for n=2:t.n
      @inbounds for n=3:t.n
  #      pofgn_M = (1+(r-b)*(r+b))*t.Mn[n+1]-t.Mn[n+3]
  #      pofgn_M = t.onemr2mb2*t.Mn[n+1]-t.Mn[n+3]
  #      pofgn_M = 2*r^2*t.Mn[n+1]-n/(n+2)*((1-r2-b2)*t.Mn[n+1]+t.kite_area2^2*t.Mn[n-1])
  #      pofgn_M = 2*r^2*t.Mn[n+1]-n/(n+2)*(t.onemr2mb2*t.Mn[n+1]+sqarea_triangle(one(T),r,b)*t.Mn[n-1])
  #      pofgn_M = 2*r2*t.Mn[n+1]-n/(n+2)*(t.onemr2mb2*t.Mn[n+1]+t.sqarea*t.Mn[n-1])
        pofgn_M = 2*r2*t.Mn[n+1]-n*t.ninv[n+2]*(t.onemr2mb2*t.Mn[n+1]+t.sqarea*t.Mn[n-1])
      # Q(G_n) is zero in this case since on limb of star z^n = 0 at the stellar
      # boundary for n > 0.
      # Compute sn[n]:
        t.sn[n+1] = -pofgn_M
        dpdr_M = 2*r*((n+2)*t.Mn[n+1]-n*t.Mn[n-1])
  #      dpdr_M = 2*r*((2-(n+2)*(b-r)^2)*t.Mn[n-1]-4*(n+2)*b*r*t.Nn[n-1])
  #      dpdr_M = 2*r*((2+n-n/t.onembmr2)*t.Mn[n+1]-n/t.k2*t.Nn[n+1])
        t.dsndr[n+1] = -dpdr_M
        if b < bcut
        # For small b, use function which doesn't involve division by b:
          dpdb_M = n*(t.Mn[n-1]*(2*r^3+b^3-b-3*r2*b)+b*t.Mn[n+1]-4*r^3*t.Nn[n-1])
        else
          dpdb_M = n*binv*((t.Mn[n+1]-t.Mn[n-1])*(r2+b2)+(b2-r2)^2*t.Mn[n-1])
        end
        t.dsndb[n+1] = -dpdb_M
      end
    end
  end
  # That's it!
  # Compute derivatives with respect to the coefficients:
  flux = zero(T)
  t.dfdrb[1]=zero(T)  # Derivative with respect to r
  t.dfdrb[2]=zero(T)  # Derivative with respect to b
  @inbounds for i=0:t.n
    # derivatives with respect to the coefficients:
    t.dfdg[i+1]= t.sn[i+1]*t.den
    # total flux:
    flux += t.g_n[i+1]*t.dfdg[i+1]
    # derivatives with respect to r and b:
    t.dfdrb[1] += t.g_n[i+1]*t.dsndr[i+1]*t.den
    t.dfdrb[2] += t.g_n[i+1]*t.dsndb[i+1]*t.den
  end
  # Include derivatives with respect to first two g_n parameters:
  t.dfdg[1] -= flux*t.den*pi
  t.dfdg[2] -= flux*t.den*pi*t.twothird
  return flux
end
