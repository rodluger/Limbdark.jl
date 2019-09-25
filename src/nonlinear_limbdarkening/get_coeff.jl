    # This routine gets the coefficients for the non-linea2
# limb-darkening law fit to a function which goes to
# zero at mu=0 and mu=1:

include("/Users/ericagol/Computer/Julia/regress.jl")

function get_coeff!(mu::Array{T,1},fnon::Array{T,1},i::Int64,nmu::Int64) where {T <: Real}
# i is the order of the fit; nmu is number of mu points which are fit.
fn = zeros(T,i-1,nmu)
for j=2:i
  fn[j-1,:] .= (j+2)*mu.^j-j*mu.^(j-2)-2mu
end
ci,cov = regress(fn,fnon,ones(BigFloat,nmu))
fmod = zeros(BigFloat,nmu)
for j=1:i-1
  fmod += ci[j]*fn[j,:]
end
return ci,fmod
end
